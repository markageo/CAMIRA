#include "ImmersedBoundary.h"
#include "../Tools/FVLookups.h"
#include "../Tools/FVTools.h"

#include "../Macros.h"
#include "../IO/ArrayIO.h"
#include "../IO/IOTools.h"
#include "../Tools/SweepTransformations.h"

#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <cmath>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/squared_distance_3.h> 
#include <CGAL/boost/graph/IO/STL.h>


namespace CFD
{

namespace 
{

// Create a mask array for cell centers
Tensor3D CreateCellMask( const Tree &tree, 
                         const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    Tensor3D mask = Tensor3D( mesh.nCells[X] + 2*CFD::nGhost, mesh.nCells[Y] + 2*CFD::nGhost, mesh.nCells[Z] + 2*CFD::nGhost ).setConstant( CellType::Fluid );

    // Identify cells with cell centeres inside the solid
    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                floatType xq = mesh.cellCenters[X](i),
                          yq = mesh.cellCenters[Y](j),
                          zq = mesh.cellCenters[Z](k);

                if ( PointInside( tree, xq, yq, zq ) )
                    mask(G(i, j, k)) = CellType::Solid;

            }
        }
    }


    return mask;
}



// Check if a cell is a fluid cell within the domain boundary
bool CellIsFluid( const TensorIndex3D &cellIndex_a,
                  const Tensor3D &mask,
                  const Mesh &mesh )
{
    using FVT::G;

    if ( static_cast<intType>( mask( G(cellIndex_a) ) ) == CellType::Solid ) 
        return false;

    for ( intType a = 0; a != Axis::count; a++ ) {
        Axis::ENUMDATA axis = static_cast<Axis::ENUMDATA>( a );

         if ( cellIndex_a[axis] < 0 )
            return false;

        if ( cellIndex_a[axis] > ( mesh.nCells[axis] - 1 ) )
            return false;

    }

    return true;
}



// Sets data for a particular source term in a particular direciton for a particular cell
void AddIBDataForDirection( IBCell &ibCell, 
                            const Axis::ENUMDATA axis,
                            const intType directionIndex,
                            const Tensor3D &mask,
                            const Mesh &mesh,
                            const Tree &tree, 
                            const InputData &inputData )
{
    using enum Axis::ENUMDATA;

    ibCell.sourceTermsData.emplace_back();
    IBCell::SourceTermData &sourceTermData = ibCell.sourceTermsData.back();

    TensorIndex3D &cellIndex = ibCell.cellIndex;

    TensorIndex3D cellIndex_g    = cellIndex;
    cellIndex_g[axis] += directionIndex;

    TensorIndex3D cellIndex_a = cellIndex;
    cellIndex_a[axis] -= directionIndex;

    sourceTermData.direction          = axis;
    sourceTermData.directionIndex     = directionIndex;
    sourceTermData.faceDirectionIndex = ( directionIndex == 1 ) ? 1 : 0 ;
    sourceTermData.cellIndex_g        = cellIndex_g;
    sourceTermData.cellIndex_a        = cellIndex_a;

    intType fidx = cellIndex[axis] + sourceTermData.faceDirectionIndex;

    // Ensure that there is enough space between the IB and other IBs and the domain boundary
    if ( !CellIsFluid( cellIndex_a, mask, mesh ) ) {
        throw std::runtime_error( "Invalid immersed boundary geometry and mesh specification: There must be at least one fluid cell between domain and solid boundaries on all grid levels!" );
    }

    switch ( inputData.geoemtryBoundaryTreatement ) {
        case GeometryBoundaryTreatement::DirectionalImmersedBoundary:
        {
            // Distance from cell center to immersed boundary along this coordinate direction
            fVector3 queryPointCoords( mesh.cellCenters[X](cellIndex[X]),
                                       mesh.cellCenters[Y](cellIndex[Y]),
                                       mesh.cellCenters[Z](cellIndex[Z]) );
            fVector3 rayDirection( 0, 0, 0 );
            rayDirection[ axis ] = static_cast<floatType>( directionIndex );
            sourceTermData.ibDistance = NearestRayIntersection(tree, queryPointCoords, rayDirection);
            break;
        }
        
        case CFD::GeometryBoundaryTreatement::Staircase:
        {   
            // Distance to the nearest cell face, approximates the immersed boundary to be on the cell face
            // i.e. staircase approximation
            sourceTermData.ibDistance = mesh.cellLengths[axis](cellIndex[axis]) / 2.0f;
            break;
        }

        default:
            break;
    }
    const floatType &ibDistance = sourceTermData.ibDistance;


    // Face area vector
    sourceTermData.faceAreaComponent = static_cast<floatType>( directionIndex )   // Gives the correct sign
                                     * mesh.cellLengths[LUT::LoOrthogonalAxis[axis]]( cellIndex[ LUT::LoOrthogonalAxis[axis] ] )
                                     * mesh.cellLengths[LUT::HiOrthogonalAxis[axis]]( cellIndex[ LUT::HiOrthogonalAxis[axis] ] );
                                    //  * mesh.cellFaceAreas[axis]( cellIndex[ LUT::LoOrthogonalAxis[axis] ], cellIndex[ LUT::HiOrthogonalAxis[axis] ] );    
    


    // Extrapolation coefficients from fluid onto face. May use further points due to stability condition
    bool meetsGhostStabilityCondition =  ibDistance >= ( mesh.cellLengths[axis](cellIndex[axis]) / 2.0f );
    if ( meetsGhostStabilityCondition ) {

        floatType dxpOn2 = mesh.cellLengths[axis](cellIndex[axis]) / 2.0f ;
        sourceTermData.faceExtrapCoeff_p  = ( ibDistance - dxpOn2 ) / ( ibDistance );
        sourceTermData.faceExtrapCoeff_a  = 0.0f;
        sourceTermData.faceExtrapCoeff_ib = dxpOn2 / ibDistance;

    } else {

        floatType dxp = mesh.cellLengths[axis](cellIndex[axis]);
        floatType dxa = mesh.cellLengths[axis](cellIndex_a[axis]);
        floatType denominator =  dxp / 2.0f   +  dxa / 2.0f  +  ibDistance;
        sourceTermData.faceExtrapCoeff_p  = 0.0f;
        sourceTermData.faceExtrapCoeff_a  = ( dxp / 2.0f - ibDistance ) / denominator;
        sourceTermData.faceExtrapCoeff_ib = ( dxa / 2.0f + dxp ) / denominator;

    }


    // Extrapolation coefficients from face to ghost cell
    if        ( directionIndex == +1 ) {    // Face on Hi side

        sourceTermData.ghostExtrapCoeff_p = - (1 - mesh.interpFactors[axis](fidx)) / mesh.interpFactors[axis](fidx);
        sourceTermData.ghostExtrapCoeff_f = 1.0f / mesh.interpFactors[axis](fidx);

    } else if ( directionIndex == -1 ) {   // Face on Lo side

        sourceTermData.ghostExtrapCoeff_p = - mesh.interpFactors[axis](fidx) / (1.0f - mesh.interpFactors[axis](fidx));
        sourceTermData.ghostExtrapCoeff_f = 1.0f / (1.0f - mesh.interpFactors[axis](fidx));

    }


    // Extrapolation coefficients from fluid to immersed boundary surface
    floatType cellInteriorDistance = abs( mesh.cellCenters[axis](cellIndex_a[axis]) - mesh.cellCenters[axis](cellIndex[axis]) ); 
    sourceTermData.ibExtrapFactor_p = ( cellInteriorDistance + ibDistance ) / cellInteriorDistance;
    sourceTermData.ibExtrapFactor_a = - ibDistance / cellInteriorDistance;


    // Far pressure ghost cell coefficients
    if        ( directionIndex == +1 ) {    // Ghost cell on Hi side

        floatType dxp      = mesh.cellLengths[axis](cellIndex[axis]),
                  dxe      = mesh.cellLengths[axis](cellIndex_g[axis]),
                  lambdaw  = mesh.interpFactors[axis](cellIndex[axis]    ),
                  lambdae  = mesh.interpFactors[axis](cellIndex[axis] + 1),
                  lambdaee = mesh.interpFactors[axis](cellIndex[axis] + 2),
                  le       = 1.0f /  mesh.cellCenterDiffInv[axis](cellIndex[axis] + 1);

        sourceTermData.farPressureCoeff_p = - dxe / (lambdae * lambdaee * le)
                                            - (dxe / dxp) * ( (1-lambdae) / lambdae ) * ( (1 - lambdae - lambdaw) / lambdaee )
                                            + (1 - lambdae) / lambdaee;

        sourceTermData.farPressureCoeff_a = (dxe / dxp) * ( (1-lambdae) / lambdae ) * ( (1 - lambdaw) / lambdaee );

        sourceTermData.farPressureCoeff_g =   dxe / (lambdae * lambdaee * le)
                                            - (dxe / dxp) * ( (1-lambdae) / lambdaee)
                                            - (1 - lambdae - lambdaee) / lambdaee;

    } else if ( directionIndex == -1) {     // Ghost cell on Lo side

        floatType dxp      = mesh.cellLengths[axis](cellIndex[axis]),
                  dxw      = mesh.cellLengths[axis](cellIndex_g[axis]),
                  lambdae  = mesh.interpFactors[axis](cellIndex[axis] + 1),
                  lambdaw  = mesh.interpFactors[axis](cellIndex[axis]    ),
                  lambdaww = mesh.interpFactors[axis](cellIndex[axis] - 1),
                  lw       = 1.0f / mesh.cellCenterDiffInv[axis](cellIndex[axis]);

        sourceTermData.farPressureCoeff_p = - dxw / ( (1 - lambdaw) * (1 - lambdaww) * lw)
                                            + (dxw / dxp) * ( lambdaw / (1 - lambdaw) ) * ( (1 - lambdae - lambdaw) / (1 - lambdaww) )
                                            + lambdaw / (1 - lambdaww);

        sourceTermData.farPressureCoeff_a = (dxw / dxp) * ( lambdaw / (1 - lambdaw) ) * ( lambdae / (1 - lambdaww) );

        sourceTermData.farPressureCoeff_g =   dxw / ( (1 - lambdaw) * (1 - lambdaww) * lw)
                                            - (dxw / dxp) * lambdaw / (1 - lambdaww)
                                            + (1 - lambdaw - lambdaww) / (1 - lambdaww);

    }

}



void SetVelocityFluxCorrectionCoefficient( IBData &ibData,
                                           const Mesh &mesh )
{

    // The denominator must be calculated seperately first
    floatType denominator = 0.0f;
    for ( auto &ibCell : ibData.ibCells ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            floatType ibCellFaceDistance = abs( sourceTermData.ibDistance - mesh.cellLengths[sourceTermData.direction](ibCell.cellIndex[sourceTermData.direction]) );

            denominator += std::pow( abs( sourceTermData.faceAreaComponent * ibCellFaceDistance) , 2.0f );

        }
    }

    // Now the coefficient for each cell
    for ( auto &ibCell : ibData.ibCells ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            floatType ibCellFaceDistance = abs( sourceTermData.ibDistance - mesh.cellLengths[sourceTermData.direction](ibCell.cellIndex[sourceTermData.direction]) );

            sourceTermData.velocityFluxCorrectionCoeff = ( std::pow( ibCellFaceDistance, 2.0f ) 
                                                        * sourceTermData.faceAreaComponent 
                                                        ) / denominator;

        }
    }

}



IBData ConstructIBData( const Polyhedron &geometry,
                        const Mesh &mesh,
                        const InputData &inputData )
{

    using enum Axis::ENUMDATA;
    using FVT::G;

    IBData ibData;

    // Create AABB tree for geometry
    Tree tree = MakeAABBTree( geometry );

    ibData.mask = CreateCellMask( tree, mesh );
    auto &mask = ibData.mask;

    // Iterate through all cells 
    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                if ( static_cast<intType>( mask(G(i, j, k)) ) == CellType::Solid )
                    continue;

                TensorIndex3D cellIndex{i, j, k};
                IBCell *ibCellPtr = nullptr;

                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                    auto CheckIBCellPtr = [&] () {
                        if ( ibCellPtr == nullptr ) {
                            ibData.ibCells.emplace_back();
                            ibCellPtr = &ibData.ibCells.back();
                            ibCellPtr->cellIndex = cellIndex;
                        } 
                    };

                    // Solid on hi side
                    TensorIndex3D hiSideCellIndex = cellIndex;
                    hiSideCellIndex[axis] += 1;
                    bool atHiBoundary = ( cellIndex[axis] == mesh.nCells[axis]-1  );
                    if ( !atHiBoundary && static_cast<intType>( mask(G(hiSideCellIndex)) ) == CellType::Solid ) {
                        CheckIBCellPtr();
                        AddIBDataForDirection( *ibCellPtr, axis, +1, mask, mesh, tree, inputData );
                    }

                    // Solid on lo side
                    TensorIndex3D loSideCellIndex = cellIndex;
                    loSideCellIndex[axis] -= 1;
                    bool atLoBoundary = ( cellIndex[axis] == 0  );
                    if ( !atLoBoundary && static_cast<intType>( mask(G(loSideCellIndex)) ) == CellType::Solid ) {
                        CheckIBCellPtr();
                        AddIBDataForDirection( *ibCellPtr, axis, -1, mask, mesh, tree, inputData );
                    }

                } );

            }
        }
    }

    // Calculate the coefficients for the velocity flux correction
    SetVelocityFluxCorrectionCoefficient( ibData, mesh );

    return ibData;
}


}   // end anonymous namespace




IBData CreateImmersedBoundaryData( const InputData &inputData, 
                                   const Mesh &mesh,
                                   const AxisTransformationMap &axisTransformation )
{
    using enum Axis::ENUMDATA;

    IBData ibData;

    if ( !inputData.hasIBGeometry ) {
        ibData.mask = Tensor3D( mesh.nCells[X] + 2*CFD::nGhost, mesh.nCells[Y] + 2*CFD::nGhost, mesh.nCells[Z] + 2*CFD::nGhost ).setConstant( CellType::Fluid );
        return ibData;
    }

    Polyhedron P = MakeGeometry( inputData );
    ibData = ConstructIBData( P, mesh, inputData );

    // Make a new new inputData object in the user coordinates to make another geometry to be output to file
    if ( inputData.outputGeometry ) {
        InputData inputDataUserCoordinates( inputData );
        TransformUserInputData( inputDataUserCoordinates, axisTransformation.Inverse() );
        Polyhedron PUserCoordinates = MakeGeometry( inputDataUserCoordinates );

        std::ofstream out(  IOTOOLS::RemoveFileExtension( inputData.geometryOutputFilename, ".stl" ) + ".stl" );
        CGAL::IO::write_STL( out, PUserCoordinates );
    }
    

    return ibData;
}


}   // end namespace CFD