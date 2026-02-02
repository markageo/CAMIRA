#include "ImmersedBoundary.h"
#include "../Core/FVLookups.h"
#include "../Core/FVTools.h"

#include "../Core/Macros.h"
#include "../IO/ArrayIO.h"
#include "../IO/IOTools.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"

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


namespace CAMIRA
{

namespace 
{

// Create a mask array for cell centers
Tensor3D CreateCellMask( const Tree &tree, 
                         const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    Tensor3D mask = Tensor3D( mesh.nCells[X] + 2*CAMIRA::nGhost, mesh.nCells[Y] + 2*CAMIRA::nGhost, mesh.nCells[Z] + 2*CAMIRA::nGhost ).setConstant( CellType::Fluid );

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



// Masks cells that form a single cell "cavity". These can cause instability and are a highly degenerate case
void RemoveMaskSingleCellCavities( Tensor3D &mask,
                                   const Mesh &mesh )
{
    using FVT::G;
    using enum Axis::ENUMDATA;

    // It is possible for the filling of one cavity to create another
    // Iterate until there are no more left. 
    const intType maxIterations = 25;
    for ( intType iters = 0; iters != maxIterations; iters++ ) {

        intType numCavities = 0;

        for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
            for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
                for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                    // Must be a fluid cell
                    if ( static_cast<intType>( mask(G(i, j, k)) ) == CellType::Solid )
                        continue;

                    EnumVector<Axis, intType> maskedCount{0, 0, 0};

                    // x 
                    maskedCount[X] += 1 - static_cast<intType>( mask( G( i + 1, j, k ) ) );
                    maskedCount[X] += 1 - static_cast<intType>( mask( G( i - 1, j, k ) ) );

                    // y
                    maskedCount[Y] += 1 - static_cast<intType>( mask( G( i, j + 1, k ) ) );
                    maskedCount[Y] += 1 - static_cast<intType>( mask( G( i, j - 1, k ) ) );

                    // z
                    maskedCount[Z] += 1 - static_cast<intType>( mask( G( i, j, k + 1 ) ) );
                    maskedCount[Z] += 1 - static_cast<intType>( mask( G( i, j, k - 1 ) ) );

                    // Check for cavity
                    intType numMaskedNeighbourCells = 0;
                    bool hasDirectionWithBothSidesMasked = false;
                    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                        numMaskedNeighbourCells += maskedCount[axis];
                        if ( maskedCount[axis] > 1 ) {
                            hasDirectionWithBothSidesMasked = true;
                        }
                    } );

                    bool hasCavity = hasDirectionWithBothSidesMasked && numMaskedNeighbourCells > 2;
                    if ( hasCavity ) {
                        mask( G(i, j, k) ) = CellType::Solid;
                        numCavities++;
                    }

                }
            }
        }

        if ( numCavities == 0 ) 
            break;

    } 
}


// Sets data for a particular source term in a particular direciton for a particular cell
void AddIBDataForDirection( IBCell &ibCell, 
                            const Axis::ENUMDATA axis,
                            const intType directionIndex,
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

    switch ( inputData.geoemtryBoundaryTreatment ) {
        case GeometryBoundaryTreatment::DirectionalImmersedBoundary:
        {
            // Distance from cell center to immersed boundary along this coordinate direction
            fVector3 queryPointCoords( mesh.cellCenters[X](cellIndex[X]),
                                       mesh.cellCenters[Y](cellIndex[Y]),
                                       mesh.cellCenters[Z](cellIndex[Z]) );
            fVector3 rayDirection( 0, 0, 0 );
            rayDirection[ axis ] = static_cast<floatType>( directionIndex );

            floatType ibDistance = NearestRayIntersection(tree, queryPointCoords, rayDirection);

            // Deal with some possible degenerate cases in the geometry
            const bool noIntersectionFound = ibDistance < 0.0f;

            TensorIndex3D cellIndex_gg = cellIndex_g;
            cellIndex_gg[axis] += directionIndex;
            const floatType maxAllowableIBDistance = mesh.cellLengths[axis](cellIndex[axis]) / 2.0f 
                                                   + mesh.cellLengths[axis](cellIndex_g[axis])
                                                   + mesh.cellLengths[axis](cellIndex_gg[axis]);
            const bool ibDistanceTooLarge  = ibDistance > maxAllowableIBDistance;

            if ( noIntersectionFound || ibDistanceTooLarge ) {
                sourceTermData.ibDistance = mesh.cellLengths[axis](cellIndex[axis]) / 2.0f;     // Take the nearest cell face
            } else {
                sourceTermData.ibDistance = ibDistance;
            }
            
            break;
        }
        
        case CAMIRA::GeometryBoundaryTreatment::Staircase:
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
                                     * mesh.cellFaceAreas[axis]( cellIndex[ LUT::LoOrthogonalAxis[axis] ], cellIndex[ LUT::HiOrthogonalAxis[axis] ] );    
    


    // Reconstruction coefficients from fluid onto face. May use further points due to stability condition
    bool meetsGhostStabilityCondition =  ibDistance >= ( mesh.cellLengths[axis](cellIndex[axis]) / 2.0f );
    if ( meetsGhostStabilityCondition ) {

        // Linear
        floatType dxpOn2 = mesh.cellLengths[axis](cellIndex[axis]) / 2.0f ;
        sourceTermData.faceReconstructionCoeff_p  = ( ibDistance - dxpOn2 ) / ( ibDistance );
        sourceTermData.faceReconstructionCoeff_a  = 0.0f;
        sourceTermData.faceReconstructionCoeff_ib = dxpOn2 / ibDistance;

        // // Quadratic
        // floatType dxp = mesh.cellLengths[axis](cellIndex[axis]);
        // floatType la  = abs( mesh.cellCenters[axis](cellIndex[axis]) - mesh.cellCenters[axis](cellIndex_a[axis]) );
        // sourceTermData.faceReconstructionCoeff_p  = - dxp * dxp / ( 4 * ibDistance * la )
        //                                           +   dxp * ( ibDistance - la ) / ( 2 * ibDistance * la );
        // sourceTermData.faceReconstructionCoeff_a  =   dxp * dxp /  ( 4 * la * ( ibDistance * la ) )
        //                                           -   dxp * ibDistance / ( 2 * la * ( ibDistance + la ) );
        // sourceTermData.faceReconstructionCoeff_ib =   dxp * dxp / ( 4 * ibDistance * ( ibDistance + la ) )
        //                                           +   dxp * la / ( 2 * ibDistance * ( ibDistance + la ) );


    } else {

        // Stencil skip
        floatType dxp = mesh.cellLengths[axis](cellIndex[axis]);
        floatType dxa = mesh.cellLengths[axis](cellIndex_a[axis]);
        floatType denominator =  dxp  +  dxa  +  2.0f * ibDistance;
        sourceTermData.faceReconstructionCoeff_p  = 0.0f;
        sourceTermData.faceReconstructionCoeff_a  = ( dxp - 2.0f * ibDistance ) / denominator;
        sourceTermData.faceReconstructionCoeff_ib = ( dxa + 2.0f * dxp ) / denominator;

        // // Nearest face
        // floatType dxp    = mesh.cellLengths[axis](cellIndex[axis]),
        //           lambda = mesh.interpFactors[axis]( fidx - directionIndex );
        // floatType denominator =  2.0f * ibDistance  + dxp;
        // if ( directionIndex == +1 ) {
        //     sourceTermData.faceReconstructionCoeff_p  = lambda          * ( 2.0f * ibDistance - dxp ) / denominator;
        //     sourceTermData.faceReconstructionCoeff_a  = (1.0f - lambda) * ( 2.0f * ibDistance - dxp ) / denominator;
        // } else {
        //     sourceTermData.faceReconstructionCoeff_p  = (1.0f - lambda) * ( 2.0f * ibDistance - dxp ) / denominator;
        //     sourceTermData.faceReconstructionCoeff_a  = lambda          * ( 2.0f * ibDistance - dxp ) / denominator;
        // }
        // sourceTermData.faceReconstructionCoeff_ib = 2.0f * dxp / denominator;

        // // Just set to zero
        // sourceTermData.faceReconstructionCoeff_p  = 0.0f;
        // sourceTermData.faceReconstructionCoeff_a  = 0.0f; 
        // sourceTermData.faceReconstructionCoeff_ib = 1.0f;
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
    sourceTermData.ibExtrapCoeff_p =   ibDistance / cellInteriorDistance  +  1.0f;
    sourceTermData.ibExtrapCoeff_a = - ibDistance / cellInteriorDistance;


    // Extrapolation coefficients from fluid to face (Does not account for presence of IB)
    sourceTermData.faceExtrapCoeff_p =   0.05f * mesh.cellLengths[axis](cellIndex[axis]) / cellInteriorDistance  +  1.0f;
    sourceTermData.faceExtrapCoeff_a = - 0.05f * mesh.cellLengths[axis](cellIndex[axis]) / cellInteriorDistance;


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



void SetVelocityFluxCorrectionCoefficient( std::vector<IBCell> &ibCellsComponent,
                                           const Mesh &mesh )
{

    // The denominator must be calculated seperately first
    floatType denominator = 0.0f;
    for ( auto &ibCell : ibCellsComponent ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            floatType ibCellFaceDistance = abs( sourceTermData.ibDistance - mesh.cellLengths[sourceTermData.direction](ibCell.cellIndex[sourceTermData.direction]) );

            denominator += std::pow( abs( sourceTermData.faceAreaComponent * ibCellFaceDistance) , 2.0f );

        }
    }

    // Now the coefficient for each cell
    for ( auto &ibCell : ibCellsComponent ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            floatType ibCellFaceDistance = abs( sourceTermData.ibDistance - mesh.cellLengths[sourceTermData.direction](ibCell.cellIndex[sourceTermData.direction]) );

            sourceTermData.velocityFluxCorrectionCoeff = ( std::pow( ibCellFaceDistance, 2.0f ) 
                                                        * sourceTermData.faceAreaComponent 
                                                        ) / denominator;

        }
    }

}



std::vector<IBCell> CreateIBCellDataForComponent( const Tensor3D &mask,
                                                  const Tree &tree,
                                                  const Mesh &mesh,
                                                  const InputData &inputData )
{

    using enum Axis::ENUMDATA;
    using FVT::G;

    std::vector<IBCell> ibCellsComponent;

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
                            ibCellsComponent.emplace_back();
                            ibCellPtr = &ibCellsComponent.back();
                            ibCellPtr->cellIndex = cellIndex;
                        } 
                    };

                    // Solid on hi side
                    TensorIndex3D hiSideCellIndex = cellIndex;
                    hiSideCellIndex[axis] += 1;
                    bool atHiBoundary = ( cellIndex[axis] == mesh.nCells[axis]-1  );
                    if ( !atHiBoundary && static_cast<intType>( mask(G(hiSideCellIndex)) ) == CellType::Solid ) {
                        CheckIBCellPtr();
                        AddIBDataForDirection( *ibCellPtr, axis, +1, mesh, tree, inputData );
                    }

                    // Solid on lo side
                    TensorIndex3D loSideCellIndex = cellIndex;
                    loSideCellIndex[axis] -= 1;
                    bool atLoBoundary = ( cellIndex[axis] == 0  );
                    if ( !atLoBoundary && static_cast<intType>( mask(G(loSideCellIndex)) ) == CellType::Solid ) {
                        CheckIBCellPtr();
                        AddIBDataForDirection( *ibCellPtr, axis, -1, mesh, tree, inputData );
                    }

                } );

            }
        }
    }

    // Calculate the coefficients for the velocity flux correction
    SetVelocityFluxCorrectionCoefficient( ibCellsComponent, mesh );

    return ibCellsComponent;
}


}   // end anonymous namespace




void SetImmersedBoundaryData( IBData &ibData,
                              const InputData &inputData,
                              const AxisTransformationMap &axisTransformation,
                              const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    // Create the initial mask. Ghost cells at domain boundaries should be masked out.
    const TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                        extents = {mesh.nCells[0], mesh.nCells[1], mesh.nCells[2]};
    ibData.mask = Tensor3D( mesh.nCells[X] + 2*CAMIRA::nGhost, mesh.nCells[Y] + 2*CAMIRA::nGhost, mesh.nCells[Z] + 2*CAMIRA::nGhost );
    SetTensorConstantParallel( ibData.mask, CellType::Solid );
    ibData.mask.slice(offsets, extents) = ibData.mask.slice(offsets, extents).constant( CellType::Fluid );

    // Leave ibData empty if there is no geometry
    if ( !inputData.hasIBGeometry )
        return;

    Polyhedron geometry = MakeGeometry( inputData, axisTransformation );

    // Separate the geometry into connected components
    std::vector<Polyhedron> polyVector = SeparatePolyhedron( geometry );

    // Set the IBcells for each one
    // ibData.ibCells.clear();
    for ( const Polyhedron &poly : polyVector ) {

        Tree tree = MakeAABBTree( poly );

        // Local mask for just this component
        Tensor3D localMask = CreateCellMask( tree, mesh );

        // Add the contribution to the global mask
        ibData.mask *= localMask;

        // Set ibCells for this component
        ibData.ibCells.emplace_back( CreateIBCellDataForComponent( localMask, tree, mesh, inputData ) );

    }

    // Remove single cell cavities in the mask, which can cause instabilities
    RemoveMaskSingleCellCavities( ibData.mask, mesh );
}



void WriteGeometryToFile( const InputData &inputData, 
                          const AxisTransformationMap &axisTransformation )
{
    InputData inputDataUserCoordinates( inputData );
    TransformUserInputData( inputDataUserCoordinates, axisTransformation.Inverse() );
    Polyhedron PUserCoordinates = MakeGeometry( inputDataUserCoordinates, axisTransformation.Identity() );

    std::ofstream out(  IOTOOLS::RemoveFileExtension( inputData.geometryOutputFilename, ".stl" ) + ".stl" );
    CGAL::IO::write_STL( out, PUserCoordinates );
}



}   // end namespace CAMIRA