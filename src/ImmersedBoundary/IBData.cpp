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



// Extend the cell mask so that face centers on the AB are always outside the solid boundary
void ExtendCellMask( Tensor3D &mask,
                     const Tree &tree, 
                     const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    // Identify cells with cell centeres inside the solid
    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                // Only changing fluid points
                if ( static_cast<intType>( mask(G(i, j, k)) ) == CellType::Solid )
                    continue;

                TensorIndex3D cellIndex = {i, j, k};

                // Mask this cell if any of its face centers are inside the solid
                bool hasFaceInsideSolid = false;

                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                    fVector3 faceCoordinates;

                    // Positive side face
                    EnumFor<Axis>( [&] (Axis::ENUMDATA a) {
                        faceCoordinates[a] = ( a == axis ) ? mesh.cellFaces[a]( cellIndex[a] + 1 ) : mesh.cellCenters[a]( cellIndex[a] );
                    } );
                    
                    if ( PointInside( tree, faceCoordinates[X], faceCoordinates[Y], faceCoordinates[Z] )  )
                        hasFaceInsideSolid = true;

                    
                    // Negative side face
                    faceCoordinates[axis] = mesh.cellFaces[axis]( cellIndex[axis] );
                    
                    if ( PointInside( tree, faceCoordinates[X], faceCoordinates[Y], faceCoordinates[Z] )  )
                        hasFaceInsideSolid = true;

                } );

                if ( hasFaceInsideSolid )   
                    mask(G(i, j, k)) = CellType::Solid;

            }
        }
    }

}



// Masks cells that form a single cell "cavity". These can cause instability and are a highly degenerate case
void RemoveMaskSingleCellCavities( Tensor3D &mask,
                                   const Mesh &mesh )
{
    using FVT::G;
    using enum Axis::ENUMDATA;

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

                if ( hasDirectionWithBothSidesMasked && numMaskedNeighbourCells > 2 )
                    mask( G(i, j, k) ) = CellType::Solid;

            }
        }
    }
}



// Check if a given index is within the domain bounds
bool OutOfBounds( const TensorIndex3D &index,
                  const Mesh &mesh )
{
    // Can't use EnumFor since return statements inside loop
    for ( int a = 0; a != Axis::count; a++ ) {  
        Axis::ENUMDATA axis = static_cast<Axis::ENUMDATA>(a);

        if ( index[axis] < 0 )
            return true;

        if ( index[axis] > mesh.nCells[axis]-1 )
            return true;
    }

    return false;
}



floatType DiagonalCellDistance( const TensorIndex3D &cellIndex,
                                const Mesh &mesh )
{
    floatType distanceSq = 0.0;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        distanceSq += pow( mesh.cellLengths[axis]( cellIndex[axis] ), 2 );
    } );
    return sqrt( distanceSq );
}



// Determine the length of the image point from the immersed boundary. We take this as the maximum 
// diagonal length of the cells surrounding the given cell
floatType GetImagePointDistance( const TensorIndex3D &cellIndex, 
                                 const Mesh &mesh )
{
    floatType maxValue = 0.0f;
    std::array<intType, 3> deltaIndex = {-1, 0, 1};

    for ( intType dk : deltaIndex ) {
        for ( intType dj : deltaIndex ) {
            for ( intType di : deltaIndex ) {

                bool isNeighbour = ( dk != 0 ) || ( dj != 0 ) || ( di != 0 );
                if ( !isNeighbour ) {
                    continue;
                }

                TensorIndex3D neighbourIndex = cellIndex;
                neighbourIndex[0] += di;
                neighbourIndex[1] += dj;
                neighbourIndex[2] += dk;

                if ( OutOfBounds( neighbourIndex, mesh ) ) {
                    continue;
                }

                floatType diagonalCellDistance = DiagonalCellDistance( neighbourIndex, mesh ); 

                if ( diagonalCellDistance >= maxValue ) {
                    maxValue = diagonalCellDistance;
                }

            }
        }
    }

    return maxValue;
}



// Sets data for Immersed Boundary Method with wall functions 
void AddWallFunctionIBDataForDirection( IBCell::SourceTermData &sourceTermData, 
                                        const TensorIndex3D &cellIndex,
                                        const Mesh &mesh,
                                        const Tree &tree )
{
    using enum Axis::ENUMDATA;

    sourceTermData.wallFunctionDataPtr = std::make_unique<IBCell::WallFunctionData>();
    auto &wallFunctionData = (*sourceTermData.wallFunctionDataPtr);

    // Face index
    TensorIndex3D faceIndex = cellIndex;
    faceIndex[sourceTermData.direction] += sourceTermData.faceDirectionIndex;

    // Face coordiantes
    fVector3 faceCoordinates;
    EnumFor<Axis>( [&] (Axis::ENUMDATA a) {
        faceCoordinates[a] = ( sourceTermData.direction == a ) ? mesh.cellFaces[a]( faceIndex[a] ) : mesh.cellCenters[a]( faceIndex[a] );
    } );


    // Nearest to the face on the immersed boundary
    fVector3 nearestIBPoint = NearestPoint( tree, faceCoordinates[X], faceCoordinates[Y], faceCoordinates[Z] );

    // Normal vector
    wallFunctionData.normalVector = faceCoordinates - nearestIBPoint;

    // Nearest point distnace
    wallFunctionData.faceCenterDistance = wallFunctionData.normalVector.norm();

    // Normalise the normal vector
    wallFunctionData.normalVector.normalize();

    // Distance to image point
    wallFunctionData.imagePointDistance = GetImagePointDistance( cellIndex, mesh );

    // Image point coordinates
    fVector3 imagePoint = faceCoordinates + wallFunctionData.imagePointDistance * wallFunctionData.normalVector; 

    // Construct field probe at image point
    wallFunctionData.fieldProbePtr = std::make_unique<FieldProbe>( mesh, imagePoint );

    // Initialise y+
    wallFunctionData.yPlusImagePoint = 0.0f;

}



// Sets data for Immersed Boundary Method without wall functions
void AddNoWallFunctionIBDataForDirection( IBCell::SourceTermData &sourceTermData, 
                                          const TensorIndex3D &cellIndex,
                                          const Axis::ENUMDATA axis,
                                          const intType directionIndex,
                                          const Mesh &mesh,
                                          const Tree &tree, 
                                          const InputData &inputData )
{
    using enum Axis::ENUMDATA;

    sourceTermData.directionalIBDataPtr = std::make_unique<IBCell::DirectionalIBData>();

    auto &directionalIBData = *sourceTermData.directionalIBDataPtr;

    TensorIndex3D &cellIndex_a = sourceTermData.cellIndex_a;
    TensorIndex3D &cellIndex_g = sourceTermData.cellIndex_g;

    switch ( inputData.geoemtryBoundaryTreatement ) {
        case GeometryBoundaryTreatement::ImmersedBoundary:
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
                directionalIBData.ibDistance = mesh.cellLengths[axis](cellIndex[axis]) / 2.0f;     // Take the nearest cell face
            } else {
                directionalIBData.ibDistance = ibDistance;
            }
            
            break;
        }
        
        case CAMIRA::GeometryBoundaryTreatement::Staircase:
        {   
            // Distance to the nearest cell face, approximates the immersed boundary to be on the cell face
            // i.e. staircase approximation
            directionalIBData.ibDistance = mesh.cellLengths[axis](cellIndex[axis]) / 2.0f;
            break;
        }

        default:
            break;
    }
    const floatType &ibDistance = directionalIBData.ibDistance;


    // Reconstruction coefficients from fluid onto face. May use further points due to stability condition
    bool meetsGhostStabilityCondition =  ibDistance >= ( mesh.cellLengths[axis](cellIndex[axis]) / 2.0f );
    if ( meetsGhostStabilityCondition ) {

        // Linear
        floatType dxpOn2 = mesh.cellLengths[axis](cellIndex[axis]) / 2.0f ;
        directionalIBData.faceReconstructionCoeff_p  = ( ibDistance - dxpOn2 ) / ( ibDistance );
        directionalIBData.faceReconstructionCoeff_a  = 0.0f;
        directionalIBData.faceReconstructionCoeff_ib = dxpOn2 / ibDistance;

        // // Quadratic
        // floatType dxp = mesh.cellLengths[axis](cellIndex[axis]);
        // floatType la  = abs( mesh.cellCenters[axis](cellIndex[axis]) - mesh.cellCenters[axis](cellIndex_a[axis]) );
        // directionalIBData.faceReconstructionCoeff_p  = - dxp * dxp / ( 4 * ibDistance * la )
        //                                           +   dxp * ( ibDistance - la ) / ( 2 * ibDistance * la );
        // directionalIBData.faceReconstructionCoeff_a  =   dxp * dxp /  ( 4 * la * ( ibDistance * la ) )
        //                                           -   dxp * ibDistance / ( 2 * la * ( ibDistance + la ) );
        // directionalIBData.faceReconstructionCoeff_ib =   dxp * dxp / ( 4 * ibDistance * ( ibDistance + la ) )
        //                                           +   dxp * la / ( 2 * ibDistance * ( ibDistance + la ) );


    } else {

        // Stencil skip
        floatType dxp = mesh.cellLengths[axis](cellIndex[axis]);
        floatType dxa = mesh.cellLengths[axis](cellIndex_a[axis]);
        floatType denominator =  dxp  +  dxa  +  2.0f * ibDistance;
        directionalIBData.faceReconstructionCoeff_p  = 0.0f;
        directionalIBData.faceReconstructionCoeff_a  = ( dxp - 2.0f * ibDistance ) / denominator;
        directionalIBData.faceReconstructionCoeff_ib = ( dxa + 2.0f * dxp ) / denominator;

        // // Nearest face
        // floatType dxp    = mesh.cellLengths[axis](cellIndex[axis]),
        //           lambda = mesh.interpFactors[axis]( fidx - directionIndex );
        // floatType denominator =  2.0f * ibDistance  + dxp;
        // if ( directionIndex == +1 ) {
        //     directionalIBData.faceReconstructionCoeff_p  = lambda          * ( 2.0f * ibDistance - dxp ) / denominator;
        //     directionalIBData.faceReconstructionCoeff_a  = (1.0f - lambda) * ( 2.0f * ibDistance - dxp ) / denominator;
        // } else {
        //     directionalIBData.faceReconstructionCoeff_p  = (1.0f - lambda) * ( 2.0f * ibDistance - dxp ) / denominator;
        //     directionalIBData.faceReconstructionCoeff_a  = lambda          * ( 2.0f * ibDistance - dxp ) / denominator;
        // }
        // directionalIBData.faceReconstructionCoeff_ib = 2.0f * dxp / denominator;

        // // Just set to zero
        // directionalIBData.faceReconstructionCoeff_p  = 0.0f;
        // directionalIBData.faceReconstructionCoeff_a  = 0.0f; 
        // directionalIBData.faceReconstructionCoeff_ib = 1.0f;
    }


    // Extrapolation coefficients from fluid to immersed boundary surface
    floatType cellInteriorDistance = abs( mesh.cellCenters[axis](cellIndex_a[axis]) - mesh.cellCenters[axis](cellIndex[axis]) ); 
    directionalIBData.ibExtrapCoeff_p =   ibDistance / cellInteriorDistance  +  1.0f;
    directionalIBData.ibExtrapCoeff_a = - ibDistance / cellInteriorDistance;


    // Extrapolation coefficients from fluid to face (Does not account for presence of IB)
    directionalIBData.faceExtrapCoeff_p =   0.05f * mesh.cellLengths[axis](cellIndex[axis]) / cellInteriorDistance  +  1.0f;
    directionalIBData.faceExtrapCoeff_a = - 0.05f * mesh.cellLengths[axis](cellIndex[axis]) / cellInteriorDistance;
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

    // Set data that is specific to whether or not wall functions are used
    if ( inputData.useWallFunctions ) {

        AddWallFunctionIBDataForDirection( sourceTermData, cellIndex, mesh, tree );

    } else {

        AddNoWallFunctionIBDataForDirection( sourceTermData, cellIndex, axis, directionIndex, mesh, tree, inputData );
    }


    // Face area vector
    sourceTermData.faceAreaComponent = static_cast<floatType>( directionIndex )   // Gives the correct sign
                                     * mesh.cellFaceAreas[axis]( cellIndex[ LUT::LoOrthogonalAxis[axis] ], cellIndex[ LUT::HiOrthogonalAxis[axis] ] );    
    
    // Extrapolation coefficients from face to ghost cell
    if        ( directionIndex == +1 ) {    // Face on Hi side

        sourceTermData.ghostExtrapCoeff_p = - (1 - mesh.interpFactors[axis](fidx)) / mesh.interpFactors[axis](fidx);
        sourceTermData.ghostExtrapCoeff_f = 1.0f / mesh.interpFactors[axis](fidx);

    } else if ( directionIndex == -1 ) {   // Face on Lo side

        sourceTermData.ghostExtrapCoeff_p = - mesh.interpFactors[axis](fidx) / (1.0f - mesh.interpFactors[axis](fidx));
        sourceTermData.ghostExtrapCoeff_f = 1.0f / (1.0f - mesh.interpFactors[axis](fidx));

    }

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

            floatType ibCellFaceDistance = 0.0f;
            if ( sourceTermData.directionalIBDataPtr ) {
                ibCellFaceDistance = abs( sourceTermData.directionalIBDataPtr->ibDistance - mesh.cellLengths[sourceTermData.direction](ibCell.cellIndex[sourceTermData.direction]) );
            } else if ( sourceTermData.wallFunctionDataPtr ) {
                ibCellFaceDistance = sourceTermData.wallFunctionDataPtr->faceCenterDistance;
            }

            denominator += std::pow( abs( sourceTermData.faceAreaComponent * ibCellFaceDistance) , 2.0f );

        }
    }

    // Now the coefficient for each cell
    for ( auto &ibCell : ibCellsComponent ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            floatType ibCellFaceDistance = 0.0f;
            if ( sourceTermData.directionalIBDataPtr ) {
                ibCellFaceDistance = abs( sourceTermData.directionalIBDataPtr->ibDistance - mesh.cellLengths[sourceTermData.direction](ibCell.cellIndex[sourceTermData.direction]) );
            } else if ( sourceTermData.wallFunctionDataPtr ) {
                ibCellFaceDistance = sourceTermData.wallFunctionDataPtr->faceCenterDistance;
            }
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
    
    ibData.useWallFunctions = inputData.useWallFunctions;
    ibData.rho = inputData.rho;
    ibData.nu  = inputData.nu;

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

        if ( inputData.useWallFunctions )
            ExtendCellMask( localMask, tree, mesh );

        RemoveMaskSingleCellCavities( localMask, mesh );

        // Add the contribution to the global mask
        ibData.mask *= localMask;

        // Set ibCells for this component
        ibData.ibCells.emplace_back( CreateIBCellDataForComponent( localMask, tree, mesh, inputData ) );

    }
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