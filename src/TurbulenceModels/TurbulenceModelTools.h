#ifndef CAMIRA_TURBULENCE_MODEL_TOOLS
#define CAMIRA_TURBULENCE_MODEL_TOOLS

#include "../Core/Types.h"
#include "../Core/FVLookups.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Geometry/Geometry.h"

#include <algorithm>

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

namespace detail
{

inline bool IsZeroArray( const Tensor2D &array )
{
    floatType testValue = 0.0f;
    const Eigen::Tensor<bool, 0> isConstant = ( array == testValue ).all();
    return isConstant(0);
}



inline bool IsWallBC( const BoundaryConditionData &bcData,
               const BoundaryPatches::ENUMDATA &bp)
{
    // Is wall BC if:
    // Pressure boundary condition is either extrapolated or zeroGradient
    // No velocity component in wall normal direction (can be a moving wall)

    using enum Axis::ENUMDATA;
    Axis::ENUMDATA wallNormal = LUT::BoundaryPatchAxis[ bp ];

    bool isValidWallPressureBC =  bcData.fields.P[bp].type == BoundaryConditions::zeroGradient
                               || bcData.fields.P[bp].type == BoundaryConditions::extrapolated;
    if ( !isValidWallPressureBC )
        return false; 

    bool allVelocityBCsFixed =  bcData.fields.U[X][bp].type == BoundaryConditions::fixed
                             && bcData.fields.U[Y][bp].type == BoundaryConditions::fixed
                             && bcData.fields.U[Z][bp].type == BoundaryConditions::fixed;

    bool wallNormalVelocityIsZero = detail::IsZeroArray( bcData.fields.U[wallNormal][bp].value );

    if ( allVelocityBCsFixed && wallNormalVelocityIsZero )
        return true;

    return false;
}

}   // end namespace detail



// Rate of strain tensor squared 
inline void CalculateVelocityDeformationRate( Tensor3D &S,
                                              const FieldData<Tensor3D> &fields, 
                                              const IBData &ibData,
                                              const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                TensorIndex3D cellP = G(i  , j  , k  ),
                              cellN = G(j  , j+1, k  ),
                              cellS = G(i  , j-1, k  ),
                              cellE = G(i+1, j  , k  ),
                              cellW = G(i-1, j  , k  ),
                              cellT = G(i  , j  , k+1),
                              cellB = G(i  , j  , k-1);

                EnumVector<Axis, floatType>  ddxU, ddyU, ddzU;

                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                    ddxU[axis] =  0.5 * (   mesh.cellCenterDiffInv[X](i+1)                                   *  fields.U[axis](cellE)
                                        + ( mesh.cellCenterDiffInv[X](i) - mesh.cellCenterDiffInv[X](i+1) )  *  fields.U[axis](cellP) 
                                        -   mesh.cellCenterDiffInv[X](i)                                     *  fields.U[axis](cellW) );


                    ddyU[axis] =  0.5 * (   mesh.cellCenterDiffInv[Y](j+1)                                   *  fields.U[axis](cellN)
                                        + ( mesh.cellCenterDiffInv[Y](j) - mesh.cellCenterDiffInv[Y](j+1) )  *  fields.U[axis](cellP) 
                                        -   mesh.cellCenterDiffInv[Y](j)                                     *  fields.U[axis](cellS) );

                    ddzU[axis] =  0.5 * (   mesh.cellCenterDiffInv[Z](k+1)                                   *  fields.U[axis](cellT)
                                        + ( mesh.cellCenterDiffInv[Z](k) - mesh.cellCenterDiffInv[Z](k+1) )  *  fields.U[axis](cellP) 
                                        -   mesh.cellCenterDiffInv[Z](k)                                     *  fields.U[axis](cellB) );

                } );

                S(cellP)    = sqrt(
                               ddxU[X] * ddxU[X]
                             + ddyU[Y] * ddyU[Y]
                             + ddzU[Z] * ddzU[Z]
                             + 0.5 * ( ddyU[X] + ddxU[Y] ) * ( ddyU[X] + ddxU[Y] )
                             + 0.5 * ( ddzU[X] + ddxU[Z] ) * ( ddzU[X] + ddxU[Z] )
                             + 0.5 * ( ddzU[Y] + ddyU[Z] ) * ( ddzU[Y] + ddyU[Z] )
                            );

            }
        }
    }    


    // Go back through and correct for the immersed boundary faces
    for ( const auto &ibCellComponent : ibData.ibCells ) {
        for ( const auto &ibCell : ibCellComponent ) {

            const intType i = ibCell.cellIndex[X];
            const intType j = ibCell.cellIndex[Y];
            const intType k = ibCell.cellIndex[Z];

            TensorIndex3D cellP = G(i  , j  , k  ),
                          cellN = G(j  , j+1, k  ),
                          cellS = G(i  , j-1, k  ),
                          cellE = G(i+1, j  , k  ),
                          cellW = G(i-1, j  , k  ),
                          cellT = G(i  , j  , k+1),
                          cellB = G(i  , j  , k-1);


            // Set the ghost cell values if they need to be
            EnumVector<Axis, floatType> velocityE, velocityW, 
                                        velocityN, velocityS, 
                                        velocityT, velocityB;

            // Set them to the field value initially
            EnumFor<Axis>( [&] (Axis::ENUMDATA velocityComponent) {
                velocityE[velocityComponent] = fields.U[velocityComponent](cellE);
                velocityW[velocityComponent] = fields.U[velocityComponent](cellW);
                velocityN[velocityComponent] = fields.U[velocityComponent](cellN);
                velocityS[velocityComponent] = fields.U[velocityComponent](cellS);
                velocityT[velocityComponent] = fields.U[velocityComponent](cellT);
                velocityB[velocityComponent] = fields.U[velocityComponent](cellB);
            } );

            // Modify any that are an IB ghost cell
            for ( const auto &sourceTermData : ibCell.sourceTermsData ) {
                const Axis::ENUMDATA faceNormal = sourceTermData.direction;

                if ( sourceTermData.directionIndex > 0 ) {  // Ghost cell on hi side
                    switch ( faceNormal ) {
                        case (Axis::X):
                            EnumFor<Axis>( [&] (Axis::ENUMDATA c) { velocityE[c] = sourceTermData.ghostCellValues.U[c]; } );
                            break;
                        case (Axis::Y):
                            EnumFor<Axis>( [&] (Axis::ENUMDATA c) { velocityN[c] = sourceTermData.ghostCellValues.U[c]; } );
                            break;
                        case (Axis::Z):
                            EnumFor<Axis>( [&] (Axis::ENUMDATA c) { velocityT[c] = sourceTermData.ghostCellValues.U[c]; } );
                            break;
                    }
                } else {                                    // Ghost cell on lo side
                    switch ( faceNormal ) {
                        case (Axis::X):
                            EnumFor<Axis>( [&] (Axis::ENUMDATA c) { velocityW[c] = sourceTermData.ghostCellValues.U[c]; } );
                            break;
                        case (Axis::Y):
                            EnumFor<Axis>( [&] (Axis::ENUMDATA c) { velocityS[c] = sourceTermData.ghostCellValues.U[c]; } );
                            break;
                        case (Axis::Z):
                            EnumFor<Axis>( [&] (Axis::ENUMDATA c) { velocityB[c] = sourceTermData.ghostCellValues.U[c]; } );
                            break;
                    }
                }

            }

            EnumVector<Axis, floatType>  ddxU, ddyU, ddzU;
            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                ddxU[axis] =  0.5 * (   mesh.cellCenterDiffInv[X](i+1)                                   *  velocityE[axis]
                                    + ( mesh.cellCenterDiffInv[X](i) - mesh.cellCenterDiffInv[X](i+1) )  *  fields.U[axis](cellP) 
                                    -   mesh.cellCenterDiffInv[X](i)                                     *  velocityW[axis]       );


                ddyU[axis] =  0.5 * (   mesh.cellCenterDiffInv[Y](j+1)                                   *  velocityN[axis]
                                    + ( mesh.cellCenterDiffInv[Y](j) - mesh.cellCenterDiffInv[Y](j+1) )  *  fields.U[axis](cellP) 
                                    -   mesh.cellCenterDiffInv[Y](j)                                     *  velocityS[axis]       );

                ddzU[axis] =  0.5 * (   mesh.cellCenterDiffInv[Z](k+1)                                   *  velocityT[axis]
                                    + ( mesh.cellCenterDiffInv[Z](k) - mesh.cellCenterDiffInv[Z](k+1) )  *  fields.U[axis](cellP) 
                                    -   mesh.cellCenterDiffInv[Z](k)                                     *  velocityB[axis]       );

            } );

            S(cellP)    = sqrt(
                            ddxU[X] * ddxU[X]
                            + ddyU[Y] * ddyU[Y]
                            + ddzU[Z] * ddzU[Z]
                            + 0.5 * ( ddyU[X] + ddxU[Y] ) * ( ddyU[X] + ddxU[Y] )
                            + 0.5 * ( ddzU[X] + ddxU[Z] ) * ( ddzU[X] + ddxU[Z] )
                            + 0.5 * ( ddzU[Y] + ddyU[Z] ) * ( ddzU[Y] + ddyU[Z] )
                        );
        }
    }


    // Set domain boundary ghost cells - set so cell face equals cell centre value
    for ( intType i = 0; i != BoundaryPatches::count; i++ ) {
        auto boundaryPatch = static_cast<BoundaryPatches::ENUMDATA>(i);
        
        Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[ boundaryPatch ];

        intType iCell_p = ( boundaryPatch == LUT::NegativePatch[axis] ) ? 0  : mesh.nCells[axis] - 1,
                iCell_g = ( boundaryPatch == LUT::NegativePatch[axis] ) ? -1 : mesh.nCells[axis];


        S.chip( G(iCell_g), axis ) = S.chip( G(iCell_p) , axis );
    }

}



inline EnumVector<Axis, Tensor3D> NearestWallDistance( const Mesh &mesh,
                                                       const Tree &tree,
                                                       const BoundaryConditionData &bcData )
{
    using enum Axis::ENUMDATA;

    EnumVector<Axis, Tensor3D> wallDistance( { Tensor3D( mesh.nFacesNormal[X][X], mesh.nFacesNormal[X][Y], mesh.nFacesNormal[X][Z] ).setZero(),
                                               Tensor3D( mesh.nFacesNormal[Y][X], mesh.nFacesNormal[Y][Y], mesh.nFacesNormal[Y][Z] ).setZero(),
                                               Tensor3D( mesh.nFacesNormal[Z][X], mesh.nFacesNormal[Z][Y], mesh.nFacesNormal[Z][Z] ).setZero() } );

    // Work out which domain boundaries are walls and store the result
    EnumVector<BoundaryPatches, bool> domainBoundaryWalls;
    
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        if ( detail::IsWallBC(bcData, bp)  ) { 
            domainBoundaryWalls[bp] = true;
        } else {
            domainBoundaryWalls[bp] = false;
        } 
    } );


    floatType inf = std::numeric_limits<floatType>::infinity();

    EnumFor<Axis>( [&] (Axis::ENUMDATA faceNormal) {

        for ( intType k = 0; k != mesh.nFacesNormal[faceNormal][Z]; k++ ) {
            for ( intType j = 0; j != mesh.nFacesNormal[faceNormal][Y]; j++ ) {
                for ( intType i = 0; i != mesh.nFacesNormal[faceNormal][X]; i++ ) {

                    TensorIndex3D faceIndex{i, j, k};
                    fVector3 faceCoords{ (faceNormal == X) ? mesh.cellFaces[X](i) : mesh.cellCenters[X](i),
                                         (faceNormal == Y) ? mesh.cellFaces[Y](j) : mesh.cellCenters[Y](j),
                                         (faceNormal == Z) ? mesh.cellFaces[Z](k) : mesh.cellCenters[Z](k), };

                    // Ignore if this point is inside the geometry
                    // Commented out since immersed boundary faces can be inside the solid, but may still need a distance to compute the eddy viscosity
                    // if ( PointInside( tree, faceCoords(X), faceCoords(Y), faceCoords(Z) ) ) 
                    //     continue; 

                    // Distance to the nearest geometry object
                    wallDistance[faceNormal](faceIndex) = inf;
                    if ( !tree.empty() ){
                        wallDistance[faceNormal](faceIndex) = NearestDistance( tree,
                                                                               faceCoords[X],
                                                                               faceCoords[Y],
                                                                               faceCoords[Z] );
                    }

                    // Check if there is a domain boundary that is closer
                    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                        
                        if ( domainBoundaryWalls[ LUT::NegativePatch[axis] ] ) {
                            wallDistance[faceNormal](faceIndex) = std::min( wallDistance[faceNormal](faceIndex),
                                                                            abs( faceCoords[axis] - mesh.cellFaces[axis]( 0 ) ) );
                        }

                        if ( domainBoundaryWalls[ LUT::PositivePatch[axis] ] ) {
                            wallDistance[faceNormal](faceIndex) = std::min( wallDistance[faceNormal](faceIndex),
                                                                            abs( mesh.cellFaces[axis]( mesh.nFacesNormal[axis](axis) - 1 ) - faceCoords[axis] ) );
                        }

                    } );

                }
            }
        }

    } );

    
    return wallDistance;
} 


}   // end namespace CAMIRA


#endif // CAMIRA_TURBULENCE_MODEL_TOOLS