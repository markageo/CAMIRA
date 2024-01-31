#include "ImmersedBoundary.h"
#include "../Tools/FVLookups.h"

#include "../IO/ArrayIO.h"

#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/squared_distance_3.h> 

#include <cmath>

#include <iostream>

namespace CFD
{

namespace 
{

// Determines if query point is inside given polyhedron geometry
bool PointInside( const Polyhedron &polyhedron, 
                  const floatType xq, 
                  const floatType yq, 
                  const floatType zq ) 
{
    using Point        = Polyhedron::Point_3;
    using Primitive    = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
    using Traits       = CGAL::AABB_traits<CGAL_Kernel, Primitive>;
    using Tree         = CGAL::AABB_tree<Traits>;
    using Point_inside = CGAL::Side_of_triangle_mesh<Polyhedron, CGAL_Kernel>;

    // Construct AABB tree with a KdTree
    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.accelerate_distance_queries();

    // Initialize the point-in-polyhedron tester
    Point_inside inside_tester(tree);
    
    // Determine the side and return true if inside!
    return inside_tester( Point( xq, yq, zq ) ) == CGAL::ON_BOUNDED_SIDE;
}



// Returns the distance squared to the nearest intersection from the given point in the direction of the given ray.
// https://stackoverflow.com/questions/69953358/to-calculate-intersections-between-a-ray-and-a-mesh
floatType GetBoundaryDistance2( const Polyhedron &polyhedron,
                                const fVector3 &queryPointCoords,
                                const fVector3 &rayDirection ) 
{
    using Point            = Polyhedron::Point_3;
    using Vector           = CGAL_Kernel::Vector_3;
    using Primitive        = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
    using Traits           = CGAL::AABB_traits<CGAL_Kernel, Primitive>;
    using Tree             = CGAL::AABB_tree<Traits>;
    using Ray              = CGAL_Kernel::Ray_3;    
    using Ray_intersection = boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type>;

    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.accelerate_distance_queries();

    const Vector rayOrientation( rayDirection(0), rayDirection(1), rayDirection(2) );
    const Point  rayOrigin( queryPointCoords(0), queryPointCoords(1), queryPointCoords(2) );
    const Ray    ray( rayOrigin, rayOrientation );
    Ray_intersection intersection = tree.first_intersection( ray );

    Point* intersectionPoint = boost::get<Point>( &(intersection->first) ); 
    floatType distance2 = static_cast<floatType>( CGAL::squared_distance( *intersectionPoint, rayOrigin ) );

    return distance2;
}



// Create a mask array for cell centers
Tensor3D CreateCellMask( const Polyhedron &geometry, 
                         const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    Tensor3D mask = Tensor3D( mesh.nCells[X], mesh.nCells[Y], mesh.nCells[Z] ).setConstant( CellType::Fluid );

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                floatType xq = mesh.cellCenters[X](i),
                          yq = mesh.cellCenters[Y](j),
                          zq = mesh.cellCenters[Z](k);

                if ( PointInside( geometry, xq, yq, zq ) ) {
                    mask(i, j, k) = CellType::Solid;
                } 

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

    if ( static_cast<intType>( mask( cellIndex_a ) ) == CellType::Solid ) 
        return false;

    for ( intType axis = 0; axis != Axis::count; axis++ ) {

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
                            const Polyhedron &geometry)
{
    using enum Axis::ENUMDATA;

    ibCell.sourceTermsData.emplace_back();
    IBCell::SourceTermData &sourceTermData = ibCell.sourceTermsData.back();

    TensorIndex3D &cellIndex = ibCell.cellIndex;

    TensorIndex3D ghostCellIndex    = cellIndex;
    ghostCellIndex[axis] += directionIndex;

    TensorIndex3D interiorCellIndex = cellIndex;
    interiorCellIndex[axis] -= directionIndex;

    sourceTermData.direction          = axis;
    sourceTermData.directionIndex     = directionIndex;
    sourceTermData.faceDirectionIndex = ( directionIndex == 1 ) ? 1 : 0 ;
    sourceTermData.cellIndex_a        = interiorCellIndex;

    intType fidx = cellIndex[axis] + sourceTermData.faceDirectionIndex;

    // Ensure that there is enough space between the IB and other IBs and the domain boundary
    if ( !CellIsFluid( interiorCellIndex, mask, mesh ) ) {
        throw std::runtime_error( "Invalid immersed boundary geometry and mesh specification. Not enough fluid cells between solid and domain boundaries!" );
    }

    // Distance from cell center to immersed boundary along this coordinate direction
    fVector3 queryPointCoords( mesh.cellCenters[X](cellIndex[X]),
                               mesh.cellCenters[Y](cellIndex[Y]),
                               mesh.cellCenters[Z](cellIndex[Z]) );
    fVector3 rayDirection( 0, 0, 0 );
    rayDirection[ axis ] = static_cast<floatType>( directionIndex );
    sourceTermData.ibDistance2 = GetBoundaryDistance2(geometry, queryPointCoords, rayDirection);
    floatType ibDistance = sqrt( sourceTermData.ibDistance2 );


    // Face area vector
    sourceTermData.faceAreaComponent = static_cast<floatType>( directionIndex )   // Gives the correct sign
                                     * mesh.cellCenters[axis]( fidx );    


    // Extrapolation coefficients onto face between ghost cell and fluid cell. May use further points due to stability condition
    bool meetsGhostStabilityCondition =  ibDistance >= ( mesh.cellLengths[axis](cellIndex[axis]) / 2.0f );
    if ( meetsGhostStabilityCondition ) {

        floatType dxpOn2 = mesh.cellLengths[axis](cellIndex[axis]) / 2.0f ;
        sourceTermData.faceExtrapCoeff_p  = ( ibDistance - dxpOn2 ) / ( ibDistance );
        sourceTermData.faceExtrapCoeff_a  = 0.0f;
        sourceTermData.faceExtrapCoeff_ib = dxpOn2 / ibDistance;

    } else {

        floatType dxp = mesh.cellLengths[axis](cellIndex[axis]);
        floatType dxa = mesh.cellLengths[axis](interiorCellIndex[axis]);
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


    // Extrapolation coefficients
    floatType cellInteriorDistance = abs( mesh.cellCenters[axis](interiorCellIndex[axis]) - mesh.cellCenters[axis](cellIndex[axis]) ); 
    sourceTermData.ibExtrapFactor_p = ( cellInteriorDistance + ibDistance ) / cellInteriorDistance;
    sourceTermData.ibExtrapFactor_a = - ibDistance / cellInteriorDistance;


    // Far pressure ghost cell coefficients
    if        ( directionIndex == +1 ) {    // Ghost cell on Hi side

        floatType dxp      = mesh.cellLengths[axis](cellIndex[axis]),
                  dxe      = mesh.cellLengths[axis](ghostCellIndex[axis]),
                  lambdaw  = mesh.interpFactors[axis](interiorCellIndex[axis]),
                  lambdae  = mesh.interpFactors[axis](cellIndex[axis] + 1),
                  lambdaee = mesh.interpFactors[axis](cellIndex[axis] + 2),
                  le       = 1.0f /  mesh.cellCenterDiffInv[axis](cellIndex[axis] + 1);

        sourceTermData.farPressureCoeff_p = - (2 * dxe) / (lambdaee * le)
                                            - (dxe / dxp) * (1 - lambdae - lambdaw) / lambdaee
                                            + (1 - lambdae) / lambdaee;

        sourceTermData.farPressureCoeff_a = (dxe / dxp) * (1 - lambdaw) / lambdaee;

        sourceTermData.farPressureCoeff_g =   (2 * dxe) / (lambdaee * le)
                                            - (dxe / dxp) * (lambdae / lambdaee)
                                            - (1 - lambdae - lambdaee) / lambdaee;

    } else if ( directionIndex == -1) {     // Ghost cell on Lo side

        floatType dxp      = mesh.cellLengths[axis](cellIndex[axis]),
                  dxw      = mesh.cellLengths[axis](ghostCellIndex[axis]),
                  lambdae  = mesh.interpFactors[axis](interiorCellIndex[axis] + 1),
                  lambdaw  = mesh.interpFactors[axis](cellIndex[axis]),
                  lambdaww = mesh.interpFactors[axis](cellIndex[axis] - 1),
                  lw       = 1.0f / mesh.cellCenterDiffInv[axis](cellIndex[axis]);

        sourceTermData.farPressureCoeff_p = - (2 * dxw) / ((1 - lambdaww) * lw)
                                            + (dxw / dxp) * (1 - lambdae - lambdaw) / (1 - lambdaww)
                                            + lambdaw / (1 - lambdaww);

        sourceTermData.farPressureCoeff_a = (dxw / dxp) * lambdae / (1 - lambdaww);

        sourceTermData.farPressureCoeff_g =   (2 * dxw) / ((1 - lambdaww) * lw)
                                            - (dxw / dxp) * (1 - lambdaw) / (1 - lambdaww)
                                            + (1 - lambdaw - lambdaww) / (1 - lambdaww);

    }

}



void SetVelocityFluxCorrectionCoefficient( IBData &ibData )
{

    // The denominator must be calculated seperately first
    floatType denominator = 0.0f;
    for ( auto &ibCell : ibData.ibCells ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            denominator += std::pow( abs( sourceTermData.faceAreaComponent) , 2.0f )
                         * sourceTermData.ibDistance2;

        }
    }

    // Now the coefficient for each cell
    for ( auto &ibCell : ibData.ibCells ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

           sourceTermData.velocityFluxCorrectionCoeff = ( sourceTermData.ibDistance2 
                                                        * sourceTermData.faceAreaComponent 
                                                        ) / denominator;

        }
    }

}



IBData ConstructIBData( const Polyhedron &geometry,
                        const Mesh &mesh )
{

    using enum Axis::ENUMDATA;

    IBData ibData;
    
    ibData.mask = CreateCellMask( geometry, mesh );
    auto &mask = ibData.mask;

    // Iterate through all cells 
    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                if ( static_cast<intType>( mask(i, j, k) ) == CellType::Solid )
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
                    if ( !atHiBoundary && static_cast<intType>( mask(hiSideCellIndex) ) == CellType::Solid ) {
                        CheckIBCellPtr();
                        AddIBDataForDirection( *ibCellPtr, axis, +1, mask, mesh, geometry );
                    }

                    // Solid on lo side
                    TensorIndex3D loSideCellIndex = cellIndex;
                    loSideCellIndex[axis] -= 1;
                    bool atLoBoundary = ( cellIndex[axis] == 0  );
                    if ( !atLoBoundary && static_cast<intType>( mask(loSideCellIndex) ) == CellType::Solid ) {
                        CheckIBCellPtr();
                        AddIBDataForDirection( *ibCellPtr, axis, -1, mask, mesh, geometry );
                    }

                } );

            }
        }
    }

    // Calculate the coefficients for the velocity flux correction
    SetVelocityFluxCorrectionCoefficient( ibData );

    return ibData;
}


}   // end anonymous namespace




IBData CreateImmersedBoundaryData( const InputData &inputData, 
                                   const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    IBData ibData;

    if ( !inputData.hasIBGeometry ) {
        ibData.mask = Tensor3D( mesh.nCells[X], mesh.nCells[Y], mesh.nCells[Z] ).setConstant( CellType::Fluid );
        return ibData;
    }

    Polyhedron P = MakeGeometry( inputData );
    ibData = ConstructIBData( P, mesh );

    return ibData;
}


}   // end namespace CFD