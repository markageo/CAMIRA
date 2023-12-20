#include "ImmersedBoundary.h"
#include "../Tools/FVLookups.h"

#include "../IO/ArrayIO.h"

#include <vector>
#include <utility>
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



// Returns the distance to the nearest intersection from the given point in the direction of the given ray.
// https://stackoverflow.com/questions/69953358/to-calculate-intersections-between-a-ray-and-a-mesh
floatType GetBoundaryDistance( const Polyhedron &polyhedron,
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

    Vector rayOrientation( rayDirection(0), rayDirection(1), rayDirection(2) );
    Point  rayOrigin( queryPointCoords(0), queryPointCoords(1), queryPointCoords(2) );
    Ray    ray( rayOrigin, rayOrientation );
    std::vector<Ray_intersection> intersections;
    tree.all_intersections( ray, std::back_inserter( intersections ) );

    if ( intersections.empty() )
        return 0.0f;

    floatType closestDistance = 0.0;
    for ( auto it = intersections.begin(); it != intersections.end(); it++ ) {
        auto intersection = (*it)->first;

        if ( Point* p = boost::get<Point>(&intersection) ) {

            floatType distance = CGAL::squared_distance( *p, rayOrigin );
            if ( it == intersections.begin() || distance < closestDistance ) {
                closestDistance = distance;
            } 

        }

    }
    return closestDistance;
}


// Create a mask array for cell centers
Tensor3D CreateCellMask( const Polyhedron &geometry, 
                         const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    Tensor3D mask = Tensor3D( mesh.nCells[X], mesh.nCells[Y], mesh.nCells[Z] ).setConstant( 1.0f );

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



// Sets data for a particular source term in a particular direciton for a particular cell
void AddIBDataForDirection( IBCell &ibCell, 
                            const Axis::ENUMDATA axis,
                            const intType directionIndex,
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
    sourceTermData.adjacentCellIndex  = interiorCellIndex;


    // Distance from cell to immersed boundary along this coordinate direction
    fVector3 queryPointCoords( mesh.cellCenters[X](cellIndex[X]),
                               mesh.cellCenters[Y](cellIndex[Y]),
                               mesh.cellCenters[Z](cellIndex[Z]) );
    fVector3 rayDirection( 0, 0, 0 );
    rayDirection[ axis ] = static_cast<floatType>( directionIndex );
    floatType ibDistance   = GetBoundaryDistance(geometry, queryPointCoords, rayDirection);


    // Interpolation coefficients onto ghost cell
    floatType cellGhostDistance = abs( mesh.cellCenters[axis](ghostCellIndex[axis]) - mesh.cellCenters[axis](cellIndex[axis]) );
    sourceTermData.cellInterpCoeff_p  = 1 - cellGhostDistance / ibDistance;
    sourceTermData.cellInterpCoeff_ib = cellGhostDistance / ibDistance;


    // Interpolation coefficients onto cell face between ghost cell
    if        ( directionIndex == +1 ) {    // Face on Hi side

        intType fidx = cellIndex[axis] + sourceTermData.faceDirectionIndex;
        sourceTermData.faceInterpCoeff_p  = 1 - mesh.interpFactors[axis](fidx);
        sourceTermData.faceInterpCoeff_ib = mesh.interpFactors[axis](fidx);

    } else if ( directionIndex == -1 ) {   // Face on Lo side

        intType fidx = cellIndex[axis] + sourceTermData.faceDirectionIndex;
        sourceTermData.faceInterpCoeff_p  = mesh.interpFactors[axis](fidx);
        sourceTermData.faceInterpCoeff_ib = 1 - mesh.interpFactors[axis](fidx);

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
                  le       = mesh.cellCenterDiffInv[axis](cellIndex[axis] + 1);

        sourceTermData.farPressureCoeff_p = - (2 * dxe) / (lambdaee * le)
                                            - (dxe / dxp) * (1 - lambdae - lambdaw) / lambdaee
                                            + (1 - lambdae) / lambdaee;

        sourceTermData.farPressureCoeff_a = (dxe / dxp) * (1 - lambdaw) / lambdaee;

        sourceTermData.farPressureCoeff_g =   (2 * dxe) / (lambdaee * le)
                                            - (dxe / dxp) * lambdae / lambdaee
                                            - (1 - lambdae - lambdaee) / lambdaee;

    } else if ( directionIndex == -1) {     // Ghost cell on Lo side

        floatType dxp      = mesh.cellLengths[axis](cellIndex[axis]),
                  dxw      = mesh.cellLengths[axis](ghostCellIndex[axis]),
                  lambdae  = mesh.interpFactors[axis](interiorCellIndex[axis]),
                  lambdaw  = mesh.interpFactors[axis](cellIndex[axis] - 0),
                  lambdaww = mesh.interpFactors[axis](cellIndex[axis] - 1),
                  lw       = mesh.cellCenterDiffInv[axis](cellIndex[axis] - 0);

        sourceTermData.farPressureCoeff_p = - (2 * dxw) / ((1 - lambdaww) * lw)
                                            + (dxw / dxp) * (1 - lambdae - lambdaw) / (1 - lambdaww)
                                            + lambdaw / (1 - lambdaww);

        sourceTermData.farPressureCoeff_a = (dxw / dxp) * lambdaw / (1 - lambdaww);

        sourceTermData.farPressureCoeff_g =   (2 * dxw) / ((1 - lambdaww) * lw)
                                            - (dxw / dxp) * (1 - lambdaw) / (1 - lambdaww)
                                            + (1 - lambdaw - lambdaww) / (1 - lambdaww);

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

                if ( static_cast<intType>( mask(i, j, k) ) != CellType::Fluid )
                    continue;

                TensorIndex3D cellIndex{i, j, k};
                bool cellHasBeenAdded = false;

                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                    // Solid on hi side
                    TensorIndex3D hiSideCellIndex = cellIndex;
                    hiSideCellIndex[axis] += 1;
                    bool atHiBoundary = ( cellIndex[axis] == mesh.nCells[axis]-1  );
                    if ( !atHiBoundary && static_cast<intType>( mask(hiSideCellIndex) ) == CellType::Solid ) {
                        
                        if ( !cellHasBeenAdded ) {
                            ibData.ibCells.emplace_back();
                            ibData.ibCells.back().cellIndex = cellIndex;
                            cellHasBeenAdded = true;
                        }

                        IBCell &ibCell = ibData.ibCells.back();
                        AddIBDataForDirection( ibCell, axis, +1, mesh, geometry );

                    }

                    // Solid on lo side
                    TensorIndex3D loSideCellIndex = cellIndex;
                    loSideCellIndex[axis] -= 1;
                    bool atLoBoundary = ( cellIndex[axis] == 0  );
                    if ( !atLoBoundary && static_cast<intType>( mask(loSideCellIndex) ) == CellType::Solid ) {
                        
                        if ( !cellHasBeenAdded ) {
                            ibData.ibCells.emplace_back();
                            ibData.ibCells.back().cellIndex = cellIndex;
                            cellHasBeenAdded = true;
                        }

                        IBCell &ibCell = ibData.ibCells.back();
                        AddIBDataForDirection( ibCell, axis, -1, mesh, geometry );

                    }

                } );

            }
        }
    }

    return ibData;
}


}   // end anonymous namespace




IBData CreateImmersedBoundaryData( const InputData &inputData, 
                                   const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    IBData ibData;

    if ( !inputData.hasIBGeometry ) {
        ibData.mask = Tensor3D( mesh.nCells[X], mesh.nCells[Y], mesh.nCells[Z] ).setConstant( 1.0f );
        return ibData;
    }

    Polyhedron P = MakeGeometry( inputData );
    ibData = ConstructIBData( P, mesh );

    return ibData;
}


}   // end namespace CFD