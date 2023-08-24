#include "ImmersedBoundary.h"

#include "../Geometry/Geometry.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>


namespace CFD
{

namespace 
{

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


void TagSolidCells( CellIDTensor3D &cellID,
                    const Mesh &mesh, 
                    const Polyhedron &geometry )
{
    using enum Axis::ENUMDATA;

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                floatType xq = mesh.cellCenters[X](i),
                          yq = mesh.cellCenters[Y](j),
                          zq = mesh.cellCenters[Z](k);

                if ( PointInside( geometry, xq, yq, zq ) ) {
                    cellID(i, j, k) = CellType::Solid;
                } 

            }
        }
    }

}


void TagGhostCells( CellIDTensor3D &cellID,
                    const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                if ( cellID(i, j, k) != CellType::Solid ) {
                    continue;
                }

                bool fluidEast   = false,
                     fluidWest   = false,
                     fluidNorth  = false,
                     fluidSouth  = false,
                     fluidTop    = false,
                     fluidBottom = false;

                if ( i != 0 ) 
                    fluidWest = cellID(i-1, j  , k  )   == CellType::Fluid;

                if ( i != mesh.nCells[X]-1 ) 
                    fluidEast = cellID(i+1, j  , k  )   == CellType::Fluid;

                if ( j != 0 ) 
                    fluidSouth  = cellID(i  ,j-1, k  )  == CellType::Fluid;

                if ( j != mesh.nCells[Y]-1 ) 
                    fluidNorth = cellID(i  , j+1, k  )  == CellType::Fluid;

                if ( k != 0 )
                    fluidBottom = cellID(i  , j  , k-1) == CellType::Fluid;

                if ( k != mesh.nCells[Z]-1 ) 
                    fluidTop = cellID(i  , j  , k+1)    == CellType::Fluid;


                bool hasNeighbouringFluidCell = fluidEast || fluidWest || fluidNorth || fluidSouth || fluidTop || fluidBottom;
                if ( hasNeighbouringFluidCell ) {
                    cellID(i, j, k) = CellType::Ghost;
                }


            }
        }
    }

}

}   // end anonymous namespace




CellIDTensor3D TagCells( const Mesh &mesh,
                         const Polyhedron &geometry )
{
    using enum Axis::ENUMDATA;

    CellIDTensor3D cellID( mesh.nCells[X], mesh.nCells[Y], mesh.nCells[Z] );
    cellID.setConstant( CellType::Fluid );

    TagSolidCells( cellID, mesh, geometry );
    
    TagGhostCells( cellID, mesh );

    return cellID;
}



}