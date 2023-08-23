#include "ImmersedBoundary.h"

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

// Find closes point on surface of polyhedron to the query point
fVector3 ClosestPoint( const Polyhedron &polyhedron, 
                       const fVector3 &queryPoint) 
{
    using Point        = Polyhedron::Point_3;
    using Primitive    = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
    using Traits       = CGAL::AABB_traits<CGAL_Kernel, Primitive>;
    using Tree         = CGAL::AABB_tree<Traits>;
    using Point_inside = CGAL::Side_of_triangle_mesh<Polyhedron, CGAL_Kernel>;

    // Construct AABB tree with a KdTree
    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.accelerate_distance_queries();

    Point closest_point = tree.closest_point( Point( queryPoint(0), queryPoint(1), queryPoint(2)) );
    return { closest_point[0], closest_point[1], closest_point[2] };
}



// Return index of nearest fluid cells to query point
std::array< TensorIndex3D, IBGhostCell::numInterpPoints > NearestFluidCellIndices( const fVector3 &queryPoint, 
                                                                                   const CellIDTensor3D &cellID,
                                                                                   const Mesh &mesh, 
                                                                                   const intType numPoints )
{
    std::array< TensorIndex3D, IBGhostCell::numInterpPoints > nearestPointIndices( numPoints );

    // TODO

    return nearestPointIndices;
}


IBGhostCell::InterpMatrix GetPointsMatrixInv( const std::array< TensorIndex3D, IBGhostCell::numInterpPoints > &fluidCellIndices,
                                              const fVector3 &boundaryPoint, 
                                              const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using MatrixRow = Eigen::Matrix<floatType, IBGhostCell::numInterpPoints, IBGhostCell::numInterpPoints>;
    IBGhostCell::InterpMatrix pointsMatrixInv;

    // First row is the boundary point
    pointsMatrixInv.row(0) = MatrixRow{1.0f, boundaryPoint(0), boundaryPoint(1), boundaryPoint(2)};

    for ( intType i = 1; i != IBGhostCell::numInterpPoints; i++ ) {

        floatType xf = mesh.cellCenters[X]( fluidCellIndices[i][X] ),
                  yf = mesh.cellCenters[Y]( fluidCellIndices[i][Y] ),
                  zf = mesh.cellCenters[Z]( fluidCellIndices[i][Z] );

         pointsMatrixInv.row(i) = MatrixRow{1.0f, xf, yf, zf};
    }

    return pointsMatrixInv.inverse();
}


}   // end anonymous namespace





IBData CreateImmersedBoundaryData( const Polyhedron &geometry,
                                   const CellIDTensor3D &cellID, 
                                   const Mesh &mesh )
{

    using enum Axis::ENUMDATA;

    IBData ibData;
    ibData.mask = Tensor3D( mesh.nCells[X], mesh.nCells[Y], mesh.nCells[Z] ).setZero();

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                if ( cellID( i, j, k ) != CellType::Ghost ) {

                    if ( cellID(i, j, k) == CellType::Fluid ) {
                        ibData.mask(i, j, k) = 1.0f;
                    }
                    continue;
                }

                fVector3 ghostPoint = { mesh.cellCenters[X](i), 
                                        mesh.cellCenters[Y](j),
                                        mesh.cellCenters[Z](k) };
                TensorIndex3D ghostPointIndex = {i, j, k};

                // Find coordinates of nearest immersed boundary point
                fVector3 boundaryPoint = ClosestPoint( geometry, ghostPoint );

                // Get the normal vector
                fVector3 normalVector = boundaryPoint - ghostPoint;

                // Determine coordinates of image point
                fVector3 imagePoint = boundaryPoint + normalVector;

                // Find nearest fluid cells for interpolation
                std::array< TensorIndex3D, IBGhostCell::numInterpPoints > fluidCellIndices = NearestFluidCellIndices( boundaryPoint, cellID, mesh, IBGhostCell::numInterpPoints );

                // Form the interpolation coefficient matrix
                IBGhostCell::InterpMatrix pointsMatrixInv = GetPointsMatrixInv( fluidCellIndices, boundaryPoint, mesh );

                // Fill up the ghost cell data struct
                ibData.ghostCells.push_back( { fluidCellIndices, ghostPointIndex, pointsMatrixInv } );

            }
        }
    }

    return ibData;

}


}