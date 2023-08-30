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

#include <cmath>

namespace CFD
{

namespace 
{

// Find closes point on surface of polyhedron to the query point
fVector3 ClosestBoundaryPoint( const Polyhedron &polyhedron, 
                               const fVector3 &queryPoint) 
{
    using Point        = Polyhedron::Point_3;
    using Primitive    = CGAL::AABB_face_graph_triangle_primitive<Polyhedron>;
    using Traits       = CGAL::AABB_traits<CGAL_Kernel, Primitive>;
    using Tree         = CGAL::AABB_tree<Traits>;

    // Construct AABB tree with a KdTree
    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.accelerate_distance_queries();

    Point closest_point = tree.closest_point( Point( queryPoint(0), queryPoint(1), queryPoint(2)) );

    return { closest_point[0], closest_point[1], closest_point[2] };
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
// diagonal length of the cells surrounging the ghost cell
floatType GetImagePointDistance( const TensorIndex3D &ghostPointIndex, 
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

                TensorIndex3D neighbourIndex = ghostPointIndex;
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



IBData ConstructIBData( const Polyhedron &geometry,
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
                fVector3 boundaryPoint = ClosestBoundaryPoint( geometry, ghostPoint );

                // Get the normal vector
                fVector3 normalVector = boundaryPoint - ghostPoint;

                // Distance from ghost point to boundary
                floatType ghostPointDistance = normalVector.norm();

                // Determine the distance of the image point from the boundary
                floatType imagePointDistance = GetImagePointDistance( ghostPointIndex, mesh );

                // Determine coordinates of image point
                fVector3 imagePoint = boundaryPoint + imagePointDistance * normalVector.normalized();

                // Create the field probe for linear interpolation
                FieldProbe fieldProbe( mesh, imagePoint.array() );

                // Coefficient for extrapolation from image point to ghost point
                floatType extrapCoeff = 1.0f - ( imagePointDistance + ghostPointDistance ) / imagePointDistance;

                // Fill up the ghost cell data struct
                ibData.ghostCells.push_back( { fieldProbe, ghostPointIndex, extrapCoeff } );

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
    CellIDTensor3D cellID =  TagCells( mesh, P);
    ibData = ConstructIBData( P, cellID, mesh );

    return ibData;
}


}   // end namespace CFD