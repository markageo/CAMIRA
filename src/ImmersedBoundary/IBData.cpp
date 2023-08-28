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



// Add fluid points and their distance squared to the foundPoints vector which are within a shell around the LoCornerIndex point in the mesh.
void SearchShell( std::vector< std::pair<TensorIndex3D, floatType> > &foundPoints,
                  const intType searchShellSize,
                  const fVector3 &queryPoint,
                  const TensorIndex3D &LoCornerIndex,
                  const CellIDTensor3D &cellID,
                  const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    auto CellCoordinate = [&] ( const TensorIndex3D &index ) -> fVector3
    { return { mesh.cellCenters[X](index[X]), mesh.cellCenters[Y](index[Y]), mesh.cellCenters[Z](index[Z]) }; };


    TensorIndex3D shellLowerBound, shellUpperBound;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        shellLowerBound[axis] = LoCornerIndex[axis] - (searchShellSize-1);
        shellUpperBound[axis] = LoCornerIndex[axis] + searchShellSize;
    } );

    for ( intType k = shellLowerBound[Z]; k <= shellUpperBound[Z]; k++ ) {
        for ( intType j = shellLowerBound[Y]; j <= shellUpperBound[Y]; j++ ) {
            for ( intType i = shellLowerBound[X]; i <= shellUpperBound[X]; i++ ) {

                bool OnShellBoundaryZ = ( k == shellLowerBound[Z] ) || ( k == shellUpperBound[Z] ),
                     OnShellBoundaryY = ( j == shellLowerBound[Y] ) || ( j == shellUpperBound[Y] ),
                     OnShellBoundaryX = ( i == shellLowerBound[X] ) || ( i == shellUpperBound[X] );
                
                if ( !(OnShellBoundaryX || OnShellBoundaryY || OnShellBoundaryZ) ) {
                    continue;
                }

                TensorIndex3D shellPointIndex = {i, j, k};

                if ( OutOfBounds( shellPointIndex, mesh ) ) {
                    continue;
                }

                if ( cellID( shellPointIndex ) != CellType::Fluid ) {
                    continue;
                }

                // Position vector of shell point
                fVector3 shellPoint = CellCoordinate( shellPointIndex );

                // Distance squared from query point to shell point
                floatType shellPointDistance2 = ( shellPoint - queryPoint ).squaredNorm();

                // Add this distance and index to the vectors
                foundPoints.push_back( std::make_pair( shellPointIndex, shellPointDistance2 ) );

            }
        }
    }

}



// Find Lo index of nodal box which surrounds the query point
// If the query point is between a node and the boundary, then the Lo index will be 
// outside the domain bounds
TensorIndex3D FindLoCornerIndex( const fVector3 &queryPoint,
                                 const TensorIndex3D &startingPointIndex, 
                                 const Mesh &mesh )
{
    TensorIndex3D LoCornerIndex = startingPointIndex;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        bool LoCornerIsOnLoSide = false,
             HiCornerIsOnHiSide = false;

        while ( !LoCornerIsOnLoSide && !HiCornerIsOnHiSide ) {

            bool axisOutOfBounds =  ( LoCornerIndex[axis] < 0 ) || ( LoCornerIndex[axis] > mesh.nCells[axis]-1 );
            if ( axisOutOfBounds ) {
                break;
            }

            LoCornerIsOnLoSide = mesh.cellCenters[axis]( LoCornerIndex[axis] ) <= queryPoint(axis);
            if ( !LoCornerIsOnLoSide ) {
                LoCornerIndex[axis]--;
            }

            HiCornerIsOnHiSide = mesh.cellCenters[axis]( LoCornerIndex[axis] + 1 ) > queryPoint(axis);
            if ( !HiCornerIsOnHiSide ) {
                LoCornerIndex[axis]++;
            }

        }

    } );

    return LoCornerIndex;
}




// Return index of nearest fluid cells to query point
std::array< TensorIndex3D, IBGhostCell::numFluidInterpPoints > NearestFluidCellIndices( const fVector3 &queryPoint, 
                                                                                   const TensorIndex3D &startingPointIndex,
                                                                                   const CellIDTensor3D &cellID,
                                                                                   const Mesh &mesh )
{
    constexpr intType numPointsToFind = IBGhostCell::numFluidInterpPoints;
    std::array< TensorIndex3D, IBGhostCell::numFluidInterpPoints > nearestPointIndices;

    TensorIndex3D LoCornerIndex = FindLoCornerIndex( queryPoint, startingPointIndex, mesh );    // This could be out of bounds

    std::vector< std::pair<TensorIndex3D, floatType> > foundPoints;

    intType maxSearchShellSize = 3;
    bool foundEnoughPoints = false;
    for ( intType searchShellSize = 1; searchShellSize < maxSearchShellSize; searchShellSize++ ) {

        SearchShell( foundPoints, searchShellSize, queryPoint, LoCornerIndex, cellID, mesh );

        if ( foundPoints.size() >= numPointsToFind ) {

            foundEnoughPoints = true;

            // Sort the points
            std::sort( foundPoints.begin(), foundPoints.end(), [](const std::pair<TensorIndex3D, floatType> &a, const std::pair<TensorIndex3D, floatType> &b) {
                return a.second < b.second;
            } );

            // Take the nearest ones
            for ( intType i = 0; i != numPointsToFind; i++ ) {
                nearestPointIndices[i] = foundPoints[i].first;
            }

            break;
        }
    }

    if ( !foundEnoughPoints ) {
        throw std::runtime_error("Immersed boundary geometry error! Not enough fluid interpolation cells found.");
    }

    return nearestPointIndices;
}



// Construct the matrix that is used to determine the polynomial coefficients for interpolating the image point
IBGhostCell::InterpMatrix GetPointsMatrixInv( const std::array< TensorIndex3D, IBGhostCell::numFluidInterpPoints > &fluidCellIndices,
                                              const fVector3 &boundaryPoint, 
                                              const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using MatrixRow = Eigen::Matrix<floatType, 1, IBGhostCell::numInterpPoints>;
    IBGhostCell::InterpMatrix pointsMatrixInv;

    // First row is the boundary point
    pointsMatrixInv.row(0) = MatrixRow{1.0f, boundaryPoint(0), boundaryPoint(1), boundaryPoint(2)};

    for ( intType i = 1; i != IBGhostCell::numInterpPoints; i++ ) {

        floatType xf = mesh.cellCenters[X]( fluidCellIndices[i-1][X] ),
                  yf = mesh.cellCenters[Y]( fluidCellIndices[i-1][Y] ),
                  zf = mesh.cellCenters[Z]( fluidCellIndices[i-1][Z] );

         pointsMatrixInv.row(i) = MatrixRow{1.0f, xf, yf, zf};
    }

    return pointsMatrixInv.inverse();
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

                // Determine coordinates of image point
                fVector3 imagePoint = boundaryPoint + normalVector;

                // Find nearest fluid cells for interpolation
                std::array< TensorIndex3D, IBGhostCell::numFluidInterpPoints > fluidCellIndices = NearestFluidCellIndices( boundaryPoint, ghostPointIndex, cellID, mesh );

                // Form the interpolation coefficient matrix
                IBGhostCell::InterpMatrix pointsMatrixInv = GetPointsMatrixInv( fluidCellIndices, boundaryPoint, mesh );

                // Fill up the ghost cell data struct
                ibData.ghostCells.push_back( { fluidCellIndices, ghostPointIndex, imagePoint, pointsMatrixInv } );

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