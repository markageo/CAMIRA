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
        return 0;

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



// Marks faces as either being in fluid or solid region
EnumVector<Axis, CellIDTensor3D> TagCellFaces( const Mesh &mesh, 
                                               const Polyhedron &geometry )
{
    using enum Axis::ENUMDATA;

    EnumVector<Axis, CellIDTensor3D> faceIDs( { CellIDTensor3D(mesh.nFacesNormal[X](X), mesh.nFacesNormal[X](Y), mesh.nFacesNormal[X](Z)),
                                                CellIDTensor3D(mesh.nFacesNormal[Y](X), mesh.nFacesNormal[Y](Y), mesh.nFacesNormal[Y](Z)),
                                                CellIDTensor3D(mesh.nFacesNormal[Z](X), mesh.nFacesNormal[Z](Y), mesh.nFacesNormal[Z](Z)) } );

    // Iterate though all cell faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA faceNormal) {
        
        for (intType k = 0; k != mesh.nFacesNormal[faceNormal](Z); k++) {
            for (intType j = 0; j != mesh.nFacesNormal[faceNormal](Y); j++) {
                for (intType i = 0; i != mesh.nFacesNormal[faceNormal](X); i++) {

                    floatType xq = mesh.cellFaces[X](i),
                              yq = mesh.cellFaces[Y](j),
                              zq = mesh.cellFaces[Z](k);

                    if ( PointInside( geometry, xq, yq, zq ) ) {
                        faceIDs[faceNormal](i, j, k) = CellType::Solid;
                    } else {
                        faceIDs[faceNormal](i, j, k) = CellType::Fluid;
                    }

                }
            }
        }

    } );

}



// Create a mask array for cell centers
Tensor3D CreateCellMask( const EnumVector<Axis, CellIDTensor3D> &cellFaceIDs, 
                         const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    Tensor3D mask = Tensor3D( mesh.nCells[X], mesh.nCells[Y], mesh.nCells[Z] ).setConstant( 1.0f );

    // Iterate through each cell
    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                // Mask it if the current cell has any faces which are tagged as solid
                for ( intType a = 0; a != Axis::count; a++ ) {
                    Axis::ENUMDATA faceNormal = static_cast<Axis::ENUMDATA>( a );

                    TensorIndex3D loFaceIndex = {i, j, k},
                                  hiFaceIndex = {i, j, k};
                    hiFaceIndex[faceNormal] += 1;

                    if ( cellFaceIDs[faceNormal](hiFaceIndex) == CellType::Solid || cellFaceIDs[faceNormal](loFaceIndex) == CellType::Solid ) {
                        mask(i, j, k) = 0.0f;
                        break;
                    }

                }

            }
        }
    }

}



// A wrapper for the find_if function which returns iterator to cell element in IBData object that corresponds to the given index. Returns
// one past end iterator if does not exist.
std::vector<IBCell>::iterator FindIBCell( IBData ibData, 
                                          const TensorIndex3D forcedCellIndex )
{
    std::vector<IBCell>::iterator ibCellsIterator = std::find_if( ibData.IBCells.begin(), ibData.IBCells.end(), 
                                                                  [&](IBCell ibCell) { return ( ibCell.cellIndex == forcedCellIndex ); } );

    return ibCellsIterator;
}



// Sets face data information for a particular boundary cell that will contain a forcing term.
void AddIBDataFace( IBData &ibData, 
                    const TensorIndex3D &positiveSideFaceIndex, 
                    const TensorIndex3D &forcedCellIndex, 
                    const intType directionIndex, 
                    const Axis::ENUMDATA faceNormal, 
                    const Mesh &mesh, 
                    const Polyhedron &geometry )
{
    using enum Axis::ENUMDATA;

    // Check if cell is already inside the IBCells vector
    std::vector<IBCell>::iterator ibCellsIterator = FindIBCell( ibData, forcedCellIndex );

    if ( ibCellsIterator == ibData.IBCells.end() ) {
        ibData.IBCells.emplace_back();
        ibCellsIterator = ibData.IBCells.end() - 1;
        ibCellsIterator->cellIndex = forcedCellIndex;
    } 


    // Create new element for this face
    ibCellsIterator->facesData.emplace_back();
    IBCell::FaceData &faceData = ibCellsIterator->facesData.back();

    // Face location data relative to forced cell
    faceData.faceNormal = faceNormal;

    if ( directionIndex == 1 ) {
        faceData.faceIndexOffset = 0;           // Zero due to the staggering of face index relative to cell index
    } else if ( directionIndex == -1 ) {
        faceData.faceIndexOffset = 1;
    }


    // Index of adject face for extrapolation
    faceData.adjacentCellIndex = forcedCellIndex;
    faceData.adjacentCellIndex[faceNormal] += directionIndex;

    // Distance of forced face and cell to boundary
    fVector3 queryFacePointCoords( mesh.cellFaces[X](positiveSideFaceIndex[X]),
                                   mesh.cellFaces[Y](positiveSideFaceIndex[Y]),
                                   mesh.cellFaces[Z](positiveSideFaceIndex[Z]) );
    fVector3 rayDirection( 0, 0, 0 );
    rayDirection[ faceNormal ] = - directionIndex;

    floatType ibFaceDistance = GetBoundaryDistance( geometry, queryFacePointCoords, rayDirection );
    floatType ibCellDistance = ibFaceDistance + mesh.cellLengths[faceNormal]( forcedCellIndex[faceNormal] ) / 2.0;


    // Interpolation coefficients
    floatType cellHalfWidth  = mesh.cellLengths[ faceNormal ]( forcedCellIndex[faceNormal]);
    faceData.interpCoeffCell = ibFaceDistance / ( cellHalfWidth + ibFaceDistance );
    faceData.interpCoeffIB   = cellHalfWidth  / ( cellHalfWidth + ibFaceDistance );

    // Extrapolation coefficients
    floatType cellCenterSpacing = mesh.cellLengths[faceNormal]( forcedCellIndex[faceNormal] ) / 2.0 
                                + mesh.cellLengths[faceNormal]( faceData.adjacentCellIndex[faceNormal] ) / 2.0;
    faceData.extrapFactor_p = ( cellCenterSpacing + ibCellDistance ) / cellCenterSpacing;
    faceData.extrapFactor_a =  - ibCellDistance / cellCenterSpacing;
}



IBData ConstructIBData( const Polyhedron &geometry,
                        const Mesh &mesh )
{

    using enum Axis::ENUMDATA;

    IBData ibData;
    
    EnumVector<Axis, CellIDTensor3D> cellFaceIDs = TagCellFaces( mesh, geometry );

    ibData.mask = CreateCellMask( cellFaceIDs, mesh );

    // Iterate though all cell faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA faceNormal) {
        
        
        for (intType k = 0; k != mesh.nFacesNormal[faceNormal](Z); k++) {
            for (intType j = 0; j != mesh.nFacesNormal[faceNormal](Y); j++) {
                for (intType i = 0; i != mesh.nFacesNormal[faceNormal](X); i++) {

                    if ( cellFaceIDs[faceNormal](i, j, k) == CellType::Fluid )
                        continue;


                    // Forced face on positive side
                    TensorIndex3D positiveSideFaceIndex = {i, j, k};
                    positiveSideFaceIndex[faceNormal] += 1;

                    TensorIndex3D positiveForcedCellIndex = positiveSideFaceIndex;

                    intType positiveDirectionIndex = 1; 

                    if ( cellFaceIDs[faceNormal](positiveSideFaceIndex) == CellType::Fluid ) {
                        AddIBDataFace( ibData, positiveSideFaceIndex, positiveForcedCellIndex, positiveDirectionIndex, faceNormal, mesh, geometry );
                    }


                    // Forced face on negative side
                    TensorIndex3D negativeSideFaceIndex = {i, j, k};
                    negativeSideFaceIndex[faceNormal] -= 1;

                    TensorIndex3D negativeForcedCellIndex = negativeSideFaceIndex;
                    negativeForcedCellIndex[faceNormal] -= 1;

                    intType negativeDirectionIndex = -1; 

                    if ( cellFaceIDs[faceNormal](negativeSideFaceIndex) == CellType::Fluid ) {
                        AddIBDataFace( ibData, negativeSideFaceIndex, negativeForcedCellIndex, negativeDirectionIndex, faceNormal, mesh, geometry );
                    }

                }
            }
        }


    } );

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