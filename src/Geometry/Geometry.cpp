#include "Geometry.h"

#include <Eigen/Dense>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Aff_transformation_3.h>

#include <boost/property_map/property_map.hpp>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/IO/STL.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace CAMIRA
{

namespace 
{

/*-------------------------------------------------------------------------------------*\
                                        Blocks
\*-------------------------------------------------------------------------------------*/

std::array< fVector3, 8> GetBlockVertices( const InputData::SolidBlockData &blockData )
{
    using enum Axis::ENUMDATA;
    std::array< fVector3, 8 > blockVertices;

    const floatType dx = blockData.dimensions(X) / 2.0f,
                    dy = blockData.dimensions(Y) / 2.0f,
                    dz = blockData.dimensions(Z) / 2.0f;

    // Block vertices, order is important!
    blockVertices[0] = blockData.centerPosition + fArray3( { dx, -dy, -dz} );   // (1, 0, 0)
    blockVertices[1] = blockData.centerPosition + fArray3( {-dx, -dy,  dz} );   // (0, 0, 1)
    blockVertices[2] = blockData.centerPosition + fArray3( {-dx, -dy, -dz} );   // (0, 0, 0)
    blockVertices[3] = blockData.centerPosition + fArray3( {-dx,  dy, -dz} );   // (1, 0, 1)
    
    blockVertices[4] = blockData.centerPosition + fArray3( { dx, -dy,  dz} );   // (1, 0, 1)
    blockVertices[5] = blockData.centerPosition + fArray3( {-dx,  dy,  dz} );   // (0, 1, 1)
    blockVertices[6] = blockData.centerPosition + fArray3( { dx,  dy, -dz} );   // (1, 1, 0)
    blockVertices[7] = blockData.centerPosition + fArray3( { dx,  dy,  dz} );   // (1, 1, 1)

    // Unitvectors
    std::array< fVector3, 3 > unitVectors = { fVector3(1.0f, 0.0f, 0.0f),
                                              fVector3(0.0f, 1.0f, 0.0f), 
                                              fVector3(0.0f, 0.0f, 1.0f) };
    fVector3 blockCenterVector = blockData.centerPosition.matrix(); 

    // Rotate each point about the block center
    for ( fVector3 &vertexPoint : blockVertices ) {

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

            // Rotate about origin
            floatType angleRadians = blockData.rotation(axis) * static_cast<floatType>( M_PI ) / 180.0f;
            Eigen::AngleAxis<floatType> R( angleRadians, unitVectors[axis] );
            vertexPoint = R * vertexPoint;

            // Translate
            fVector3 translationVector = blockCenterVector - R * blockCenterVector;
            Eigen::Translation<floatType, 3> T( translationVector );
            vertexPoint = T * vertexPoint;

        } );

    }
    
    return blockVertices;
}



// Breaks CGAL Polyhedron block created through Euler operations into tetrahedrons
void TriangulateBlockPolyhedron( Polyhedron &P)
{
    using Halfedge_handle = Polyhedron::Halfedge_handle; 

    // Break this up into tetrahedrons by splitting facets
    
    // First corner
    Halfedge_handle h0 = P.edges_begin();
    P.split_facet( h0, h0->next()->next() );

    Halfedge_handle h1 = h0->opposite()->prev();
    P.split_facet( h1, h1->next()->next() );

    Halfedge_handle h2 = h1->opposite()->prev();
    P.split_facet( h2, h2->next()->next() );


    // Opposite corner
    Halfedge_handle h3 = h0->next()->next()->opposite()->next()->opposite()->prev()->prev();
    P.split_facet( h3, h3->next()->next() );

    Halfedge_handle h4 = h0->opposite()->next()->opposite()->prev()->opposite()->prev();
    P.split_facet( h4, h4->next()->next() );

    // Across
    Halfedge_handle h5 = h0->opposite()->next()->opposite()->prev();
    P.split_facet( h5, h5->next()->next() );

}



// Makes block given vertices using Euler operations. Vertex order is must be correct.
Polyhedron MakeBlockPolyhedron( const InputData::SolidBlockData &solidBlockData )
{
    using enum Axis::ENUMDATA;
    using Point           = Polyhedron::Point_3;
    using Halfedge_handle = Polyhedron::Halfedge_handle; 

    // Get block vertices
    std::array< fVector3 , 8> userBlockVertices = GetBlockVertices( solidBlockData );

    // Convert the input points
    constexpr size_t nVertices = 8;
    std::array< Point, nVertices > blockVertices;
    for ( size_t i = 0; i != nVertices; i++ ) {
        blockVertices[i] = Point( userBlockVertices[i](X),
                                  userBlockVertices[i](Y),
                                  userBlockVertices[i](Z) );
    }

     // See polyhedron_prog_cube.cpp from CGAL library examples
    Polyhedron P;
    Halfedge_handle h = P.make_tetrahedron( blockVertices[0],
                                            blockVertices[1],
                                            blockVertices[2],
                                            blockVertices[3] );

    Halfedge_handle g = h->next()->opposite()->next();      
    P.split_edge( h->next());
    P.split_edge( g->next());
    P.split_edge( g);                                              
    h->next()->vertex()->point()     = blockVertices[4];
    g->next()->vertex()->point()     = blockVertices[5];
    g->opposite()->vertex()->point() = blockVertices[6];      
    Halfedge_handle f = P.split_facet( g->next(),
                                       g->next()->next()->next());
    Halfedge_handle e = P.split_edge( f);
    e->vertex()->point() = blockVertices[7];                   
    P.split_facet( e, f->next()->next());                  


    // Turn it into tetrahedrons
    TriangulateBlockPolyhedron( P );

    CGAL_postcondition( P.is_valid());
    return P;
}



void AddBlocks( std::vector< Polyhedron > &geometryPolyhedra,
                const InputData &inputData )
{
    for ( const InputData::SolidBlockData &solidBlockData : inputData.solidBlocks ) {

        geometryPolyhedra.push_back( MakeBlockPolyhedron( solidBlockData ) );

    }
}



/*-------------------------------------------------------------------------------------*\
                                        Spheres
\*-------------------------------------------------------------------------------------*/


// Makes triangulated sphere surface
Polyhedron MakeSpherePolyhedron( const InputData::SolidSphereData &solidSphereData )
{
    using enum Axis::ENUMDATA;
    
    using Tr       = CGAL::Surface_mesh_default_triangulation_3;
    using C2t3     = CGAL::Complex_2_in_triangulation_3<Tr>;
    using GT       = Tr::Geom_traits;
    using Sphere_3 = GT::Sphere_3;
    using Point_3  = GT::Point_3;
    using FT       = GT::FT;

    typedef FT (*Function)(Point_3);
    using Surface_3 = CGAL::Implicit_surface_3<GT, Function>;
    using Surface_mesh = CGAL::Surface_mesh<Point_3>;

    // Implicit function for sphere surface
    floatType radius2 = std::pow( solidSphereData.diameter / 2.0f, 2.0f );
    Point_3 centerPosition = Point_3( solidSphereData.centerPosition(X), 
                                      solidSphereData.centerPosition(Y), 
                                      solidSphereData.centerPosition(Z) );
    auto sphere_function = [&] (Point_3 p) -> FT
    {
        const FT x2 = std::pow( ( p.x() - centerPosition.x() ), 2 ), 
                 y2 = std::pow( ( p.y() - centerPosition.y() ), 2 ),
                 z2 = std::pow( ( p.z() - centerPosition.z() ), 2 );
        return x2 + y2 + z2 - radius2;
    };

    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

    // defining the surface
    Surface_3 surface( sphere_function,                            // pointer to function
                       Sphere_3(centerPosition, 2.0f * radius2) ); // bounding sphere

    // Meshing criteria
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria( 5.0f,                       // lower bound on minimum angle of facets (degrees)
                                                        0.025f * sqrt( radius2 ),    // upper bound on radius of Delaunay balls
                                                        0.01f * sqrt( radius2 ) );  // upper bound on center-center facet distances

    // meshing surface
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
    Surface_mesh sm;
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

    // Convert to Polyhedron
    Polyhedron P;
    CGAL::copy_face_graph(sm, P);

    return P;
}


void AddSpheres( std::vector< Polyhedron > &geometryPolyhedra,
                 const InputData &inputData )
{
    for ( const InputData::SolidSphereData &solidSphereData : inputData.solidSpheres ) {

        geometryPolyhedra.push_back( MakeSpherePolyhedron( solidSphereData ) );

    }
}



/*-------------------------------------------------------------------------------------*\
                                       From file
\*-------------------------------------------------------------------------------------*/


void RotatePolyhedron( Polyhedron &P, 
                       const fArray3 &rotationAnglesDegrees )
{
    using Transformation = CGAL::Aff_transformation_3<CGAL_Kernel>;
    using enum Axis::ENUMDATA;

    // Convert to radians
    const fArray3 rotationAnglesRadians = rotationAnglesDegrees * M_PI / 180.0; 
    const floatType thetaX = rotationAnglesRadians(X),
                    thetaY = rotationAnglesRadians(Y),
                    thetaZ = rotationAnglesRadians(Z);  

    // Rotate X
    Transformation rotationX( 1.0, 0.0, 0.0,
                              0.0, std::cos( thetaX ), - std::sin( thetaX ),
                              0.0, std::sin( thetaX ),   std::cos( thetaX ) );

    // Rotate Y
    Transformation rotationY(  std::cos( thetaY ), 0.0, std::sin( thetaY ),
                               0.0, 1.0, 0.0, 
                              -std::sin( thetaY ), 0.0, std::cos( thetaY ) );

    // Rotate Z
    Transformation rotationZ( std::cos( thetaZ ), -std::sin( thetaZ ), 0.0,
                              std::sin( thetaZ ),  std::cos( thetaZ ), 0.0,
                              0.0, 0.0, 1.0 );

    // Apply transformation to all points                              
    for (auto p = P.points_begin(); p != P.points_end(); p++) {
        *p = rotationX.transform(*p);
        *p = rotationY.transform(*p);
        *p = rotationZ.transform(*p);
    }

}


void TransformPolyhedronFromUserToCodeCoordinates( Polyhedron &P,
                                                   const AxisTransformationMap &axisTransformation )
{
    using enum Axis::ENUMDATA;
    using Point = CGAL_Kernel::Point_3;

    for (auto p = P.points_begin(); p != P.points_end(); p++) {

        *p = Point( axisTransformation.CodeAxisReverseSign(X) * p->cartesian( axisTransformation.UserAxis(X) ),
                    axisTransformation.CodeAxisReverseSign(Y) * p->cartesian( axisTransformation.UserAxis(Y) ),
                    axisTransformation.CodeAxisReverseSign(Z) * p->cartesian( axisTransformation.UserAxis(Z) ) );

    }
}


void AddSTLFiles( std::vector< Polyhedron > &geometryPolyhedra,
                  const InputData &inputData,
                  const AxisTransformationMap &axisTransformation )
{

    for ( const InputData::STLGeometryData &stlGeometryData : inputData.stlGeometries ) {

        Polyhedron P;
        bool success = CGAL::IO::read_STL( stlGeometryData.filename, P );
        if ( !success ) {
            throw std::runtime_error( "Failed reading STL geometry file '" + stlGeometryData.filename + "'." );
        }

        // Transform the geometry to code coordinates
        TransformPolyhedronFromUserToCodeCoordinates( P, axisTransformation );

        // Rotate the geometry
        RotatePolyhedron( P, stlGeometryData.rotation );

        geometryPolyhedra.push_back( P );

    }

}



}   // end anonymous namespace


/*-------------------------------------------------------------------------------------*\
                                   Complete Geometry
\*-------------------------------------------------------------------------------------*/


Polyhedron MakeGeometry( const InputData &inputData,
                         const AxisTransformationMap &axisTransformation )
{
    std::vector< Polyhedron > geometryPolyhedra;
    AddBlocks( geometryPolyhedra, inputData );
    AddSpheres( geometryPolyhedra, inputData );
    AddSTLFiles( geometryPolyhedra, inputData, axisTransformation );

    Polyhedron P;
    for ( Polyhedron & polyhedron : geometryPolyhedra ) {
        CGAL::Polygon_mesh_processing::corefine_and_compute_union( polyhedron, P, P );
    }

    return P;
}



/*-------------------------------------------------------------------------------------*\
                               General Geometry Functions
\*-------------------------------------------------------------------------------------*/


Tree MakeAABBTree( const Polyhedron &polyhedron )
{
    // Construct AABB tree with a KdTree
    Tree tree(faces(polyhedron).first, faces(polyhedron).second, polyhedron);
    tree.accelerate_distance_queries();
    return tree;
}


std::vector<Polyhedron> SeparatePolyhedron( const Polyhedron &polyhedron )
{
    using face_descriptor     = boost::graph_traits<Polyhedron>::face_descriptor;
    using FaceComponentMap   = std::map<face_descriptor, std::size_t>;

    std::vector<Polyhedron> components;
    if ( polyhedron.empty() )
        return components;

    // FaceComponentMap holds references to the original polyhedron, a new one needs to be created for each sub polyhedron if keep_connected_components is to be used. 
    // Need to create the map for the original poly so we know the number of components
    FaceComponentMap faceComponentMap;
    std::size_t nComponents = CGAL::Polygon_mesh_processing::connected_components( polyhedron, 
                                                                                   boost::associative_property_map<FaceComponentMap>(faceComponentMap));

    // Make a new poly for each component
    for ( std::size_t i = 0; i != nComponents; i++ ) {

        Polyhedron subPolyhedron = polyhedron;
            
        // Make FaceComponentMap for this poly
        CGAL::Polygon_mesh_processing::connected_components( subPolyhedron, 
                                                             boost::associative_property_map<FaceComponentMap>(faceComponentMap));

        std::vector<std::size_t> componentsToKeep = {i};
        
        CGAL::Polygon_mesh_processing::keep_connected_components( subPolyhedron, 
                                                                  componentsToKeep, 
                                                                  boost::make_assoc_property_map(faceComponentMap));
        components.push_back( std::move( subPolyhedron ) );
    }
    
    return components;
}



// Determines if query point is inside given geometry tree
bool PointInside( const Tree &tree, 
                  const floatType xq, 
                  const floatType yq, 
                  const floatType zq )
{
    using Point        = Polyhedron::Point_3;
    using Point_inside = CGAL::Side_of_triangle_mesh<Polyhedron, CGAL_Kernel>;

    if ( tree.empty() )
        return false;

    // Initialize the point-in-polyhedron tester
    Point_inside inside_tester(tree);
    
    // Determine the side and return true if inside!
    return inside_tester( Point( xq, yq, zq ) ) == CGAL::ON_BOUNDED_SIDE;
}



// Returns the distance to the nearest intersection from the given point in the direction of the given ray.
floatType NearestRayIntersection( const Tree &tree,
                                  const fVector3 &queryPointCoords,
                                  const fVector3 &rayDirection )
{
    using Point            = Polyhedron::Point_3;
    using Vector           = CGAL_Kernel::Vector_3;
    using Ray              = CGAL_Kernel::Ray_3;    
    using Ray_intersection = boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type>;

    const Vector rayOrientation( rayDirection(0), rayDirection(1), rayDirection(2) );
    const Point  rayOrigin( queryPointCoords(0), queryPointCoords(1), queryPointCoords(2) );
    const Ray    ray( rayOrigin, rayOrientation );
    Ray_intersection intersection = tree.first_intersection( ray );

    if ( !intersection ) 
        return -1.0f;   // No intersection was found

    Point* intersectionPoint = boost::get<Point>( &(intersection->first) ); 
    floatType distance2 = static_cast<floatType>( CGAL::squared_distance( *intersectionPoint, rayOrigin ) );

    return sqrt( distance2 );
}



// Gets the nearest distance of point to geometry tree
floatType NearestDistance( const Tree &tree, 
                           const floatType xq, 
                           const floatType yq, 
                           const floatType zq  )
{
    using Point        = Polyhedron::Point_3;

    Point pointq(xq, yq, zq);
    CGAL_Kernel::FT distance2 = tree.squared_distance( pointq );

    return sqrt( static_cast<floatType>( distance2 ) );
}



}