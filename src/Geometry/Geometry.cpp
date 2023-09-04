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

#include <CGAL/boost/graph/helpers.h>

#include <cmath>

namespace CFD
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
            floatType angleRadians = blockData.rotation(axis) * M_PI / 180.0f;
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


}   // end anonymous namespace



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
    floatType radius2 = std::pow( solidSphereData.diameter / 2.0f, 2 );
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
    Surface_3 surface( sphere_function,                     // pointer to function
                       Sphere_3(centerPosition, 2.0f * radius2) ); // bounding sphere

    // Meshing criteria
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria( 30.0f,  // angular bound
                                                        0.1f,   // radius bound
                                                        0.1f);  // distance bound

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
                                   Complete Geometry
\*-------------------------------------------------------------------------------------*/


Polyhedron MakeGeometry( const InputData &inputData )
{
    std::vector< Polyhedron > geometryPolyhedra;
    AddBlocks( geometryPolyhedra, inputData );
    AddSpheres( geometryPolyhedra, inputData );

    Polyhedron P;
    for ( Polyhedron & polyhedron : geometryPolyhedra ) {
        CGAL::Polygon_mesh_processing::corefine_and_compute_union( polyhedron, P, P );
    }

    return P;
}


}