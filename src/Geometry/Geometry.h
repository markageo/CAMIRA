#ifndef CFD_GEOMETRY
#define CFD_GEOMETRY

#include "../Core/Types.h"
#include "../IO/InputProcessing.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>


namespace CFD
{

using CGAL_Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;    // This uses doubles, but is needed for correct distance queries
using Polyhedron = CGAL::Polyhedron_3< CGAL_Kernel >;
using Tree = CGAL::AABB_tree< CGAL::AABB_traits<CGAL_Kernel, CGAL::AABB_face_graph_triangle_primitive<Polyhedron>> >;

Polyhedron MakeGeometry( const InputData & );

Tree MakeAABBTree( const Polyhedron & );


bool PointInside( const Tree &, 
                  const floatType, 
                  const floatType, 
                  const floatType );


floatType NearestRayIntersection( const Tree &,
                                  const fVector3 &,
                                  const fVector3 & );


floatType NearestDistance( const Tree &, 
                           const floatType, 
                           const floatType, 
                           const floatType ); 


}


#endif  // CFD_GEOMETRY