#ifndef CFD_GEOMETRY
#define CFD_GEOMETRY

#include "../IO/InputProcessing.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


namespace CFD
{

using CGAL_Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;    // This uses doubles, but is needed for correct distance queries
using Polyhedron = CGAL::Polyhedron_3< CGAL_Kernel >;

// Make CGAL polyhedron from use input block data
Polyhedron MakeGeometry( const InputData & );

}


#endif  // CFD_GEOMETRY