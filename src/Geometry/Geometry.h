#ifndef CFD_GEOMETRY
#define CFD_GEOMETRY

#include "../Types.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../IO/InputProcessing.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

namespace CFD
{

using CGAL_Kernel = CGAL::Simple_cartesian< floatType >;
using Polyhedron = CGAL::Polyhedron_3< CGAL_Kernel >;

// Make CGAL polyhedron from use input block data
Polyhedron MakeGeometry( const InputData & );

}


#endif  // CFD_GEOMETRY