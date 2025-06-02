#ifndef CFD_PARALLEL 
#define CFD_PARALLEL

#include "../Core/Types.h"
#include <vector>

namespace CFD
{

std::vector< std::vector<intType> > CreateForward1DColourSet( const intType );

std::vector< std::vector<intType> > CreateReverse1DColourSet( const intType );


} // end namespace CFD

#endif // CFD_PARALLEL