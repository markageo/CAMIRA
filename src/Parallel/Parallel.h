#ifndef CAMIRA_PARALLEL 
#define CAMIRA_PARALLEL

#include "../Core/Types.h"
#include <vector>

namespace CAMIRA
{

std::vector< std::vector<intType> > CreateForward1DColourSet( const intType );

std::vector< std::vector<intType> > CreateReverse1DColourSet( const intType );


} // end namespace CAMIRA

#endif // CAMIRA_PARALLEL