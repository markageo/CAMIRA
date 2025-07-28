#ifndef CAMIRA_PARALLEL 
#define CAMIRA_PARALLEL

#include "../Core/Types.h"
#include <vector>

namespace CAMIRA
{

std::vector< std::vector<intType> > CreateForward1DColorSet( const intType );

std::vector< std::vector<intType> > CreateReverse1DColorSet( const intType );


std::vector< std::vector<intType> > CreateForward2DColorSet( const iArray3 & );

std::vector< std::vector<intType> > CreateReverse2DColorSet( const iArray3 & );


std::vector< std::vector<intType> > CreateForward3DColorSet( const iArray3 & );

std::vector< std::vector<intType> > CreateReverse3DColorSet( const iArray3 & );



} // end namespace CAMIRA

#endif // CAMIRA_PARALLEL