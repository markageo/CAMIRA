#ifndef CFD_PARALLEL 
#define CFD_PARALLEL

#include "../Core/Types.h"
#include "RAJA/RAJA.hpp"
#include "camp/resource.hpp"


namespace CFD
{

RAJA::TypedIndexSet<RAJA::ListSegment> CreateForward3ColorSet( const iArray3 &,
                                                               camp::resources::Resource );

RAJA::TypedIndexSet<RAJA::ListSegment> CreateReverse3ColorSet( const iArray3 &,
                                                               camp::resources::Resource );


RAJA::TypedIndexSet<RAJA::ListSegment> Create3ColorSetColumns( const iArray3 &,
                                                               camp::resources::Resource );                                         


RAJA::TypedIndexSet<RAJA::ListSegment> Create3ColorSetPlane( const iArray3 &,
                                                             camp::resources::Resource );                                         

} // end namespace CFD

#endif // CFD_PARALLEL