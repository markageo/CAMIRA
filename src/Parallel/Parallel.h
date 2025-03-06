#ifndef CFD_PARALLEL 
#define CFD_PARALLEL

#include "../Core/Types.h"
#include "RAJA/RAJA.hpp"
#include "camp/resource.hpp"


namespace CFD
{

RAJA::TypedIndexSet<RAJA::ListSegment> Create3ColorSet( const iArray3 &,
                                                        camp::resources::Resource );

} // end namespace CFD

#endif // CFD_PARALLEL