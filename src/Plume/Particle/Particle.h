#ifndef CAMIRA_PLUME_PARTICLE
#define CAMIRA_PLUME_PARTICLE

#include "Core/Types.h"

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{

// Struct to hold particle data
struct Particle {
    fArray3 position;
    floatType mass;
};


}   // end namespace PLUME

}   // end namespace CAMIRA

#endif  // CAMIRA_PLUME_PARTICLE