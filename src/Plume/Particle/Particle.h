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
    fVector3 position = {0, 0, 0};
    floatType mass = 0;
    TensorIndex3D positionIndex = {0, 0, 0};    // Of nearest vertex point in the lo direction
    fVector3 velocity = {0, 0, 0};
    bool active = true;
};


}   // end namespace PLUME

}   // end namespace CAMIRA

#endif  // CAMIRA_PLUME_PARTICLE