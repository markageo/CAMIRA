#ifndef CAMIRA_PLUME_SOURCES
#define CAMIRA_PLUME_SOURCES

#include "Core/Types.h"
#include "Core/Mesh/Mesh.h"
#include "Plume/InputProcessing/InputProcessing.h"
#include "Plume/Particles/Particles.h"

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{

void AddInstantaneousReleasePointParticles( Particles &,
                                            const Mesh &,
                                            const InputData & );



void AddContinuousReleasePointParticles( Particles &,
                                         const Mesh &,
                                         const InputData & );


}   // end namespace PLUME

}   // end namespace CAMIRA

#endif  // CAMIRA_PLUME_SOURCES