#ifndef CAMIRA_PLUME_SOURCES
#define CAMIRA_PLUME_SOURCES

#include "Core/Types.h"
#include "Plume/InputProcessing/InputProcessing.h"
#include "Plume/Particle/Particle.h"

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{

void AddInstantaneousReleasePointParticles( std::vector<Particle> &,
                                            const InputData & );



void AddContinuousReleasePointParticles( std::vector<Particle> &,
                                         const InputData & );


}   // end namespace PLUME

}   // end namespace CAMIRA

#endif  // CAMIRA_PLUME_SOURCES