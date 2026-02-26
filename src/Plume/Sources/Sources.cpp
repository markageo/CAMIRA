#include "Plume/Sources/Sources.h"

namespace CAMIRA
{
    
using namespace CORE;

namespace PLUME
{


namespace 
{


void AddParticle( std::vector<Particle> &particles,
                  const InputData &inputData,
                  const intType nParticles,
                  const floatType initialMassRemaining,
                  const fArray3 &releaseLocation )
{
    floatType massRemaining = initialMassRemaining;
    floatType massPerParticle = 1.0f / inputData.initialParticlesPerUnitMass;
    for ( intType i = 0; i != nParticles; i++ ) {

        floatType mass = ( massRemaining >= massPerParticle ) ? massPerParticle
                                                              : massRemaining; 

        particles.emplace_back( releaseLocation, mass );

        massRemaining -= inputData.initialParticlesPerUnitMass;
    }

}


}   // end anonymous namespace



void AddInstantaneousReleasePointParticles( std::vector<Particle> &particles,
                                            const InputData &inputData )
{

    for ( const auto &release : inputData.instantaneousReleasePoints ) {

        // Total number of particles release all at once
        intType nParticles = static_cast<intType>( std::ceil( release.totalMass * inputData.initialParticlesPerUnitMass ) );
        floatType massRemaining = release.totalMass;

        AddParticle( particles, inputData, nParticles, massRemaining, release.location );
    }

}



void AddContinuousReleasePointParticles( std::vector<Particle> &particles,
                                         const InputData &inputData)
{

    for ( const auto &release : inputData.continuousReleasePoints ) {

        // Total number of particles that will be released in this timestep
        intType nParticles = static_cast<intType>( std::ceil( release.massFlowRate * inputData.timeStepSize * inputData.initialParticlesPerUnitMass ) );
        floatType massRemaining = release.massFlowRate * inputData.timeStepSize;

        AddParticle( particles, inputData, nParticles, massRemaining, release.location );
    }

}


}   // end namespace PLUME

}   // end namespace CAMIRA