#ifndef CAMIRA_PLUME_LAGRANGIAN
#define CAMIRA_PLUME_LAGRANGIAN

#include "Core/Types.h"
#include "Core/FVLookups.h"
#include "Core/Mesh/Mesh.h"
#include "Plume/Particle/Particle.h"
#include "Plume/InputProcessing/InputProcessing.h"
#include "Plume/ConfigEnums.h"

#include <vector>
#include <algorithm>
#include <cmath>

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{


namespace 
{

floatType signum( floatType f) {
    if ( f > 0.0f ) return  1.0f;
    if ( f < 0.0f ) return -1.0f;
    return 0.0f;
}



void UpdateParticleVelocity( Particle &particle,
                             const Mesh &mesh,
                             const EnumVector<Axis, Tensor3D> &velocityField )
{
    UpdateParticlePositionIndexLinearSearch( particle, mesh );
    UpdateParticleVelocityTrilinearInterp( particle, mesh, velocityField );
}



fVector3 GetParticleStep( const fVector3 &velocity,
                          const floatType &timeStep )
{
    return timeStep * velocity;
}



void StepParticle( Particle &particle,
                   const fVector3 &delta,
                   const Mesh &mesh, 
                   const EnumVector<BoundaryPatches, InputData::BoundaryConditionInputData> &boundaryConditions )
{
    using enum Axis::ENUMDATA;

    // New position
    fVector3 newPosition = particle.position + delta;

    // Check where it is out of bounds
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        floatType hiBounds = mesh.cellFaces[axis]( mesh.nFacesNormal[axis][axis]-1 ),
                  loBounds = mesh.cellFaces[axis]( 0 );

        BoundaryPatches::ENUMDATA boundaryPatch;
        bool exitedBounds = false;
        
        if        ( newPosition(axis) >= hiBounds ) {
            boundaryPatch = LUT::PositivePatch[axis];
            exitedBounds = true;
        } else if ( newPosition(axis) <= loBounds ) {
            boundaryPatch = LUT::NegativePatch[axis];
            exitedBounds = true;
        }


        if ( exitedBounds ) {

            switch ( boundaryConditions[ boundaryPatch ].type ) {

                case BoundaryConditions::outflow:
                    particle.active = false;
                    break;

                case BoundaryConditions::reflection:
                    newPosition(axis) = particle.position(axis) - delta(axis);  // Can be improved
                    break;

                case BoundaryConditions::periodic:
                    newPosition(axis) = signum( delta(axis) ) * ( loBounds - hiBounds ) + particle.position(axis) + delta(axis); 
                    break;
            }

        }

    } );

    particle.position = newPosition;

}



}   // end anonymous namespace



inline void UpdateParticles( std::vector<Particle> &particles,
                      const Mesh &mesh, 
                      const EnumVector<Axis, Tensor3D> &velocityField,
                      const InputData &inputData )
{

    for ( Particle &particle : particles ) {

        UpdateParticleVelocity( particle, mesh, velocityField );

        fVector3 delta = GetParticleStep( particle.velocity, inputData.timeStepSize );

        StepParticle( particle, delta, mesh, inputData.boundaryConditions );

    }

    // Remove any inactive particles
    // This might not be the most efficient since it preserves order, and we don't need that
    particles.erase( 
        std::remove_if( particles.begin(), 
                        particles.end(), 
                        [] (const Particle &p) { return !p.active; } 
                      ),
        particles.end()
                    );

}




}   // end namespace PLUME

}   // end namespace CAMIRA    


#endif // CAMIRA_PLUME_LAGRANGIAN