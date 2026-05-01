#ifndef CAMIRA_PLUME_LAGRANGIAN
#define CAMIRA_PLUME_LAGRANGIAN

#include "Core/Types.h"
#include "Core/FVLookups.h"
#include "Core/Mesh/Mesh.h"
#include "Core/Geometry/Geometry.h"
#include "Plume/Particle/Particle.h"
#include "Plume/InputProcessing/InputProcessing.h"
#include "Plume/ConfigEnums.h"

#include <vector>
#include <algorithm>
#include <cmath>
#include <random>

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



fVector3 GetParticleStep( const fVector3 &velocity,
                          const floatType turbulentDiffusivity,
                          const fVector3 &gradTurbulentDiffusivity,
                          const InputData &inputData )
{
    // Need to be sure these static variables are thread safe!
    static std::random_device rd{};
    static std::mt19937 gen{rd()};
    static std::normal_distribution d{0.0f, 1.0f};

    const fVector3 advection = inputData.timeStepSize * velocity;

    const fVector3 eta = { d(gen), d(gen), d(gen) };

    const fVector3 molecularDiffusion = sqrt( 2.0f * inputData.diffusionCoeff * inputData.timeStepSize ) * eta;

    const fVector3 turbulentDiffusion = gradTurbulentDiffusivity * inputData.timeStepSize
                                      + sqrt( 2.0f * turbulentDiffusivity * inputData.timeStepSize ) * eta;

    return advection + molecularDiffusion + turbulentDiffusion;
}



fVector3 RecursiveReflection( const Tree &tree,
                              const fVector3 oldPosition,
                              const fVector3 newPosition )
{
    const fVector3 currentDelta = newPosition - oldPosition;
    const auto intersectionData = SegmentIntersectionAndNormal( tree, oldPosition, currentDelta );
    const fVector3 intersectionPoint = intersectionData.point;
    const fVector3 normalVector = intersectionData.normal; 

    fVector3 reflectedPosition = newPosition - 2.0f * ( normalVector.dot( newPosition - intersectionPoint )  ) * normalVector;

    // Shift the intersection point along delta to prevent an unintended intersection from initial intersection point
    const floatType shiftScale = 1e-8;
    const fVector3 shiftedIntersectionPoint = intersectionPoint + shiftScale * (reflectedPosition - intersectionPoint);

    if ( SegmentIntersects( tree, shiftedIntersectionPoint, reflectedPosition ) ) {
        reflectedPosition = RecursiveReflection( tree, shiftedIntersectionPoint, reflectedPosition );
    }

    return reflectedPosition;
}



void StepParticle( Particle &particle,
                   const fVector3 &delta,
                   const Mesh &mesh, 
                   const Tree &tree,
                   const EnumVector<BoundaryPatches, InputData::BoundaryConditionInputData> &boundaryConditions )
{
    using enum Axis::ENUMDATA;

    // New position
    fVector3 newPosition = particle.position + delta;

    // Domain boundary intersection
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
                {
                    const floatType wallPosition = ( newPosition(axis) >= hiBounds ) ? hiBounds 
                                                                                     : loBounds;
                    newPosition(axis) = 2.0f * wallPosition - newPosition(axis);
                    break;
                }

                case BoundaryConditions::periodic:
                    newPosition(axis) = signum( delta(axis) ) * ( loBounds - hiBounds ) + particle.position(axis) + delta(axis); 
                    break;
            }

        }

    } );


    // Solid geometry intersection and reflection 
    if ( SegmentIntersects( tree, particle.position, newPosition ) ) {
        newPosition = RecursiveReflection( tree, particle.position, newPosition );
    }

    particle.position = newPosition;

}



}   // end anonymous namespace



inline void UpdateParticles( std::vector<Particle> &particles,
                             const Mesh &mesh, 
                             const EnumVector<Axis, Tensor3D> &velocityField,
                             const Tensor3D &nuTurbField,
                             const Tree &tree,
                             const InputData &inputData )
{

    for ( Particle &particle : particles ) {

        UpdateParticlePositionIndexLinearSearch( particle, mesh );

        fVector3 localVelocity;
        EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {
            localVelocity[axis] = GetFieldQuantityTrilinearInterp( particle, mesh, velocityField[axis] );
        } );

        const floatType turbulentDiffusivity = GetFieldQuantityTrilinearInterp( particle, mesh, nuTurbField );

        const fVector3 gradTurbulentDiffisivity = GetFieldQuantityGradient( particle, mesh, nuTurbField );

        const fVector3 delta = GetParticleStep( localVelocity, turbulentDiffusivity, gradTurbulentDiffisivity, inputData );

        StepParticle( particle, delta, mesh, tree, inputData.boundaryConditions );

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


inline void SplitParticles( std::vector<Particle> &particles )
{

    const intType initialNumberOfParticles = particles.size();

    for ( intType i = 0; i != initialNumberOfParticles; i++ ) {

        // Half the mass of the particle
        particles[i].mass /= 2.0f;

        // Copy the particle to the end
        particles.push_back( particles[i] );

    }

}




}   // end namespace PLUME

}   // end namespace CAMIRA    


#endif // CAMIRA_PLUME_LAGRANGIAN