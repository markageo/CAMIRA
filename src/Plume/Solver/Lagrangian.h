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


[[maybe_unused]]
void UpdateParticlePositionIndexLinearSearch( Particle &particle,
                                              const Mesh &mesh )
{

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        
        while ( true ) {

            bool found = particle.position(axis) >= mesh.cellFaces[axis]( particle.positionIndex[axis] )
                      && particle.position(axis) <= mesh.cellFaces[axis]( particle.positionIndex[axis] + 1 );

            if ( found )
                break;

            if ( particle.position(axis) < mesh.cellFaces[axis]( particle.positionIndex[axis] ) ) {
                particle.positionIndex[axis]--;
            } else {
                particle.positionIndex[axis]++;
            }

        }

    } );
    
}


[[maybe_unused]]
void UpdateParticlePositionIndexBinarySearch( Particle &particle,
                                              const Mesh &mesh )
{
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        
        intType loBound = 0,
                hiBound = mesh.nFacesNormal[axis][axis] - 1;

        intType index;

        while ( loBound <= hiBound ) {

            index = loBound + static_cast<intType>( std::floor( static_cast<floatType>( hiBound - loBound ) / 2.0f ) );

            if        ( mesh.cellFaces[axis](index) > particle.position(axis) ) {
                hiBound = index;
            } else if ( mesh.cellFaces[axis](index) < particle.position(axis) ) {
                loBound = index;
            } else if ( hiBound == loBound + 1 ) {
                break;
            }

        }

        particle.positionIndex[axis] = index;

    } );
}



void UpdateParticleVelocityTrilinearInterp( Particle &particle, 
                                            const Mesh &mesh,
                                            const EnumVector<Axis, Tensor3D> &velocityField )
{
    using enum Axis::ENUMDATA;

    // Lattice coordinates
    const intType i = particle.positionIndex[X],
                  j = particle.positionIndex[Y],
                  k = particle.positionIndex[Z];

    
    const floatType xd = ( particle.position(X) - mesh.cellFaces[X]( i ) )
                       / ( mesh.cellFaces[X]( i + 1 ) - mesh.cellFaces[X]( i ) ),

                    yd = ( particle.position(Y) - mesh.cellFaces[Y]( j ) )
                       / ( mesh.cellFaces[Y]( j + 1 ) - mesh.cellFaces[Y]( j ) ),

                    zd = ( particle.position(Z) - mesh.cellFaces[Z]( k ) )
                       / ( mesh.cellFaces[Z]( k + 1 ) - mesh.cellFaces[Z]( k ) );

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        const Tensor3D &field = velocityField[axis];

        // Velocity at lattice points
        const floatType c000 = field( i  , j  , k   ),
                        c100 = field( i+1, j  , k   ),
                        c010 = field( i  , j+1, k   ),
                        c001 = field( i  , j  , k+1 ),
                        c101 = field( i+1, j  , k+1 ),
                        c011 = field( i  , j+1, k+1 ),
                        c110 = field( i+1, j+1, k   ),
                        c111 = field( i+1, j+1, k+1 );

        // Interpolate x direction
        const floatType c00 = c000 * ( 1 - xd )  +  c100 * xd,
                        c01 = c001 * ( 1 - xd )  +  c101 * xd,
                        c10 = c010 * ( 1 - xd )  +  c110 * xd,
                        c11 = c011 * ( 1 - xd )  +  c111 * xd;

        // Interpolate y direction
        const floatType c0 = c00 * ( 1 - yd )  +  c10 * yd,
                        c1 = c01 * ( 1 - yd )  +  c11 * yd;

        // Interpolate z direction
        const floatType c = c0 * ( 1 - zd )  +  c1 * zd;

        particle.velocity(axis) = c;

    } );
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

        floatType hiBounds = mesh.cellFaces[axis]( mesh.nFacesNormal[axis][axis] ),
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



void UpdateParticles( std::vector<Particle> &particles,
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
                      ) 
                    );

}


}   // end namespace PLUME

}   // end namespace CAMIRA    


#endif // CAMIRA_PLUME_LAGRANGIAN