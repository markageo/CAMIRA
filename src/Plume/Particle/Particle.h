#ifndef CAMIRA_PLUME_PARTICLE
#define CAMIRA_PLUME_PARTICLE

#include "Core/Types.h"
#include "Core/Mesh/Mesh.h"

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



inline void UpdateParticlePositionIndexLinearSearch( Particle &particle,
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


inline void UpdateParticlePositionIndexBinarySearch( Particle &particle,
                                                     const Mesh &mesh )
{
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        
        intType loBound = 0,
                hiBound = mesh.nFacesNormal[axis][axis] - 1;

        intType &index = particle.positionIndex[axis];

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

    } );
}



inline void UpdateParticleVelocityTrilinearInterp( Particle &particle, 
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


}   // end namespace PLUME

}   // end namespace CAMIRA

#endif  // CAMIRA_PLUME_PARTICLE