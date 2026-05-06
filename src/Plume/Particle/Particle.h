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
    TensorIndex3D positionIndex = {0, 0, 0};    // Of nearest vertex point in the lo direction
    floatType mass = 0;
    bool active = true;
};




inline void UpdateParticlePositionIndexLinearSearch( Particle &particle,
                                                     const Mesh &mesh )
{

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        
        while ( true ) {

            bool found = ( particle.position(axis) >= mesh.cellFaces[axis]( particle.positionIndex[axis] )     )
                      && ( particle.position(axis) <= mesh.cellFaces[axis]( particle.positionIndex[axis] + 1 ) );

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
            } 
            
            if ( hiBound == loBound + 1 ) {
                break;
            }

        }

    } );
}



inline floatType GetFieldQuantityTrilinearInterp( const Particle &particle, 
                                                  const Mesh &mesh,
                                                  const Tensor3D &field )
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


    // Quantity at lattice points
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

    return c;
}



inline fVector3 GetFieldQuantityGradient( const Particle &particle, 
                                          const Mesh &mesh,
                                          const Tensor3D &field )
{
    using enum Axis::ENUMDATA;

    fVector3 grad = {0, 0, 0};

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


    const floatType dx = mesh.cellLengths[X](i),
                    dy = mesh.cellLengths[Y](j),
                    dz = mesh.cellLengths[Z](k);


    // Quantity at lattice points
    const floatType c000 = field( i  , j  , k   ),
                    c100 = field( i+1, j  , k   ),
                    c010 = field( i  , j+1, k   ),
                    c001 = field( i  , j  , k+1 ),
                    c101 = field( i+1, j  , k+1 ),
                    c011 = field( i  , j+1, k+1 ),
                    c110 = field( i+1, j+1, k   ),
                    c111 = field( i+1, j+1, k+1 );


                    
    // x gradient
    const floatType gx00 = ( c100 - c000 ) / dx,
                    gx10 = ( c110 - c010 ) / dx,
                    gx01 = ( c101 - c001 ) / dx,
                    gx11 = ( c111 - c011 ) / dx;

    // Interpolate y direction
    const floatType gx0 = gx00 * ( 1 - yd )  +  gx10 * yd,
                    gx1 = gx01 * ( 1 - yd )  +  gx11 * yd;

    // Interpolate z direction
    const floatType gx = gx0 * ( 1 - zd )  +  gx1 * zd;
    grad(X) = gx;



    // y gradient
    const floatType gy00 = ( c010 - c000 ) / dy,
                    gy10 = ( c110 - c100 ) / dy,
                    gy01 = ( c011 - c001 ) / dy,
                    gy11 = ( c111 - c101 ) / dy;

    // Interpolate x direction
    const floatType gy0 = gy00 * ( 1 - xd )  +  gy10 * xd,
                    gy1 = gy01 * ( 1 - xd )  +  gy11 * xd;

    // Interpolate z direction
    const floatType gy = gy0 * ( 1 - zd )  +  gy1 * zd;
    grad(Y) = gy;



    // z gradient
    const floatType gz00 = ( c001 - c000 ) / dz,
                    gz10 = ( c101 - c100 ) / dz,
                    gz01 = ( c011 - c010 ) / dz,
                    gz11 = ( c111 - c110 ) / dz;

    // Interpolate x direction
    const floatType gz0 = gz00 * ( 1 - xd )  +  gz10 * xd,
                    gz1 = gz01 * ( 1 - xd )  +  gz11 * xd;

    // Interpolate y direction
    const floatType gz = gz0 * ( 1 - yd )  +  gz1 * yd;
    grad(Z) = gz;

    return grad;
}


}   // end namespace PLUME

}   // end namespace CAMIRA

#endif  // CAMIRA_PLUME_PARTICLE