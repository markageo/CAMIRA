#ifndef CAMIRA_PLUME_PARTICLES
#define CAMIRA_PLUME_PARTICLES

#include "Core/Types.h"
#include "Core/Mesh/Mesh.h"

#include <vector>
#include <algorithm>

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{

// Struct to hold particle data
struct Particles {

    intType N;
    std::vector<floatType> x, y, z;
    std::vector<floatType> mass;
    std::vector<intType> i, j, k;   //index of nearest vertex point in the lo direction
    std::vector<bool> active;

    Particles() : N(0) {};

    void Reserve( intType n )
    {
        x.reserve(n);
        y.reserve(n);
        z.reserve(n);

        i.reserve(n);
        j.reserve(n);
        k.reserve(n);

        mass.reserve(n);

        active.reserve(n);
    }


    void AddParticleBack()
    {
        x.push_back(0.0f);
        y.push_back(0.0f);
        z.push_back(0.0f);

        i.push_back(0);
        j.push_back(0);
        k.push_back(0);

        mass.push_back(0.0f);

        active.push_back(true);

        N++;
    }


    void AddIdenticalParticleBack( intType idx )
    {
        x.push_back( x[idx] );
        y.push_back( y[idx] );
        z.push_back( z[idx] );

        i.push_back( i[idx] );
        j.push_back( j[idx] );
        k.push_back( k[idx] );

        mass.push_back( mass[idx] );

        active.push_back( active[idx] );

        N++;
    }

    
    void RemoveInactiveParticles()
    {

        // This does no preserve order
        for ( intType p = 0; p < N; /* NULL */ ) {

            if ( !active[p] ) {

                SwapAndPop(x, p);
                SwapAndPop(y, p);
                SwapAndPop(z, p);

                SwapAndPop(i, p);
                SwapAndPop(j, p);
                SwapAndPop(k, p);

                SwapAndPop(mass, p);

                SwapAndPop(active, p);

                N--;

            } else {
                p++;
            }

        }

    }

    private:

        template<typename T>
        void SwapAndPop( std::vector<T> &vec,
                         intType idx )
        {
            std::iter_swap( vec.begin() + idx, vec.end() - 1 );     // Can't use std::swap due to bool specialisation of std::vector
            vec.pop_back();
        }

};


namespace 
{

void UpdateCoordinateIndexLinearSearch( intType &index, 
                                        const floatType x,
                                        const Tensor1D &coords )
{

    while ( true ) {

        bool outOfBounds = index > coords.size()-1
                        || index < 0;

        if ( outOfBounds ) {
            std::cout << "x: " << x
                    << ", index: " << index
                    << ", coords.size(): " << coords.size()
                    << std::endl;
        }

        bool found = ( x >= coords( index )     )
                  && ( x <= coords( index + 1 ) );

        if ( found )
            break;

        if ( x < coords( index ) ) {
            index--;
        } else {
            index++;
        }

    }

} 



void UpdateCoordinateIndexBinarySearch( intType &index, 
                                        const floatType x,
                                        const Tensor1D &coords )
{
    intType loBound = 0,
            hiBound = coords.size() - 1;

    while ( loBound <= hiBound ) {

        index = loBound + static_cast<intType>( std::floor( static_cast<floatType>( hiBound - loBound ) / 2.0f ) );

        if        ( coords(index) > x ) {
            hiBound = index;
        } else if ( coords(index) < x ) {
            loBound = index;
        } 
        
        if ( hiBound == loBound + 1 ) {
            break;
        }

    }
}

}   // end namespace anonymous



inline void UpdateParticlePositionIndexLinearSearch( Particles &particles,
                                                     const intType idx,
                                                     const Mesh &mesh )
{
    UpdateCoordinateIndexLinearSearch( particles.i[idx], particles.x[idx], mesh.cellFaces[0] );
    UpdateCoordinateIndexLinearSearch( particles.j[idx], particles.y[idx], mesh.cellFaces[1] );
    UpdateCoordinateIndexLinearSearch( particles.k[idx], particles.z[idx], mesh.cellFaces[2] );    
}



inline void UpdateParticlePositionIndexBinarySearch( Particles &particles,
                                                     const intType idx,
                                                     const Mesh &mesh )
{
    UpdateCoordinateIndexBinarySearch( particles.i[idx], particles.x[idx], mesh.cellFaces[0] );
    UpdateCoordinateIndexBinarySearch( particles.j[idx], particles.y[idx], mesh.cellFaces[1] );
    UpdateCoordinateIndexBinarySearch( particles.k[idx], particles.z[idx], mesh.cellFaces[2] );   
}



inline floatType GetFieldQuantityTrilinearInterp( const Particles &particles,
                                                  const intType idx, 
                                                  const Mesh &mesh,
                                                  const Tensor3D &field )
{
    using enum Axis::ENUMDATA;

    // Lattice coordinates
    const intType &i = particles.i[idx],
                  &j = particles.j[idx],
                  &k = particles.k[idx];

    
    const floatType xd = ( particles.x[idx] - mesh.cellFaces[X]( i ) )
                       / ( mesh.cellFaces[X]( i + 1 ) - mesh.cellFaces[X]( i ) ),

                    yd = ( particles.y[idx] - mesh.cellFaces[Y]( j ) )
                       / ( mesh.cellFaces[Y]( j + 1 ) - mesh.cellFaces[Y]( j ) ),

                    zd = ( particles.z[idx] - mesh.cellFaces[Z]( k ) )
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



inline fVector3 GetFieldQuantityGradient( const Particles &particles, 
                                          const intType idx,
                                          const Mesh &mesh,
                                          const Tensor3D &field )
{
    using enum Axis::ENUMDATA;

    fVector3 grad = {0, 0, 0};

    // Lattice coordinates
    const intType &i = particles.i[idx],
                  &j = particles.j[idx],
                  &k = particles.k[idx];

    
    const floatType xd = ( particles.x[idx] - mesh.cellFaces[X]( i ) )
                       / ( mesh.cellFaces[X]( i + 1 ) - mesh.cellFaces[X]( i ) ),

                    yd = ( particles.y[idx] - mesh.cellFaces[Y]( j ) )
                       / ( mesh.cellFaces[Y]( j + 1 ) - mesh.cellFaces[Y]( j ) ),

                    zd = ( particles.z[idx] - mesh.cellFaces[Z]( k ) )
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

#endif  // CAMIRA_PLUME_PARTICLES