#include "Concentration.h"
#include <algorithm>
#include <omp.h>

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{


void UpdateConcentrationField( Tensor3D &concentrationField,
                               const std::vector<Particle> &particles,
                               const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    // We will add concentration to each cell cumulatively
    SetTensorZeroParallel( concentrationField );

    #pragma omp parallel for
    for ( const auto &particle : particles ) {

        // Weighting is by cell center, we have cell face index for the face immediately to the lo side of the particle
        // Need to find index of cell center immediately to the lo side
        const intType i = ( particle.position(X) > mesh.cellCenters[X]( particle.positionIndex[X] ) ) ? particle.positionIndex[X] : std::max( static_cast<intType>(0), particle.positionIndex[X] - 1 ),
                      j = ( particle.position(Y) > mesh.cellCenters[Y]( particle.positionIndex[Y] ) ) ? particle.positionIndex[Y] : std::max( static_cast<intType>(0), particle.positionIndex[Y] - 1 ),
                      k = ( particle.position(Z) > mesh.cellCenters[Z]( particle.positionIndex[Z] ) ) ? particle.positionIndex[Z] : std::max( static_cast<intType>(0), particle.positionIndex[Z] - 1 );

        // These are values that are used if at domain boundary, will get overwritten if not.
        intType ip1 = i,
                jp1 = j,
                kp1 = k;

        floatType wx = 0.5f,
                  wy = 0.5f,  
                  wz = 0.5f;  

        // Weighting factors by distance to cell center
        if ( i > 0 && i < mesh.nCells[X]-1 ) {
            ip1 = i+1;
            wx = ( mesh.cellCenters[X]( ip1 ) - particle.position(X) ) * mesh.cellCenterDiffInv[X]( ip1 );
        } 

        if ( j > 0 && j < mesh.nCells[Y]-1 ) {
            jp1 = j+1;
            wy = ( mesh.cellCenters[Y]( jp1 ) - particle.position(Y) ) * mesh.cellCenterDiffInv[Y]( jp1 );
        }

        if ( k > 0 && k < mesh.nCells[Z]-1 ) {
            kp1 = k+1;
            wz = ( mesh.cellCenters[Z]( kp1 ) - particle.position(Z) ) * mesh.cellCenterDiffInv[Z]( kp1 );
        }

        // Calculate concentrations
        const floatType c000 = wx * wy * wz * particle.mass 
                             / ( mesh.cellLengths[X](i) * mesh.cellLengths[Y](j) * mesh.cellLengths[Z](k) ),
                        
                        c100 = ( 1.0f - wx ) * wy * wz * particle.mass
                             / ( mesh.cellLengths[X](ip1) * mesh.cellLengths[Y](j) * mesh.cellLengths[Z](k) ),

                        c010 = wx * ( 1.0f - wy ) * wz * particle.mass
                             / ( mesh.cellLengths[X](i) * mesh.cellLengths[Y](jp1) * mesh.cellLengths[Z](k) ),

                        c001 = wx * wy * ( 1.0f - wz ) * particle.mass
                             / ( mesh.cellLengths[X](i) * mesh.cellLengths[Y](j) * mesh.cellLengths[Z](kp1) ),

                        c101 = ( 1.0f - wx ) * wy * ( 1.0f - wz ) * particle.mass
                             / ( mesh.cellLengths[X](ip1) * mesh.cellLengths[Y](j) * mesh.cellLengths[Z](kp1) ),

                        c110 = ( 1.0f - wx ) * ( 1.0f - wy ) * wz * particle.mass
                             / ( mesh.cellLengths[X](ip1) * mesh.cellLengths[Y](jp1) * mesh.cellLengths[Z](k) ),

                        c011 = wx * ( 1.0f - wy ) * ( 1.0f - wz ) * particle.mass
                             / ( mesh.cellLengths[X](i) * mesh.cellLengths[Y](jp1) * mesh.cellLengths[Z](kp1) ),

                        c111 = ( 1.0f - wx ) * ( 1.0f - wy ) * ( 1.0f - wz ) * particle.mass
                             / ( mesh.cellLengths[X](ip1) * mesh.cellLengths[Y](jp1) * mesh.cellLengths[Z](kp1) );
        

        // Add concentrations
        #pragma omp atomic
        concentrationField(i  , j  , k  ) += c000;

        #pragma omp atomic
        concentrationField(ip1, j  , k  ) += c100;

        #pragma omp atomic
        concentrationField(i  , jp1, k  ) += c010;

        #pragma omp atomic
        concentrationField(i  , j  , kp1) += c001;

        #pragma omp atomic
        concentrationField(ip1, j  , kp1) += c101;

        #pragma omp atomic
        concentrationField(ip1, jp1, k  ) += c110;

        #pragma omp atomic
        concentrationField(i  , jp1, kp1) += c011;

        #pragma omp atomic
        concentrationField(ip1, jp1, kp1) += c111;

    }

}



}   // end namespace PLUME

}   // end namespace CAMIRA
