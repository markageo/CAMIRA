#include "Concentration.h"


namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{

namespace
{
}   // end anonymous namespace


void UpdateConcentrationField( Tensor3D &concentrationField,
                               const std::vector<Particle> &particles,
                               const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    // We will add concentration to each cell cumulatively
    concentrationField.setZero();

    for ( const auto &particle : particles ) {

        const intType i = particle.positionIndex[X],
                      j = particle.positionIndex[Y],
                      k = particle.positionIndex[Z];

        // Weighting factors by distance to cell center
        const floatType wx = ( mesh.cellCenters[X]( i+1 ) - particle.position(X) ) * mesh.cellCenterDiffInv[X]( i+1 ),  
                        wy = ( mesh.cellCenters[Y]( j+1 ) - particle.position(Y) ) * mesh.cellCenterDiffInv[Y]( j+1 ),  
                        wz = ( mesh.cellCenters[Z]( k+1 ) - particle.position(Z) ) * mesh.cellCenterDiffInv[Z]( k+1 );  

        // Add concentrations
        concentrationField(i  , j  , k  ) += wx * wy * wz * particle.mass 
                                           / ( mesh.cellLengths[X](i) * mesh.cellLengths[Y](j) * mesh.cellLengths[Z](k) );

        concentrationField(i+1, j  , k  ) += ( 1.0f - wx ) * wy * wz * particle.mass
                                           / ( mesh.cellLengths[X](i+1) * mesh.cellLengths[Y](j) * mesh.cellLengths[Z](k) );


        concentrationField(i  , j+1, k  ) += wx * ( 1.0f - wy ) * wz * particle.mass
                                           / ( mesh.cellLengths[X](i) * mesh.cellLengths[Y](j+1) * mesh.cellLengths[Z](k) );

        concentrationField(i  , j  , k+1) += wx * wy * ( 1.0f - wz ) * particle.mass
                                           / ( mesh.cellLengths[X](i) * mesh.cellLengths[Y](j) * mesh.cellLengths[Z](k+1) );


        concentrationField(i+1, j  , k+1) += ( 1.0f - wx ) * wy * ( 1.0f - wz ) * particle.mass
                                           / ( mesh.cellLengths[X](i+1) * mesh.cellLengths[Y](j) * mesh.cellLengths[Z](k+1) );


        concentrationField(i+1, j+1, k  ) += ( 1.0f - wx ) * ( 1.0f - wy ) * wz * particle.mass
                                           / ( mesh.cellLengths[X](i+1) * mesh.cellLengths[Y](j+1) * mesh.cellLengths[Z](k) );


        concentrationField(i  , j+1, k+1) += wx * ( 1.0f - wy ) * ( 1.0f - wz ) * particle.mass
                                           / ( mesh.cellLengths[X](i) * mesh.cellLengths[Y](j+1) * mesh.cellLengths[Z](k+1) );


        concentrationField(i+1, j+1, k+1) += ( 1.0f - wx ) * ( 1.0f - wy ) * ( 1.0f - wz ) * particle.mass
                                           / ( mesh.cellLengths[X](i+1) * mesh.cellLengths[Y](j+1) * mesh.cellLengths[Z](k+1) );

    }

}



}   // end namespace PLUME

}   // end namespace CAMIRA
