#include "Concentration.h"
#include <algorithm>
#include <omp.h>

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{


void UpdateConcentrationField( Tensor3D &concentrationField,
                               const Particles &particles,
                               const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    // We will add concentration to each cell cumulatively
    SetTensorZeroParallel( concentrationField );

    for ( intType idx = 0; idx != particles.N; idx++ ) {

        const floatType &mass = particles.mass[idx];

        const floatType &x = particles.x[idx],
                        &y = particles.y[idx],
                        &z = particles.z[idx];

        // Weighting is by cell center, we have cell face index for the face immediately to the lo side of the particle
        // Need to find index of cell center immediately to the lo side
        intType i = particles.i[idx],
                j = particles.j[idx],
                k = particles.k[idx];

        i = ( x > mesh.cellCenters[X]( i ) ) ? i : std::max( static_cast<intType>(0), i - 1 ),
        j = ( y > mesh.cellCenters[Y]( j ) ) ? j : std::max( static_cast<intType>(0), j - 1 ),
        k = ( z > mesh.cellCenters[Z]( k ) ) ? k : std::max( static_cast<intType>(0), k - 1 );

          const floatType concentration = mass / ( mesh.cellLengths[X](i) * mesh.cellLengths[Y](j) * mesh.cellLengths[Z](k) );

          // Add concentrations
          #pragma omp atomic
          concentrationField(i, j, k) += concentration;

    }
    
}



}   // end namespace PLUME

}   // end namespace CAMIRA
