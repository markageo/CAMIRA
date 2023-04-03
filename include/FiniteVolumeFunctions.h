#ifndef FV_FUNCTIONS
#define FV_FUNCTIONS

#include "FiniteVolumeStructures.h"
#include "Types.h"
#include <vector>

namespace CFD
{

// Update face velocities
void UpdateFaceVelocities( ArrayAllocator<Fields>  &, const Mesh &, const ArrayAllocator<Fields> &, 
    const std::vector< std::vector< InputData::BoundaryConditionStruct > > &);

// Initialise finite volume coefficients
void InitialiseFVCoefficients(FVCoefficients &, const Mesh &, const ArrayAllocator<Fields> &, const InputData &);

// Update finite volume coefficients (Picard linearisation)
void UpdateFVCoefficients(FVCoefficients &, const Mesh &, const ArrayAllocator<Fields> &);

} // end namespace CFD

#endif // FV_FUNCTIONS