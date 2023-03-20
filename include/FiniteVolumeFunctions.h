#ifndef FV_FUNCTIONS
#define FV_FUNCTIONS

#include "FiniteVolumeStructures.h"
#include "Types.h"
#include "Tensor"

namespace CFD
{

// Update face velocities
void UpdateFaceVelocities( ArrayAllocator<CFD::Fields::ENUMDATA>  &, const Mesh &, const ArrayAllocator<CFD::Fields::ENUMDATA> &, 
    const std::vector< std::vector< InputData::BoundaryConditionStruct > > &);

} // end namespace CFD

#endif // FV_FUNCTIONS