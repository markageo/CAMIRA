#ifndef FV_FUNCTIONS
#define FV_FUNCTIONS

#include "FiniteVolumeStructures.h"
#include "Types.h"
#include <vector>

namespace CFD
{

// Update face velocities
void UpdateFaceVelocities( ArrayAllocator<Fields, CFD::array3D>  &, const Mesh &, const ArrayAllocator<Fields, CFD::array3D> &, 
    const std::vector< std::vector< InputData::BoundaryConditionStruct > > &);

// Initialise finite volume coefficients
void InitialiseFVCoefficients(FVCoefficients &, const Mesh &, const ArrayAllocator<Fields, CFD::array3D> &, const InputData &);

// Update finite volume coefficients (Picard linearisation)
void UpdateFVCoefficients(FVCoefficients &, const Mesh &, const ArrayAllocator<Fields, CFD::array3D> &);

} // end namespace CFD

#endif // FV_FUNCTIONS