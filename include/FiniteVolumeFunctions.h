#ifndef FV_FUNCTIONS
#define FV_FUNCTIONS

#include "FiniteVolumeStructures.h"
#include "Types.h"

namespace CFD
{

// Structure to store finite volume discrete equation coefficients (Picard linearisation)
struct FVCoefficients
{
    // Naming convention:
    //  aev - coefficient for equation 'e' multiplying with variable 'v'
    //  be  - source term for equation 'e'
    // 
    // 'e' can take the values
    //  u: U momentum, v: V momentum, w: W momentum, c: Conitnuity
    //
    // 'v' can take the values
    //  u: x velocity, v: v velocity, w: w velocity, p: pressure

    FVCoefficients(const CFD::indexVector3 &);
    CFD::ArrayAllocator<CFD::TransportCoefficients::ENUMDATA, CFD::array3D> auu, avv, aww;          // Momentum velocity coefficients
    CFD::ArrayAllocator<CFD::TransportCoefficients::ENUMDATA, CFD::array1D> aup, avp, awp;          // Momentum pressure coefficients
    CFD::ArrayAllocator<CFD::TransportCoefficients::ENUMDATA, CFD::array1D> acu, acv, acw; //acp;     // Continuity velocity and pressure coefficients
    CFD::array3D                                                            bu, bv, bw, bc;         // Source terms for momentum and continuity equations 
};


// Update face velocities
void UpdateFaceVelocities( ArrayAllocator<CFD::Fields::ENUMDATA>  &, const Mesh &, const ArrayAllocator<CFD::Fields::ENUMDATA> &, 
    const std::vector< std::vector< InputData::BoundaryConditionStruct > > &);

} // end namespace CFD

#endif // FV_FUNCTIONS