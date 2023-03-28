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

    // In the finite volume formulation, momentum equations are divided by the cell face area 
    // normal to the direction of the momentum equation. This means the pressure coefficients
    // can be stored in 1D arrays, and cell face areas only appear in the diffusion equations.

    FVCoefficients(const CFD::indexVector3 &);
    ArrayAllocator<TransportCoefficients::ENUMDATA, array3D> auu, avv, aww;          // Momentum velocity coefficients
    ArrayAllocator<TransportCoefficients::ENUMDATA, array1D> aup, avp, awp;          // Momentum pressure coefficients
    ArrayAllocator<TransportCoefficients::ENUMDATA, array1D> acu, acv, acw;          // Continuity velocity coefficients
    ArrayAllocator<TransportCoefficients::ENUMDATA, array3D> acp;                    // Continuity pressure coefficients
    array3D                                                  bu, bv, bw, bc;         // Momentum and continuity source terms (appear on the right hand size)     
};


// Update face velocities
void UpdateFaceVelocities( ArrayAllocator<Fields::ENUMDATA>  &, const Mesh &, const ArrayAllocator<Fields::ENUMDATA> &, 
    const std::vector< std::vector< InputData::BoundaryConditionStruct > > &);

// Initialise finite volume coefficients
void InitialiseFVCoefficients(FVCoefficients &, const Mesh &, const ArrayAllocator<Fields::ENUMDATA> &);

// Update finite volume coefficients (Picard linearisation)
void UpdateFVCoefficients(FVCoefficients &, const Mesh &, const ArrayAllocator<Fields::ENUMDATA> &);



} // end namespace CFD

#endif // FV_FUNCTIONS