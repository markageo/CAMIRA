#ifndef SOLVER
#define SOLVER

#include "InputProcessing.h"
#include "Types.h"
#include "FiniteVolumeStructures.h"

namespace CFD
{

    // // Structure to store finite volume discrete equation coefficients (Picard linearisation)
    // struct FVCoefficients
    // {
    //     // Naming convention:
    //     //  aev - coefficient for equation 'e' multiplying with variable 'v'
    //     //  be  - source term for equation 'e'
    //     // 
    //     // 'e' can take the values
    //     //  u: U momentum, v: V momentum, w: W momentum, c: Conitnuity
    //     //
    //     // 'v' can take the values
    //     //  u: x velocity, v: v velocity, w: w velocity, p: pressure

    //     FVCoefficients(const CFD::indexVector3 &);
    //     CFD::ArrayAllocator<CFD::TransportCoefficients::ENUMDATA> auu, aup, bu, 
    //                                                               avv, avp, bv, 
    //                                                               aww, awp, bw, 
    //                                                               acu, acv, acw, acp, bc;
    // };

    void SweepSolve(CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &, const CFD::Mesh &, const CFD::InputData &);

}   // end namespace CFD    


#endif // SOLVER