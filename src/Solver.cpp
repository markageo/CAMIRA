#include "Solver.h"
#include "InputProcessing.h"
#include "Types.h"
#include "FiniteVolumeStructures.h"


using C = CFD::TransportCoefficients::ENUMDATA;

// CFD::FVCoefficients::FVCoefficients(const CFD::indexVector3 &dims) :
//     auu({C::p, C::n, C::e, C::s, C::w, C::t, C::b}, dims),
//     aup({C::p, C::e, C::w}, dims),

//     avv({C::p, C::n, C::e, C::s, C::w, C::t, C::b}, dims),
//     avp({C::p, C::n, C::s}, dims),

//     aww({C::p, C::n, C::e, C::s, C::w, C::t, C::b}, dims),
//     awp({C::p, C::t, C::b}, dims),
//     {}



void CFD::SweepSolve(ArrayAllocator<CFD::Fields::ENUMDATA>  &fields, const Mesh &mesh, const InputData &inputData) 
{

    // Initialise cell face velocities


    // Initialise transport equation coefficients

    
    // Sweep

        // Update cell face velocities


        // Update transport equation coefficients


        // Nonlinear iterations


            // Coupled sweep for RED nodes in +z direction


            // Coupled sweep for BLACK nodes in -z direction


        // Check residual

}