#include "Solver.h"
#include "InputProcessing.h"
#include "SimulationParameters.h"
#include "FiniteVolumeStructures.h"

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