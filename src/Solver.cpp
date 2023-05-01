#include "Types.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"
#include "Solver.h"


void CFD::SweepSolve(ArrayAllocator<CFD::Fields, CFD::array3D>  &fields, 
                     const Mesh &mesh, 
                     const InputData &inputData) 
{
    // Initialise
    ArrayAllocator<Fields, array3D> faceVelocities = InitialiseFaceVelocities( mesh, fields, inputData );
    FVCoefficients fvCoeffs = InitialiseFVCoefficients( mesh, fields, inputData );
    ArrayAllocator<Fields, array3D> fieldsOld = fields;

    
    // Outer iterations
    while ( true )      // If less than max number of sweeps
    {

        // Coupled sweep for RED nodes in +z direction

        // Coupled sweep for BLACK nodes in -z direction

        
        // Check residual
        if ( true ) {
            break;
        }

        // Update nonlinear coefficients
        if ( true ) {
            UpdateFaceVelocities( faceVelocities, mesh, fields, inputData );
            UpdateFVCoefficients( fvCoeffs, mesh, fields, inputData );
        }

        // Update fields
        fieldsOld = fields;

    }
        

}