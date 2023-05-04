#include "Types.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"
#include "Solver.h"



namespace
{

using namespace CFD;

// Update the residual of the fields. Calculated as the L1 norm between current and previous iterations
void UpdateResiduals( EnumVector<Fields, floatType> &residuals,
                     const ArrayAllocator<Fields, array3D> &fields,
                     const ArrayAllocator<Fields, array3D> &fieldsOld,
                     const EnumVector<Fields, floatType> &residualsInitial = EnumVector<Fields, floatType>({1.0, 1.0, 1.0, 1.0}))
{
    Fields::ENUMDATA field;
    for (int f = 0; f != Fields::count; f++) {
        field = static_cast<Fields::ENUMDATA>(f);

        auto fieldDiff = fields[field] - fieldsOld[field];  // auto delays calculation
        Eigen::Tensor<floatType, 0> temp = fieldDiff.abs().mean();
        residuals[field] = temp(0) / residualsInitial[field];
    }
}


// Check it residual tolerence is met
bool MetResidualTolerence( const EnumVector<Fields, floatType> &residuals,
                           const EnumVector<Fields, floatType> &residualsTarget)
{
    Fields::ENUMDATA field;
    for (int f = 0; f != Fields::count; f++) {
        field = static_cast<Fields::ENUMDATA>(f);

        if ( residuals[field] > residualsTarget[field] )
            return false;

    }
    return true;
}

}   // end anonymous namespace


void CFD::SweepSolve(ArrayAllocator<CFD::Fields, CFD::array3D>  &fields, 
                     const Mesh &mesh, 
                     const InputData &inputData) 
{
    // Extract from input data
    const InputData::PlaneSweepSettings  planeSweepSettings  = inputData.planeSweepSettings;
    const InputData::PlaneSolverSettings planeSolverSettings = inputData.planeSolverSettings;
    const InputData::LineSolverSettings  lineSolverSettings  = inputData.lineSolverSettings;

    // Initialise
    ArrayAllocator<Fields, array3D> faceVelocities = InitialiseFaceVelocities( mesh, fields, inputData );
    ArrayAllocator<Fields, array3D> fieldsOld( fields );
    FVCoefficients fvCoeffs = InitialiseFVCoefficients( mesh, fields, inputData );
    
    // Counters and residuals
    intType nOuterIterations = 0;
    intType nInnerIterations = 0;
    EnumVector<Fields, floatType> residualsInner, residualsOuter, residualsOuterInitial, residualsInnerInitial;
    std::vector< EnumVector<Fields, floatType> > residualsHistory( planeSweepSettings.maxOuterIterations );


    // Outer iterations
    while ( nOuterIterations < planeSweepSettings.maxOuterIterations ) 
    {
        
        // Inner iterations
        nInnerIterations = 0;
        while ( nInnerIterations < planeSweepSettings.maxInnerIterations ) {
        
            // Coupled sweep for RED nodes in +z direction

            // Coupled sweep for BLACK nodes in -z direction
            

            // Update residuals
            if ( nInnerIterations == 0 ) {
                UpdateResiduals(residualsInnerInitial, fields, fieldsOld);
            }
            UpdateResiduals(residualsInner, fields, fieldsOld);
            nInnerIterations++;

            // Check residual tolerence
            if ( MetResidualTolerence(residualsInner, planeSweepSettings.maxInnerResiduals) ) {
                break;
            }
        }
        

        // Update residuals
        if ( nOuterIterations == 0 ) {
            UpdateResiduals(residualsOuterInitial, fields, fieldsOld);
        }
        UpdateResiduals(residualsOuter, fields, fieldsOld, residualsOuterInitial);
        residualsHistory[nOuterIterations] = residualsOuter;
        nOuterIterations++;

        // Check residual tolerence
        if ( MetResidualTolerence(residualsOuter, planeSweepSettings.maxOuterResiduals) ) {
            break;
        }

        // Update nonlinear coefficients
        UpdateFaceVelocities( faceVelocities, mesh, fields, inputData );
        UpdateFVCoefficients( fvCoeffs, mesh, fields, inputData );

        // Update fields
        fieldsOld = fields;
    }
        
}