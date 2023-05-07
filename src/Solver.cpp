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



// class PlaneSolver
// {
//     public:

//         PlaneSolver( ArrayAllocator<Fields, array3D> &fields,
//                      FVCoefficients &fvCoeffs,
//                      const InputData::PlaneSolverSettings &planeSolverSettings,
//                      const InputData::LineSolverSettings &lineSolverSettings) :
//         m_fields( fields ),
//         m_fvCoeffs( fvCoeffs ),
//         m_maxIterations( planeSolverSettings.maxIterations ),
//         m_maxResiduals( planeSolverSettings.maxResiduals ),
//         m_lineSolver( fields, fvCoeffs, lineSolverSettings )
//         {}

//         void SolvePlaneForwardStagger(const intType k)
//         {
//             using enum Axis::ENUMDATA;
//             intType nLines = m_fields[Fields::ENUMDATA::P].dimension(Z) - 2*nGhost;    // All fields should have the same dimensions

//             // Forwards sweep RED cells
//             for (intType j = 0; j != nLines-1; j = j + 2) {
//                 m_lineSolver.SolveLineForwardStagger(j, k);
//             }


//             // Backwards sweep BLACK cells
//             for (intType j = nLines-1; j != 0; j = j - 2) {
//                 m_lineSolver.SolveLineBackwardStagger(j, k);
//             }
//         }

//         void SolvePlaneBackwardStagger(const intType k)
//         {
//             using enum Axis::ENUMDATA;
//             intType nLines = m_fields[Fields::ENUMDATA::P].dimension(Z) - 2*nGhost;    // All fields should have the same dimensions

//             // Forwards sweep RED cells
//             for (intType j = 0; j != nLines-1; j = j + 2) {
                
                

//             }


//             // Backwards sweep BLACK cells
//             for (intType j = nLines-1; j != 0; j = j - 2) {
                


//             }
//         }

//     private:
//         intType m_maxIterations;
//         EnumVector<Fields, floatType> m_maxResiduals;
//         ArrayAllocator<Fields, array3D> &m_fields;
//         const FVCoefficients &m_fvCoeffs;
//         LineSolver m_lineSolver;
// };



// class LineSolver
// {
//     public:

//         LineSolver( ArrayAllocator<Fields, array3D> &fields,
//                     FVCoefficients &fvCoeffs,
//                     const InputData::LineSolverSettings &lineSolverSettings) :
//         m_fields( fields ),
//         m_fvCoeffs( fvCoeffs ),
//         m_maxIterations( lineSolverSettings.maxIterations ),
//         m_maxResiduals( lineSolverSettings.maxResiduals )
//         {}

//         void SolveLineForwardStagger(const intType j, const intType k)
//         {

//         }

//         void SolveLineBackwardStagger(const intType j, const intType k)
//         {

//         }

//     private:
//         intType m_maxIterations;
//         EnumVector<Fields, floatType> m_maxResiduals;
//         ArrayAllocator<Fields, array3D> &m_fields;
//         const FVCoefficients &m_fvCoeffs;

// };



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
    intType nOuterIterations, nInnerIterations;
    EnumVector<Fields, floatType> residualsInner, residualsOuter, residualsOuterInitial, residualsInnerInitial;
    std::vector< EnumVector<Fields, floatType> > residualsHistory( planeSweepSettings.maxOuterIterations );


    // Outer iterations
    nOuterIterations = 0;
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


    // Strip unused data
    residualsHistory.resize( nOuterIterations );

        
}