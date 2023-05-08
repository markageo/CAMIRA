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




// // Performs a single local update of block coupled equations
// class BlockSolver
// {

//     public:

//         void UpdateBlock(const intType i, const intType j, const intType k)
//         {
//             using enum Fields::ENUMDATA;
//             using enum TransportCoefficients::ENUMDATA;

//             // Precompute momentum RHS, divided by AP coefficients
//             floatType bu = ( m_fvCoeffs.Umom.B(i+1, j, k)

//                            - m_fvCoeffs.Umom.AU[n](i+1, j, k) * m_fields[U]( G(i+1, j+1, k  ) )
//                            - m_fvCoeffs.Umom.AU[e](i+1, j, k) * m_fields[U]( G(i+2, j  , k  ) )
//                            - m_fvCoeffs.Umom.AU[s](i+1, j, k) * m_fields[U]( G(i+1, j-1, k  ) )
//                            - m_fvCoeffs.Umom.AU[w](i+1, j, k) * m_fields[U]( G(i+1, j+1, k  ) )
//                            - m_fvCoeffs.Umom.AU[t](i+1, j, k) * m_fields[U]( G(i+1, j+1, k+1) )
//                            - m_fvCoeffs.Umom.AU[b](i+1, j, k) * m_fields[U]( G(i+1, j+1, k-1) )

//                            - m_fvCoeffs.Umom.AP[p](i+1, j, k) * m_fields[P]( G(i+1, j  , k  ) )
//                            - m_fvCoeffs.Umom.AP[e](i+1, j, k) * m_fields[P]( G(i+2, j  , k  ) )
//                            //- m_fvCoeffs.Umom.AP[w](i+1, j, k) * m_fields[P]( G(i  , j  , k  ) )

//                            ) / m_fvCoeffs.Umom.AU[p](i+1, j, k);



//             // Continuity for pressure
//             floatType bp = (m_fvCoeffs.Cont.B(i, j, k)
                            
//                           //- m_fvCoeffs.Cont.AV[n](i, j, k) * m_fields[V]( G(i, j+1, k) )
//                           - m_fvCoeffs.Cont.AV[p](i, j, k) * m_fields[V]( G(i, j  , k) )
//                           - m_fvCoeffs.Cont.AV[s](i, j, k) * m_fields[V]( G(i, j-1, k) )
                          
//                           //- m_fvCoeffs.Cont.AU[e](i, j, k) * m_fields[U]( G(i+1, j, k) )
//                           - m_fvCoeffs.Cont.AU[p](i, j, k) * m_fields[U]( G(i  , j, k) )
//                           - m_fvCoeffs.Cont.AU[w](i, j, k) * m_fields[U]( G(i-1, j, k) )
                          
//                           //- m_fvCoeffs.Cont.AW[t](i, j, k) * m_fields[W]( G(i, j, k+1) )
//                           - m_fvCoeffs.Cont.AW[p](i, j, k) * m_fields[W]( G(i, j, k  ) )
//                           - m_fvCoeffs.Cont.AW[b](i, j, k) * m_fields[W]( G(i, j, k-1) )
                          
//                           - m_fvCoeffs.Cont.AP[n](i, j, k) * m_fields[P]( G(i  , j+1, k  ) )
//                           - m_fvCoeffs.Cont.AP[e](i, j, k) * m_fields[P]( G(i+1, j  , k  ) )
//                           - m_fvCoeffs.Cont.AP[s](i, j, k) * m_fields[P]( G(i  , j-1, k  ) )
//                           - m_fvCoeffs.Cont.AP[w](i, j, k) * m_fields[P]( G(i-1, j  , k  ) )
//                           - m_fvCoeffs.Cont.AP[t](i, j, k) * m_fields[P]( G(i  , j  , k+1) )
//                           - m_fvCoeffs.Cont.AP[b](i, j, k) * m_fields[P]( G(i  , j  , k-1) )

//                           - m_fvCoeffs.Cont.AP[nn](i, j, k) * m_fields[P]( G(i  , j+2, k  ) )
//                           - m_fvCoeffs.Cont.AP[ee](i, j, k) * m_fields[P]( G(i+2, j  , k  ) )
//                           - m_fvCoeffs.Cont.AP[ss](i, j, k) * m_fields[P]( G(i  , j-2, k  ) )
//                           - m_fvCoeffs.Cont.AP[ww](i, j, k) * m_fields[P]( G(i-2, j  , k  ) )
//                           - m_fvCoeffs.Cont.AP[tt](i, j, k) * m_fields[P]( G(i  , j  , k+2) )
//                           - m_fvCoeffs.Cont.AP[bb](i, j, k) * m_fields[P]( G(i  , j  , k-2) )

//                           ) / m_fvCoeffs.Cont.AP[p](i, j, k);

//             // This only needs to be updated at linearisation
//             floatType K1 = m_fvCoeffs.Cont.AP[p](i, j, k)
//                          - m_fvCoeffs.Cont.AU[e](i, j, k) * m_fvCoeffs.Umom.AP[w](i+1, j  , k  ) / m_fvCoeffs.Umom.AU[p](i+1, j  , k  )
//                          - m_fvCoeffs.Cont.AV[n](i, j, k) * m_fvCoeffs.Vmom.AP[s](i  , j+1, k  ) / m_fvCoeffs.Vmom.AV[p](i  , j+1, k  )
//                          - m_fvCoeffs.Cont.AW[t](i, j, k) * m_fvCoeffs.Wmom.AP[b](i  , j  , k+1) / m_fvCoeffs.Umom.AW[p](i  , j  , k+1);


//             // Update P continuity
//             m_fields[P]( G(i, j, k) ) = ( bp
//                                         - m_fvCoeffs.Cont.AU[e](i, j, k) * bu
//                                         //- m_fvCoeffs.Cont.AV[n](i, j, k) * bv
//                                         //- m_fvCoeffs.Cont.AW[t](i, j, k) * bw
//                                         ) / K1;



//             // Update U momentum 
//             m_fields[U]( G(i+1, j, k) ) = bu 
//                                         - m_fvCoeffs.Umom.AP[w](i+1, j, k) * m_fields[P]( G(i  , j  , k  ) ) / m_fvCoeffs.Umom.AU[p](i+1, j, k);


//             // Update V momentum

//             // Update W momentum

//         }

//     private:
//         ArrayAllocator<Fields, array3D> &m_fields;
//         const FVCoefficients &m_fvCoeffs;
// };



// Solves line
class LineSolver
{
    public:

        LineSolver( ArrayAllocator<Fields, array3D> &fields,
                    FVCoefficients &fvCoeffs,
                    const InputData::LineSolverSettings &lineSolverSettings) :
        m_fields( fields ),
        m_fvCoeffs( fvCoeffs ),
        m_maxIterations( lineSolverSettings.maxIterations ),
        m_maxResiduals( lineSolverSettings.maxResiduals )
        {}

        void SolveLine(const intType j, const intType k)
        {

        }


    private:
        intType m_maxIterations;
        EnumVector<Fields, floatType> m_maxResiduals;
        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;

};


// Solves plane
class PlaneSolver
{
    public:

        PlaneSolver( ArrayAllocator<Fields, array3D> &fields,
                     FVCoefficients &fvCoeffs,
                     const InputData::PlaneSolverSettings &planeSolverSettings,
                     const InputData::LineSolverSettings &lineSolverSettings) :
        m_fields( fields ),
        m_fvCoeffs( fvCoeffs ),
        m_maxIterations( planeSolverSettings.maxIterations ),
        m_maxResiduals( planeSolverSettings.maxResiduals ),
        m_lineSolver( fields, fvCoeffs, lineSolverSettings )
        {}

        void SolvePlane(const intType k)
        {

        }


    private:
        intType m_maxIterations;
        EnumVector<Fields, floatType> m_maxResiduals;
        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;
        LineSolver m_lineSolver;
};




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