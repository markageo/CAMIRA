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




// Performs a single local update of block coupled equations
class BlockSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    public:

        template<TC Ustag, TC Vstag, TC Wstag >
        void UpdateBlock(const intType i, const intType j, const intType k)
        {
            using enum Fields::ENUMDATA;
            using enum TransportCoefficients::ENUMDATA;

            // Staggering must be valid
            static_assert( (Ustag == e) || (Ustag == w), "Invalid U momentum staggering" );
            static_assert( (Vstag == n) || (Vstag == s), "Invalid V momentum staggering" );
            static_assert( (Wstag == t) || (Wstag == b), "Invalid W momentum staggering" );

            // For indexing the staggered cells
            intType iU(i), jU(j), kU(k); // U momentum
            intType iV(i), jV(j), kV(k); // V momentum
            intType iW(i), jW(j), kW(k); // W momentum

            // U momentum 
            if        constexpr ( Ustag == e ) {
                iU++;
            } else if constexpr ( Ustag == w ) {
                iU--;
            }

            // V momentum
            if        constexpr ( Vstag == n ) {
                jV++;
            } else if constexpr ( Vstag == s) {
                jV--;
            }

            // W momentum
            if        constexpr ( Wstag == t ) {
                kW++;
            } else if constexpr ( Wstag == b ) {
                kW--;
            }


            // Precompute momentum RHS, divided by AP coefficients
            // U momentum
            floatType bU = ( m_fvCoeffs.Umom.B(iU, jU, kU)

                           - m_fvCoeffs.Umom.AU[n](iU, jU, kU) * m_fields[U]( G(iU  , jU+1, kU  ) )
                           - m_fvCoeffs.Umom.AU[e](iU, jU, kU) * m_fields[U]( G(iU+1, jU  , kU  ) )
                           - m_fvCoeffs.Umom.AU[s](iU, jU, kU) * m_fields[U]( G(iU  , jU-1, kU  ) )
                           - m_fvCoeffs.Umom.AU[w](iU, jU, kU) * m_fields[U]( G(iU-1, jU  , kU  ) )
                           - m_fvCoeffs.Umom.AU[t](iU, jU, kU) * m_fields[U]( G(iU  , jU  , kU+1) )
                           - m_fvCoeffs.Umom.AU[b](iU, jU, kU) * m_fields[U]( G(iU  , jU  , kU-1) )

                           - m_fvCoeffs.Umom.AP[p](iU, jU, kU) * m_fields[P]( G(iU  , jU  , kU  ) )
                           - m_fvCoeffs.Umom.AP[e](iU, jU, kU) * m_fields[P]( G(iU+1, jU  , kU  ) )
                           //- m_fvCoeffs.Umom.AP[w](iU, jU, kU) * m_fields[P]( G(iU-1, jU  , kU  ) )

                           ) / m_fvCoeffs.Umom.AU[p](iU, jU, kU);


            // V momentum 
            floatType bV = ( m_fvCoeffs.Vmom.B(iV, jV, kV)

                           - m_fvCoeffs.Vmom.AV[n](iV, jV, kV) * m_fields[V]( G(iV  , jV+1, kV  ) )
                           - m_fvCoeffs.Vmom.AV[e](iV, jV, kV) * m_fields[V]( G(iV+1, jV  , kV  ) )
                           - m_fvCoeffs.Vmom.AV[s](iV, jV, kV) * m_fields[V]( G(iV  , jV-1, kV  ) )
                           - m_fvCoeffs.Vmom.AV[w](iV, jV, kV) * m_fields[V]( G(iV-1, jV  , kV  ) )
                           - m_fvCoeffs.Vmom.AV[t](iV, jV, kV) * m_fields[V]( G(iV  , jV  , kV+1) )
                           - m_fvCoeffs.Vmom.AV[b](iV, jV, kV) * m_fields[V]( G(iV  , jV  , kV-1) )

                           - m_fvCoeffs.Vmom.AP[p](iV, jV, kV) * m_fields[P]( G(iV  , jV  , kV  ) )
                           - m_fvCoeffs.Vmom.AP[n](iV, jV, kV) * m_fields[P]( G(iV  , jV+1, kV  ) )
                           //- m_fvCoeffs.Umom.AP[s](iV, jV, kV) * m_fields[P]( G(iV  , jV-1, kV  ) )

                           ) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV);


            // W momentum
            floatType bW = ( m_fvCoeffs.Wmom.B(iW, jW, kW)

                           - m_fvCoeffs.Wmom.AW[n](iW, jW, kW) * m_fields[W]( G(iW  , jW+1, kW  ) )
                           - m_fvCoeffs.Wmom.AW[e](iW, jW, kW) * m_fields[W]( G(iW+1, jW  , kW  ) )
                           - m_fvCoeffs.Wmom.AW[s](iW, jW, kW) * m_fields[W]( G(iW  , jW-1, kW  ) )
                           - m_fvCoeffs.Wmom.AW[w](iW, jW, kW) * m_fields[W]( G(iW-1, jW  , kW  ) )
                           - m_fvCoeffs.Wmom.AW[t](iW, jW, kW) * m_fields[W]( G(iW  , jW  , kW+1) )
                           - m_fvCoeffs.Wmom.AW[b](iW, jW, kW) * m_fields[W]( G(iW  , jW  , kW-1) )

                           - m_fvCoeffs.Wmom.AP[p](iW, jW, kW) * m_fields[P]( G(iW  , jW  , kW  ) )
                           - m_fvCoeffs.Wmom.AP[t](iW, jW, kW) * m_fields[P]( G(iW  , jW  , kW  ) )
                           //- m_fvCoeffs.Wmom.AP[b](iW, jW, kW) * m_fields[P]( G(iW  , jW  , kW  ) )

                           ) / m_fvCoeffs.Wmom.AW[p](iW, jW, kW);



            // Continuity for pressure
            floatType bP = (m_fvCoeffs.Cont.B(i, j, k)
                            
                          //- m_fvCoeffs.Cont.AV[n](i, j, k) * m_fields[V]( G(i, j+1, k) )
                          - m_fvCoeffs.Cont.AV[p](i, j, k) * m_fields[V]( G(i, j  , k) )
                          - m_fvCoeffs.Cont.AV[s](i, j, k) * m_fields[V]( G(i, j-1, k) )
                          
                          //- m_fvCoeffs.Cont.AU[e](i, j, k) * m_fields[U]( G(i+1, j, k) )
                          - m_fvCoeffs.Cont.AU[p](i, j, k) * m_fields[U]( G(i  , j, k) )
                          - m_fvCoeffs.Cont.AU[w](i, j, k) * m_fields[U]( G(i-1, j, k) )
                          
                          //- m_fvCoeffs.Cont.AW[t](i, j, k) * m_fields[W]( G(i, j, k+1) )
                          - m_fvCoeffs.Cont.AW[p](i, j, k) * m_fields[W]( G(i, j, k  ) )
                          - m_fvCoeffs.Cont.AW[b](i, j, k) * m_fields[W]( G(i, j, k-1) )
                          
                          - m_fvCoeffs.Cont.AP[n](i, j, k) * m_fields[P]( G(i  , j+1, k  ) )
                          - m_fvCoeffs.Cont.AP[e](i, j, k) * m_fields[P]( G(i+1, j  , k  ) )
                          - m_fvCoeffs.Cont.AP[s](i, j, k) * m_fields[P]( G(i  , j-1, k  ) )
                          - m_fvCoeffs.Cont.AP[w](i, j, k) * m_fields[P]( G(i-1, j  , k  ) )
                          - m_fvCoeffs.Cont.AP[t](i, j, k) * m_fields[P]( G(i  , j  , k+1) )
                          - m_fvCoeffs.Cont.AP[b](i, j, k) * m_fields[P]( G(i  , j  , k-1) )

                          - m_fvCoeffs.Cont.AP[nn](i, j, k) * m_fields[P]( G(i  , j+2, k  ) )
                          - m_fvCoeffs.Cont.AP[ee](i, j, k) * m_fields[P]( G(i+2, j  , k  ) )
                          - m_fvCoeffs.Cont.AP[ss](i, j, k) * m_fields[P]( G(i  , j-2, k  ) )
                          - m_fvCoeffs.Cont.AP[ww](i, j, k) * m_fields[P]( G(i-2, j  , k  ) )
                          - m_fvCoeffs.Cont.AP[tt](i, j, k) * m_fields[P]( G(i  , j  , k+2) )
                          - m_fvCoeffs.Cont.AP[bb](i, j, k) * m_fields[P]( G(i  , j  , k-2) )

                          ) / m_fvCoeffs.Cont.AP[p](i, j, k);

            // This only needs to be updated at linearisation
            floatType K = m_fvCoeffs.Cont.AP[p](i, j, k)
                        - m_fvCoeffs.Cont.AU[e](i, j, k) * m_fvCoeffs.Umom.AP[w](iU, jU, kU) / m_fvCoeffs.Umom.AU[p](iU, jU, kU)
                        - m_fvCoeffs.Cont.AV[n](i, j, k) * m_fvCoeffs.Vmom.AP[s](iV, jV, kV) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV)
                        - m_fvCoeffs.Cont.AW[t](i, j, k) * m_fvCoeffs.Wmom.AP[b](iW, jW, kW) / m_fvCoeffs.Umom.AW[p](iW, jW, kW);
            K = 1.0 / K;


            // Update P continuity
            m_fields[P]( G(i, j, k) ) = ( bP
                                        - m_fvCoeffs.Cont.AU[e](i, j, k) * bU
                                        - m_fvCoeffs.Cont.AV[n](i, j, k) * bV
                                        - m_fvCoeffs.Cont.AW[t](i, j, k) * bW
                                        ) * K;



            // Update U momentum 
            m_fields[U]( G(iU, jU, kU) ) = bU 
                                         - m_fvCoeffs.Umom.AP[w](iU, jU, kU) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Umom.AU[p](iU, jU, kU);


            // Update V momentum
            m_fields[V]( G(iV, jV, kV) ) = bV 
                                         - m_fvCoeffs.Vmom.AP[s](iV, jV, kV) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV);

            // Update W momentum
            m_fields[W]( G(iW, jW, kW) ) = bW 
                                         - m_fvCoeffs.Wmom.AP[b](iW, jW, kW) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Wmom.AW[p](iW, jW, kW);

        }




    private:
        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;
};



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