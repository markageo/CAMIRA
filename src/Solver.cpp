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

        BlockSolver(ArrayAllocator<Fields, array3D> &fields, const FVCoefficients &fvCoeffs) :
            m_fields( fields ), m_fvCoeffs( fvCoeffs ) {};


        template<TC Ustag, TC Vstag, TC Wstag >
        void UpdateBlock(const intType i, const intType j, const intType k)
        {
            using enum Fields::ENUMDATA;
            using enum TransportCoefficients::ENUMDATA;

            // Staggering must be valid
            static_assert( (Ustag == e) || (Ustag == w), "Invalid U momentum staggering" );
            static_assert( (Vstag == n) || (Vstag == s), "Invalid V momentum staggering" );
            static_assert( (Wstag == t) || (Wstag == b), "Invalid W momentum staggering" );

            static constexpr intType iUoffset = Ustaggering(Ustag), 
                                     jVoffset = Vstaggering(Vstag), 
                                     kWoffset = Wstaggering(Wstag);

            // Forces compile time evaluation
            static_assert( (iUoffset == 1) || (iUoffset == -1));        
            static_assert( (jVoffset == 1) || (jVoffset == -1));  
            static_assert( (kWoffset == 1) || (kWoffset == -1));      

            // For indexing the staggered cells
            intType iU(i + iUoffset), jU(j           ), kU(k           ); // U momentum
            intType iV(i           ), jV(j + jVoffset), kV(k           ); // V momentum
            intType iW(i           ), jW(j           ), kW(k + kWoffset); // W momentum


            // Precompute momentum RHS, divided by AP coefficients
            // U momentum
            floatType bU = ( m_fvCoeffs.Umom.B(iU, jU, kU)

                           - m_fvCoeffs.Umom.AU[n](iU, jU, kU) * m_fields[U]( G(iU  , jU+1, kU  ) )
                           - m_fvCoeffs.Umom.AU[e](iU, jU, kU) * m_fields[U]( G(iU+1, jU  , kU  ) )
                           - m_fvCoeffs.Umom.AU[s](iU, jU, kU) * m_fields[U]( G(iU  , jU-1, kU  ) )
                           - m_fvCoeffs.Umom.AU[w](iU, jU, kU) * m_fields[U]( G(iU-1, jU  , kU  ) )
                           - m_fvCoeffs.Umom.AU[t](iU, jU, kU) * m_fields[U]( G(iU  , jU  , kU+1) )
                           - m_fvCoeffs.Umom.AU[b](iU, jU, kU) * m_fields[U]( G(iU  , jU  , kU-1) )

                           - m_fvCoeffs.Umom.AP[p](iU) * m_fields[P]( G(iU  , jU  , kU  ) )
                           - m_fvCoeffs.Umom.AP[e](iU) * m_fields[P]( G(iU+1, jU  , kU  ) )
                           //- m_fvCoeffs.Umom.AP[w](iU) * m_fields[P]( G(iU-1, jU  , kU  ) )

                           ) / m_fvCoeffs.Umom.AU[p](iU, jU, kU);


            // V momentum 
            floatType bV = ( m_fvCoeffs.Vmom.B(iV, jV, kV)

                           - m_fvCoeffs.Vmom.AV[n](iV, jV, kV) * m_fields[V]( G(iV  , jV+1, kV  ) )
                           - m_fvCoeffs.Vmom.AV[e](iV, jV, kV) * m_fields[V]( G(iV+1, jV  , kV  ) )
                           - m_fvCoeffs.Vmom.AV[s](iV, jV, kV) * m_fields[V]( G(iV  , jV-1, kV  ) )
                           - m_fvCoeffs.Vmom.AV[w](iV, jV, kV) * m_fields[V]( G(iV-1, jV  , kV  ) )
                           - m_fvCoeffs.Vmom.AV[t](iV, jV, kV) * m_fields[V]( G(iV  , jV  , kV+1) )
                           - m_fvCoeffs.Vmom.AV[b](iV, jV, kV) * m_fields[V]( G(iV  , jV  , kV-1) )

                           - m_fvCoeffs.Vmom.AP[p](jV) * m_fields[P]( G(iV  , jV  , kV  ) )
                           - m_fvCoeffs.Vmom.AP[n](jV) * m_fields[P]( G(iV  , jV+1, kV  ) )
                           //- m_fvCoeffs.Umom.AP[s](jV) * m_fields[P]( G(iV  , jV-1, kV  ) )

                           ) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV);


            // W momentum
            floatType bW = ( m_fvCoeffs.Wmom.B(iW, jW, kW)

                           - m_fvCoeffs.Wmom.AW[n](iW, jW, kW) * m_fields[W]( G(iW  , jW+1, kW  ) )
                           - m_fvCoeffs.Wmom.AW[e](iW, jW, kW) * m_fields[W]( G(iW+1, jW  , kW  ) )
                           - m_fvCoeffs.Wmom.AW[s](iW, jW, kW) * m_fields[W]( G(iW  , jW-1, kW  ) )
                           - m_fvCoeffs.Wmom.AW[w](iW, jW, kW) * m_fields[W]( G(iW-1, jW  , kW  ) )
                           - m_fvCoeffs.Wmom.AW[t](iW, jW, kW) * m_fields[W]( G(iW  , jW  , kW+1) )
                           - m_fvCoeffs.Wmom.AW[b](iW, jW, kW) * m_fields[W]( G(iW  , jW  , kW-1) )

                           - m_fvCoeffs.Wmom.AP[p](kW) * m_fields[P]( G(iW  , jW  , kW  ) )
                           - m_fvCoeffs.Wmom.AP[t](kW) * m_fields[P]( G(iW  , jW  , kW+1) )
                           //- m_fvCoeffs.Wmom.AP[b](kW) * m_fields[P]( G(iW  , jW  , kW-1) )

                           ) / m_fvCoeffs.Wmom.AW[p](iW, jW, kW);



            // Continuity for pressure
            floatType bP = (m_fvCoeffs.Cont.B(i, j, k)
                            
                          //- m_fvCoeffs.Cont.AV[n](j) * m_fields[V]( G(i, j+1, k) )
                          - m_fvCoeffs.Cont.AV[p](j) * m_fields[V]( G(i, j  , k) )
                          - m_fvCoeffs.Cont.AV[s](j) * m_fields[V]( G(i, j-1, k) )
                          
                          //- m_fvCoeffs.Cont.AU[e](i) * m_fields[U]( G(i+1, j, k) )
                          - m_fvCoeffs.Cont.AU[p](i) * m_fields[U]( G(i  , j, k) )
                          - m_fvCoeffs.Cont.AU[w](i) * m_fields[U]( G(i-1, j, k) )
                          
                          //- m_fvCoeffs.Cont.AW[t](k) * m_fields[W]( G(i, j, k+1) )
                          - m_fvCoeffs.Cont.AW[p](k) * m_fields[W]( G(i, j, k  ) )
                          - m_fvCoeffs.Cont.AW[b](k) * m_fields[W]( G(i, j, k-1) )
                          
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
                        - m_fvCoeffs.Cont.AU[e](i) * m_fvCoeffs.Umom.AP[w](iU) / m_fvCoeffs.Umom.AU[p](iU, jU, kU)
                        - m_fvCoeffs.Cont.AV[n](j) * m_fvCoeffs.Vmom.AP[s](jV) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV)
                        - m_fvCoeffs.Cont.AW[t](k) * m_fvCoeffs.Wmom.AP[b](kW) / m_fvCoeffs.Umom.AW[p](iW, jW, kW);
            K = 1.0 / K;


            // Update P from continuity
            m_fields[P]( G(i, j, k) ) = ( bP
                                        - m_fvCoeffs.Cont.AU[e](i) * bU
                                        - m_fvCoeffs.Cont.AV[n](j) * bV
                                        - m_fvCoeffs.Cont.AW[t](k) * bW
                                        ) * K;


            // Update U from momentum 
            m_fields[U]( G(iU, jU, kU) ) = bU 
                                         - m_fvCoeffs.Umom.AP[w](iU) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Umom.AU[p](iU, jU, kU);


            // Update V from momentum
            m_fields[V]( G(iV, jV, kV) ) = bV 
                                         - m_fvCoeffs.Vmom.AP[s](jV) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV);

            // Update W from momentum
            m_fields[W]( G(iW, jW, kW) ) = bW 
                                         - m_fvCoeffs.Wmom.AP[b](kW) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Wmom.AW[p](iW, jW, kW);

        }


    private:
        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;

        // Return the index by which the momentum is staggered
        static constexpr intType Ustaggering(TC coeff)
        {
            if        ( coeff == TC::e ) {
                return 1;
            } else if ( coeff ==  TC::w ) {
                return -1;
            }
        }

        static constexpr intType Vstaggering(TC coeff)
        {
            if        ( coeff == TC::n ) {
                return 1;
            } else if ( coeff ==  TC::s ) {
                return -1;
            }
        }

        static constexpr intType Wstaggering(TC coeff)
        {
            if        ( coeff == TC::t ) {
                return 1;
            } else if ( coeff ==  TC::b ) {
                return -1;
            }
        }
    
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
            BlockSolver blockSolver(m_fields, m_fvCoeffs);
            using enum TransportCoefficients::ENUMDATA;
            blockSolver.UpdateBlock<e, n, t>(1, 2, 3);

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