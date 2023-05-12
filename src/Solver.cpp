#include "Types.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"
#include "Solver.h"


namespace CFD
{

namespace
{

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


// Helper class for staggered indexing
class StaggerIndexing 
{   
    using TC = TransportCoefficients::ENUMDATA;
    using F = Fields::ENUMDATA;

    public:

        // Pressure coupled in the momentum equation
        TC cPcoupled, cPleft, cPright;
        intType iPcoupled, iPleft, iPright;

        // Momentum coupled in the continuity equation
        TC cMcoupled, cMleft, cMright;
        intType iMcoupled, iMleft, iMright;

        constexpr StaggerIndexing(F field, TC staggeredCoeff)
        { 
           
            if        ( field == F::U ) {
                SetCompassU(staggeredCoeff);

            } else if ( field == F::V ) {
                SetCompassV(staggeredCoeff);

            } else if ( field == F::W ) {
                SetCompassW(staggeredCoeff);

            }
            SetIndex();
        } 

    private:

        constexpr void SetCompassU(TC staggeredCoeff)
        { 
            if        ( staggeredCoeff == TC::w ) {
                cPcoupled = TC::e;
                cPleft = TC::w;
                cPright = TC::p;

                cMcoupled = TC::w;
                cMleft = TC::p;
                cMright = TC::e;

            } else if ( staggeredCoeff == TC::p ) {
                cPcoupled = TC::p;
                cPleft = TC::w;
                cPright = TC::e;

                cMcoupled = TC::p;
                cMleft = TC::w;
                cMright = TC::e;

            } else if ( staggeredCoeff == TC::e) {
                cPcoupled = TC::w;
                cPleft = TC::p;
                cPright = TC::e;

                cMcoupled = TC::e;
                cMleft = TC::w;
                cMright = TC::p;

            }

        } 

        constexpr void SetCompassV(TC staggeredCoeff)
        { 
            if        ( staggeredCoeff == TC::s ) {
                cPcoupled = TC::n;
                cPleft = TC::s;
                cPright = TC::p;

                cMcoupled = TC::s;
                cMleft = TC::p;
                cMright = TC::n;

            } else if ( staggeredCoeff == TC::p ) {
                cPcoupled = TC::p;
                cPleft = TC::s;
                cPright = TC::n;

                cMcoupled = TC::p;
                cMleft = TC::s;
                cMright = TC::n;

            } else if ( staggeredCoeff == TC::n) {
                cPcoupled = TC::s;
                cPleft = TC::p;
                cPright = TC::n;

                cMcoupled = TC::n;
                cMleft = TC::s;
                cMright = TC::p;

            }
        } 

        constexpr void SetCompassW(TC staggeredCoeff)
        { 
            if        ( staggeredCoeff == TC::b ) {
                cPcoupled = TC::t;
                cPleft = TC::b;
                cPright = TC::p;

                cMcoupled = TC::b;
                cMleft = TC::p;
                cMright = TC::t;

            } else if ( staggeredCoeff == TC::p ) {
                cPcoupled = TC::p;
                cPleft = TC::b;
                cPright = TC::t;

                cMcoupled = TC::p;
                cMleft = TC::b;
                cMright = TC::t;

            } else if ( staggeredCoeff == TC::t) {
                cPcoupled = TC::b;
                cPleft = TC::p;
                cPright = TC::t;

                cMcoupled = TC::t;
                cMleft = TC::b;
                cMright = TC::p;
            }
        } 


        constexpr void SetIndex()
        {
            iPcoupled = CoeffIndex[cPcoupled];
            iPleft    = CoeffIndex[cPleft];
            iPright   = CoeffIndex[cPright];

            iMcoupled = CoeffIndex[cMcoupled];
            iMleft    = CoeffIndex[cMleft];
            iMright   = CoeffIndex[cMright];
        }  

};


// Performs a single local update of block coupled equations
class BlockSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    public:

        BlockSolver(ArrayAllocator<Fields, array3D> &fields, const FVCoefficients &fvCoeffs) :
            m_fields( fields ), m_fvCoeffs( fvCoeffs ) {};


        // Core function which updates the local coupled system. Templated by staggering direction.
        template<TC Ustag, TC Vstag, TC Wstag >
        EnumVector<Fields, floatType> UpdateBlock(const intType i, const intType j, const intType k)
        {
            using enum Fields::ENUMDATA;
            using enum TransportCoefficients::ENUMDATA;

            // Staggering must be valid
            static_assert( (Ustag == e) || (Ustag == w) || (Ustag == p), "Invalid U momentum staggering" );
            static_assert( (Vstag == n) || (Vstag == s) || (Vstag == p), "Invalid V momentum staggering" );
            static_assert( (Wstag == t) || (Wstag == b) || (Wstag == p), "Invalid W momentum staggering" );

            // Indexing variables to take care of staggering
            static constexpr StaggerIndexing sU( U, Ustag );
            static constexpr StaggerIndexing sV( V, Vstag );
            static constexpr StaggerIndexing sW( W, Wstag );
            AssertIndexing<sU, sV, sW>();


            // New values to return
            EnumVector<Fields, floatType> newBlock;

            // For indexing the staggered cells
            intType iU(i + sU.iMcoupled), jU(j               ), kU(k               ); // U momentum
            intType iV(i               ), jV(j + sV.iMcoupled), kV(k               ); // V momentum
            intType iW(i               ), jW(j               ), kW(k + sW.iMcoupled); // W momentum


            // Precompute momentum RHS divided by AP coefficients
            // U momentum
            floatType bU = ( m_fvCoeffs.Umom.B(iU, jU, kU)

                           - m_fvCoeffs.Umom.AU[n](iU, jU, kU) * m_fields[U]( G(iU  , jU+1, kU  ) )
                           - m_fvCoeffs.Umom.AU[e](iU, jU, kU) * m_fields[U]( G(iU+1, jU  , kU  ) )
                           - m_fvCoeffs.Umom.AU[s](iU, jU, kU) * m_fields[U]( G(iU  , jU-1, kU  ) )
                           - m_fvCoeffs.Umom.AU[w](iU, jU, kU) * m_fields[U]( G(iU-1, jU  , kU  ) )
                           - m_fvCoeffs.Umom.AU[t](iU, jU, kU) * m_fields[U]( G(iU  , jU  , kU+1) )
                           - m_fvCoeffs.Umom.AU[b](iU, jU, kU) * m_fields[U]( G(iU  , jU  , kU-1) )

                           - m_fvCoeffs.Umom.AP[sU.cPleft ](iU) * m_fields[P]( G(iU+sU.iPleft , jU  , kU  ) )
                           - m_fvCoeffs.Umom.AP[sU.cPright](iU) * m_fields[P]( G(iU+sU.iPright, jU  , kU  ) )

                           ) / m_fvCoeffs.Umom.AU[p](iU, jU, kU);


            // V momentum 
            floatType bV = ( m_fvCoeffs.Vmom.B(iV, jV, kV)

                           - m_fvCoeffs.Vmom.AV[n](iV, jV, kV) * m_fields[V]( G(iV  , jV+1, kV  ) )
                           - m_fvCoeffs.Vmom.AV[e](iV, jV, kV) * m_fields[V]( G(iV+1, jV  , kV  ) )
                           - m_fvCoeffs.Vmom.AV[s](iV, jV, kV) * m_fields[V]( G(iV  , jV-1, kV  ) )
                           - m_fvCoeffs.Vmom.AV[w](iV, jV, kV) * m_fields[V]( G(iV-1, jV  , kV  ) )
                           - m_fvCoeffs.Vmom.AV[t](iV, jV, kV) * m_fields[V]( G(iV  , jV  , kV+1) )
                           - m_fvCoeffs.Vmom.AV[b](iV, jV, kV) * m_fields[V]( G(iV  , jV  , kV-1) )

                           - m_fvCoeffs.Vmom.AP[sV.cPleft ](jV) * m_fields[P]( G(iV  , jV+sV.iPleft , kV  ) )
                           - m_fvCoeffs.Vmom.AP[sV.cPright](jV) * m_fields[P]( G(iV  , jV+sV.iPright, kV  ) )

                           ) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV);


            // W momentum
            floatType bW = ( m_fvCoeffs.Wmom.B(iW, jW, kW)

                           - m_fvCoeffs.Wmom.AW[n](iW, jW, kW) * m_fields[W]( G(iW  , jW+1, kW  ) )
                           - m_fvCoeffs.Wmom.AW[e](iW, jW, kW) * m_fields[W]( G(iW+1, jW  , kW  ) )
                           - m_fvCoeffs.Wmom.AW[s](iW, jW, kW) * m_fields[W]( G(iW  , jW-1, kW  ) )
                           - m_fvCoeffs.Wmom.AW[w](iW, jW, kW) * m_fields[W]( G(iW-1, jW  , kW  ) )
                           - m_fvCoeffs.Wmom.AW[t](iW, jW, kW) * m_fields[W]( G(iW  , jW  , kW+1) )
                           - m_fvCoeffs.Wmom.AW[b](iW, jW, kW) * m_fields[W]( G(iW  , jW  , kW-1) )

                           - m_fvCoeffs.Wmom.AP[sW.cPleft ](kW) * m_fields[P]( G(iW  , jW  , kW+sW.iPleft ) )
                           - m_fvCoeffs.Wmom.AP[sW.cPright](kW) * m_fields[P]( G(iW  , jW  , kW+sW.iPright) )

                           ) / m_fvCoeffs.Wmom.AW[p](iW, jW, kW);



            // Continuity for pressure
            floatType bP = (m_fvCoeffs.Cont.B(i, j, k)
                            
                          - m_fvCoeffs.Cont.AU[sU.cMleft ](i) * m_fields[U]( G(i+sU.iMleft , j, k) )
                          - m_fvCoeffs.Cont.AU[sU.cMright](i) * m_fields[U]( G(i+sU.iMright, j, k) )

                          - m_fvCoeffs.Cont.AV[sV.cMleft ](j) * m_fields[V]( G(i, j+sV.iMleft , k) )
                          - m_fvCoeffs.Cont.AV[sV.cMright](j) * m_fields[V]( G(i, j+sV.iMright, k) )
                          
                          - m_fvCoeffs.Cont.AW[sW.cMleft ](k) * m_fields[W]( G(i, j, k+sW.iMleft ) )
                          - m_fvCoeffs.Cont.AW[sW.cMright](k) * m_fields[W]( G(i, j, k+sW.iMright) )
                          
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
                        - m_fvCoeffs.Cont.AU[sU.cMcoupled](i) * m_fvCoeffs.Umom.AP[sU.cPcoupled](iU) / m_fvCoeffs.Umom.AU[p](iU, jU, kU)
                        - m_fvCoeffs.Cont.AV[sV.cMcoupled](j) * m_fvCoeffs.Vmom.AP[sV.cPcoupled](jV) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV)
                        - m_fvCoeffs.Cont.AW[sW.cMcoupled](k) * m_fvCoeffs.Wmom.AP[sW.cPcoupled](kW) / m_fvCoeffs.Umom.AW[p](iW, jW, kW);
            K = 1.0f / K;


            // Update P from continuity
            newBlock[P] = ( bP
                        - m_fvCoeffs.Cont.AU[sU.cMcoupled](i) * bU
                        - m_fvCoeffs.Cont.AV[sV.cMcoupled](j) * bV
                        - m_fvCoeffs.Cont.AW[sW.cMcoupled](k) * bW
                          ) * K;


            // Update U from momentum 
            newBlock[U]= bU 
                        - m_fvCoeffs.Umom.AP[sU.cPcoupled](iU) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Umom.AU[p](iU, jU, kU);


            // Update V from momentum
            newBlock[V] = bV 
                        - m_fvCoeffs.Vmom.AP[sV.cPcoupled](jV) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV);

            // Update W from momentum
            newBlock[W] = bW 
                        - m_fvCoeffs.Wmom.AP[sW.cPcoupled](kW) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Wmom.AW[p](iW, jW, kW);

            return newBlock;
        }


    private:

        const ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;

        // Static assert on all members of StaggerIndexing class to check compile time evaluation
        template< StaggerIndexing sU, StaggerIndexing sV, StaggerIndexing sW >
        constexpr void AssertIndexing()
        {
            // sU
            static_assert( sU.cPcoupled == TC::w || sU.cPcoupled == TC::p || sU.cPcoupled == TC::e);
            static_assert( sU.cPleft    == TC::w || sU.cPleft    == TC::p || sU.cPleft    == TC::e);
            static_assert( sU.cPright   == TC::w || sU.cPright   == TC::p || sU.cPright   == TC::e);
            static_assert( sU.cMcoupled == TC::w || sU.cMcoupled == TC::p || sU.cMcoupled == TC::e);
            static_assert( sU.cMleft    == TC::w || sU.cMleft    == TC::p || sU.cMleft    == TC::e);
            static_assert( sU.cMright   == TC::w || sU.cMright   == TC::p || sU.cMright   == TC::e);

            static_assert( sU.iPcoupled == -1    || sU.iPcoupled == 0     || sU.iPcoupled == 1);
            static_assert( sU.iPleft    == -1    || sU.iPleft    == 0     || sU.iPleft    == 1);
            static_assert( sU.iPright   == -1    || sU.iPright   == 0     || sU.iPright   == 1);
            static_assert( sU.iMcoupled == -1    || sU.iMcoupled == 0     || sU.iMcoupled == 1);
            static_assert( sU.iMleft    == -1    || sU.iMleft    == 0     || sU.iMleft    == 1);
            static_assert( sU.iMright   == -1    || sU.iMright   == 0     || sU.iMright   == 1);


            // sV
            static_assert( sV.cPcoupled == TC::s || sV.cPcoupled == TC::p || sV.cPcoupled == TC::n);
            static_assert( sV.cPleft    == TC::s || sV.cPleft    == TC::p || sV.cPleft    == TC::n);
            static_assert( sV.cPright   == TC::s || sV.cPright   == TC::p || sV.cPright   == TC::n);
            static_assert( sV.cMcoupled == TC::s || sV.cMcoupled == TC::p || sV.cMcoupled == TC::n);
            static_assert( sV.cMleft    == TC::s || sV.cMleft    == TC::p || sV.cMleft    == TC::n);
            static_assert( sV.cMright   == TC::s || sV.cMright   == TC::p || sV.cMright   == TC::n);

            static_assert( sV.iPcoupled == -1    || sV.iPcoupled == 0     || sV.iPcoupled == 1);
            static_assert( sV.iPleft    == -1    || sV.iPleft    == 0     || sV.iPleft    == 1);
            static_assert( sV.iPright   == -1    || sV.iPright   == 0     || sV.iPright   == 1);
            static_assert( sV.iMcoupled == -1    || sV.iMcoupled == 0     || sV.iMcoupled == 1);
            static_assert( sV.iMleft    == -1    || sV.iMleft    == 0     || sV.iMleft    == 1);
            static_assert( sV.iMright   == -1    || sV.iMright   == 0     || sV.iMright   == 1);


            // sW
            static_assert( sW.cPcoupled == TC::b || sW.cPcoupled == TC::p || sW.cPcoupled == TC::t);
            static_assert( sW.cPleft    == TC::b || sW.cPleft    == TC::p || sW.cPleft    == TC::t);
            static_assert( sW.cPright   == TC::b || sW.cPright   == TC::p || sW.cPright   == TC::t);
            static_assert( sW.cMcoupled == TC::b || sW.cMcoupled == TC::p || sW.cMcoupled == TC::t);
            static_assert( sW.cMleft    == TC::b || sW.cMleft    == TC::p || sW.cMleft    == TC::t);
            static_assert( sW.cMright   == TC::b || sW.cMright   == TC::p || sW.cMright   == TC::t);

            static_assert( sW.iPcoupled == -1    || sW.iPcoupled == 0     || sW.iPcoupled == 1);
            static_assert( sW.iPleft    == -1    || sW.iPleft    == 0     || sW.iPleft    == 1);
            static_assert( sW.iPright   == -1    || sW.iPright   == 0     || sW.iPright   == 1);
            static_assert( sW.iMcoupled == -1    || sW.iMcoupled == 0     || sW.iMcoupled == 1);
            static_assert( sW.iMleft    == -1    || sW.iMleft    == 0     || sW.iMleft    == 1);
            static_assert( sW.iMright   == -1    || sW.iMright   == 0     || sW.iMright   == 1);
        }
};





// Solves line
class LineSolver
{
    using TC = TransportCoefficients::ENUMDATA; 

    public:

        LineSolver( ArrayAllocator<Fields, array3D> &fields,
                    FVCoefficients &fvCoeffs,
                    const InputData::LineSolverSettings &lineSolverSettings) :
        m_fields( fields ),
        m_fvCoeffs( fvCoeffs ),
        m_maxIterations( lineSolverSettings.maxIterations ),
        m_maxResiduals( lineSolverSettings.maxResiduals ),
        m_blockSolver( fields, fvCoeffs )
        {}


        template< TC Vstag, TC Wstag >
        void SolveLine(const intType j, const intType k)
        {
            using enum TransportCoefficients::ENUMDATA;
            using enum Fields::ENUMDATA;
            intType ni = m_fields[Fields::ENUMDATA::U].dimension(Axis::ENUMDATA::X);

            // Temporary for storing new block update
            EnumVector<Fields, floatType> newBlock;
            TC Ustag = e;

            intType nIterations = 0;
            while ( nIterations < m_maxIterations ) 
            {
                // Forward sweep
                Ustag = e;
                for (intType i = 0; i != ni-1; i++) {

                    newBlock = m_blockSolver.UpdateBlock<Ustag, Vstag, Wstag>(i, j, k);
                    m_fields[U](i+CoeffIndex[Ustag], j, k) = newBlock[U];
                    m_fields[V](i, j+CoeffIndex[Vstag], k) = newBlock[V];
                    m_fields[W](i, j, k+CoeffIndex[Wstag]) = newBlock[W];
                    m_fields[P](i, j, k) = newBlock[P];
                    // ***Update with relaxation

                }

                // Backward sweep
                Ustag = w;
                for (intType i = ni-1; i != 0; i--) {

                    newBlock = m_blockSolver.UpdateBlock<Ustag, Vstag, Wstag>(i, j, k);
                    m_fields[U](i+CoeffIndex[Ustag], j, k) = newBlock[U];
                    m_fields[V](i, j+CoeffIndex[Vstag], k) = newBlock[V];
                    m_fields[W](i, j, k+CoeffIndex[Wstag]) = newBlock[W];
                    m_fields[P](i, j, k) = newBlock[P];
                    // ***Update with relaxation

                }

                nIterations++;
                // Check residual
                if ( MetLineResidualTolerence() ) {

                }
            }

        }


    private:

        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;
        intType m_maxIterations;
        EnumVector<Fields, floatType> m_maxResiduals;
        BlockSolver m_blockSolver;

        void UpdateLineResidual()
        {

        }

        bool MetLineResidualTolerence()
        {
            return false;
        }


};





// Solves plane
class PlaneSolver
{
    using TC = TransportCoefficients::ENUMDATA;

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

        template<TC Wstag>
        void SolvePlane(const intType k)
        {

        }


    private:
        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;
        intType m_maxIterations;
        EnumVector<Fields, floatType> m_maxResiduals;
        LineSolver m_lineSolver;
};





void SweepSolve( ArrayAllocator<CFD::Fields, CFD::array3D>  &fields, 
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
        
            // Forward sweep

            // Backward sweep


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


}   // end namespace CFD