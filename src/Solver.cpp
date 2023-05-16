#include "Types.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"
#include "Solver.h"

#include <type_traits>

namespace CFD
{

namespace
{

// Update the residual of the fields. Calculated as the L1 norm between current and previous iterations
// Template so different dimension arrays can be used
template< typename arrayType >
void UpdateResiduals( EnumVector<Fields, floatType> &residuals,
                      const ArrayAllocator<Fields, arrayType> &fields,
                      const ArrayAllocator<Fields, arrayType> &fieldsOld,
                      const EnumVector<Fields, floatType> &residualsInitial = EnumVector<Fields, floatType>({1.0f, 1.0f, 1.0f, 1.0f}))
{
    static_assert( std::is_same<arrayType, array1D>::value || 
                   std::is_same<arrayType, array2D>::value ||
                   std::is_same<arrayType, array3D>::value );

    Fields::ENUMDATA field;
    for (int f = 0; f != Fields::count; f++) {
        field = static_cast<Fields::ENUMDATA>(f);

        auto fieldDiff = fields[field] - fieldsOld[field];  // auto delays calculation
        array0D temp = fieldDiff.abs().mean();
        residuals[field] = temp(0) / residualsInitial[field];
    }
}


// Check if residual tolerence is met
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




// Static assert on all members of StaggerIndexing class to check compile time evaluation
template< StaggerIndexing sI, Fields::ENUMDATA field >
constexpr void AssertStaggerIndexing()
{
    using TC = TransportCoefficients::ENUMDATA;
    using F = Fields::ENUMDATA;

    if        ( field == F::U ) {
        static_assert( sI.cPcoupled == TC::w || sI.cPcoupled == TC::p || sI.cPcoupled == TC::e );
        static_assert( sI.cPleft    == TC::w || sI.cPleft    == TC::p || sI.cPleft    == TC::e );
        static_assert( sI.cPright   == TC::w || sI.cPright   == TC::p || sI.cPright   == TC::e );
        static_assert( sI.cMcoupled == TC::w || sI.cMcoupled == TC::p || sI.cMcoupled == TC::e );
        static_assert( sI.cMleft    == TC::w || sI.cMleft    == TC::p || sI.cMleft    == TC::e );
        static_assert( sI.cMright   == TC::w || sI.cMright   == TC::p || sI.cMright   == TC::e );

    } else if (field == F::V) {
        static_assert( sI.cPcoupled == TC::s || sI.cPcoupled == TC::p || sI.cPcoupled == TC::n );
        static_assert( sI.cPleft    == TC::s || sI.cPleft    == TC::p || sI.cPleft    == TC::n );
        static_assert( sI.cPright   == TC::s || sI.cPright   == TC::p || sI.cPright   == TC::n );
        static_assert( sI.cMcoupled == TC::s || sI.cMcoupled == TC::p || sI.cMcoupled == TC::n );
        static_assert( sI.cMleft    == TC::s || sI.cMleft    == TC::p || sI.cMleft    == TC::n );
        static_assert( sI.cMright   == TC::s || sI.cMright   == TC::p || sI.cMright   == TC::n );

    } else if (field == F::W) {
        static_assert( sI.cPcoupled == TC::b || sI.cPcoupled == TC::p || sI.cPcoupled == TC::t );
        static_assert( sI.cPleft    == TC::b || sI.cPleft    == TC::p || sI.cPleft    == TC::t );
        static_assert( sI.cPright   == TC::b || sI.cPright   == TC::p || sI.cPright   == TC::t );
        static_assert( sI.cMcoupled == TC::b || sI.cMcoupled == TC::p || sI.cMcoupled == TC::t );
        static_assert( sI.cMleft    == TC::b || sI.cMleft    == TC::p || sI.cMleft    == TC::t );
        static_assert( sI.cMright   == TC::b || sI.cMright   == TC::p || sI.cMright   == TC::t );

    }

    static_assert( sI.iPcoupled == -1    || sI.iPcoupled == 0     || sI.iPcoupled == 1 );
    static_assert( sI.iPleft    == -1    || sI.iPleft    == 0     || sI.iPleft    == 1 );
    static_assert( sI.iPright   == -1    || sI.iPright   == 0     || sI.iPright   == 1 );
    static_assert( sI.iMcoupled == -1    || sI.iMcoupled == 0     || sI.iMcoupled == 1 );
    static_assert( sI.iMleft    == -1    || sI.iMleft    == 0     || sI.iMleft    == 1 );
    static_assert( sI.iMright   == -1    || sI.iMright   == 0     || sI.iMright   == 1 );
}





// Performs a single local update of block coupled equations
class BlockSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    public:

        BlockSolver(ArrayAllocator<Fields, array3D> &fields, const FVCoefficients &fvCoeffs) :
            m_fields( fields ), 
            m_fvCoeffs( fvCoeffs ) 
        {};


        // Core function which updates the local coupled system. Templated by staggering direction.
        template<TC Ustag, TC Vstag, TC Wstag >
        void UpdateBlock(const intType i, const intType j, const intType k)
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
            AssertStaggerIndexing<sU, U>();
            AssertStaggerIndexing<sV, V>();
            AssertStaggerIndexing<sW, W>();

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
            m_fields[P]( G(i, j, k) ) = ( bP
                                        - m_fvCoeffs.Cont.AU[sU.cMcoupled](i) * bU
                                        - m_fvCoeffs.Cont.AV[sV.cMcoupled](j) * bV
                                        - m_fvCoeffs.Cont.AW[sW.cMcoupled](k) * bW
                                        ) * K;


            // Update U from momentum 
            m_fields[U]( G(iU, jU, kU) ) = bU 
                                         - m_fvCoeffs.Umom.AP[sU.cPcoupled](iU) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Umom.AU[p](iU, jU, kU);


            // Update V from momentum
            m_fields[V]( G(iV, jV, kV) ) = bV 
                                         - m_fvCoeffs.Vmom.AP[sV.cPcoupled](jV) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV);

            // Update W from momentum
            m_fields[W]( G(iW, jW, kW) ) = bW 
                                         - m_fvCoeffs.Wmom.AP[sW.cPcoupled](kW) * m_fields[P]( G(i, j, k) ) / m_fvCoeffs.Wmom.AW[p](iW, jW, kW);
        }


    private:

        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;
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
        m_relaxation( lineSolverSettings.relaxation ),
        m_blockSolver( fields, fvCoeffs )
        {}


        template< TC Vstag, TC Wstag >
        void SolveLine(const intType j, const intType k)
        {
            using enum TransportCoefficients::ENUMDATA;
            using enum Fields::ENUMDATA;
            using enum Axis::ENUMDATA;
            intType ni = m_fields[U].dimension(X) - 2*nGhost;

            // Staggering must be valid
            static_assert( (Wstag == n) || (Wstag == s) || (Wstag == p), "Invalid V momentum staggering" );
            static_assert( (Wstag == t) || (Wstag == b) || (Wstag == p), "Invalid W momentum staggering" );

            // Temporary for storing new block update
            EnumVector<Fields, array1D> oldLine( array1D( m_fields[U].dimension(X) ) );

            // Residuals 
            EnumVector<Fields, floatType> delta, residuals, residualsInitialInv;

            // Staggered indexing
            EnumVector<Fields, intType> iS( {0, 0, 0, 0} ), 
                                        jS( {j, j+coeffIndex[Vstag], j, j} ), 
                                        kS( {k, k, k+coeffIndex[Wstag], k} );

            intType nIterations = 0;
            while ( nIterations < m_maxIterations ) 
            {
                residuals = {0.0f, 0.0f, 0.0f, 0.0f};

                // Set the old line values
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { oldLine[field] = m_fields[field].chip(k, Z).chip(j, Y); } );


                // Update in place and relax
                for (intType i = 0; i != ni-1; i++) {   // Forward sweep

                    m_blockSolver.UpdateBlock<e, Vstag, Wstag>(i, j, k);

                    EnumFor<Fields>( [&] (Fields::ENUMDATA f) { iS[field] = i } );
                    iS[U] = i + coeffIndex[e];

                    EnumFor<Fields>( [&] (Fields::ENUMDATA f) 
                    { 
                         delta[f] = m_relaxation[f] * ( m_fields[f]( G(iS[f], jS[f], kS[f]) ) - oldLine[f]( G(iS[f]) ) );   
                         m_fields[f]( G(iS[f], jS[f], kS[f]) ) = oldLine[f]( G(iS[f]) ) + delta[f];                        
                         residuals[f] += abs( delta[f] );                                                                
                    } );
                }


                for (intType i = ni-1; i != 0; i--) {   // Backward sweep
                    
                    m_blockSolver.UpdateBlock<w, Vstag, Wstag>(i, j, k);

                    EnumFor<Fields>( [&] (Fields::ENUMDATA f) { iS[field] = i } );
                    iS[U] = i + coeffIndex[w];

                    EnumFor<Fields>( [&] (Fields::ENUMDATA f) 
                    { 
                         delta[f] = m_relaxation[f] * ( m_fields[f]( G(iS[f], jS[f], kS[f]) ) - oldLine[f]( G(iS[f]) ) );   
                         m_fields[f]( G(iS[f], jS[f], kS[f]) ) = oldLine[f]( G(iS[f]) ) + delta[f];                         
                         residuals[f] += abs( delta[f] );                                                                  
                    } );
                }


                // Normalise residuals
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] /= ni; } );
                if ( nIterations == 0 ) {
                    EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residualsInitialInv[field] = 1.0f / residuals[field]; } );
                }
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] *= residualsInitialInv[field];} );
                nIterations++;

                // Check residual
                if ( MetResidualTolerence( residuals, m_maxResiduals ) ) {
                    break;
                }

            }
        }

    private:

        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;
        intType m_maxIterations;
        EnumVector<Fields, floatType> m_maxResiduals;
        EnumVector<Fields, floatType> m_relaxation;
        BlockSolver m_blockSolver;
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
        m_relaxation( planeSolverSettings.relaxation ),
        m_lineSolver( fields, fvCoeffs, lineSolverSettings )
        {}


        template<TC Wstag>
        ArrayAllocator<Fields, array2D> SolvePlane(const intType k)
        {
            using enum Fields::ENUMDATA;
            using enum Axis::ENUMDATA;
            using enum TransportCoefficients::ENUMDATA;
            intType ni = m_fields[U].dimension(X) - 2*nGhost;
            intType nj = m_fields[U].dimension(Y) - 2*nGhost;

            // Staggering must be valid
            static_assert( (Wstag == t) || (Wstag == b) || (Wstag == p), "Invalid W momentum staggering" );

            // Old plane solution
            EnumVector<Fields, array2D> oldPlane( array2D( m_fields[U].dimension(X) , m_fields[U].dimension(Y) ) );

            // Residuals
            EnumVector<Fields, floatType> residuals, residualsInitialInv;
            EnumVector<Fields, array1D> delta( array1D( m_fields[U].dimension(X) ) );

            intType nIterations = 0;
            while ( nIterations < m_maxIterations ) 
            {
                residuals = {0.0f, 0.0f, 0.0f, 0.0f};

                // Set value of old plane
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { oldPlane[field] = m_fields[field].chip(k, Z); } );


                for (intType j = 0; j != nj-1; j++) {   // Forward sweep

                    m_lineSolver.SolveLine<n, Wstag>(j, k);

                    delta[U] = delta[U].constant( m_relaxation[U] ) * ( m_fields[U].chip(k, Z).chip(j  , Y) - oldPlane[U].chip(j, Y) );
                    delta[V] = delta[V].constant( m_relaxation[V] ) * ( m_fields[V].chip(k, Z).chip(j+1, Y) - oldPlane[V].chip(j+1, Y) );
                    delta[W] = delta[W].constant( m_relaxation[W] ) * ( m_fields[W].chip(k, Z).chip(j  , Y) - oldPlane[W].chip(j, Y) );
                    delta[P] = delta[P].constant( m_relaxation[P] ) * ( m_fields[P].chip(k, Z).chip(j  , Y) - oldPlane[P].chip(j, Y) );

                    m_fields[U].chip(k, Z).chip(j  , Y) = oldPlane[U].chip(j  , Y) + delta[U];
                    m_fields[V].chip(k, Z).chip(j+1, Y) = oldPlane[U].chip(j+1, Y) + delta[V];
                    m_fields[W].chip(k, Z).chip(j  , Y) = oldPlane[U].chip(j  , Y) + delta[W];
                    m_fields[P].chip(k, Z).chip(j  , Y) = oldPlane[U].chip(j  , Y) + delta[P];

                    EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] += static_cast<array0D>( delta[field].abs().sum() )(0); } );

                }

                for (intType j = nj-1; j != 0; j--) {   // Backward sweep

                    m_lineSolver.SolveLine<s, Wstag>(j, k);
                    
                    delta[U] = delta[U].constant( m_relaxation[U] ) * ( m_fields[U].chip(k, Z).chip(j  , Y) - oldPlane[U].chip(j, Y) );
                    delta[V] = delta[V].constant( m_relaxation[V] ) * ( m_fields[V].chip(k, Z).chip(j-1, Y) - oldPlane[V].chip(j-1, Y) );
                    delta[W] = delta[W].constant( m_relaxation[W] ) * ( m_fields[W].chip(k, Z).chip(j  , Y) - oldPlane[W].chip(j, Y) );
                    delta[P] = delta[P].constant( m_relaxation[P] ) * ( m_fields[P].chip(k, Z).chip(j  , Y) - oldPlane[P].chip(j, Y) );

                    m_fields[U].chip(k, Z).chip(j  , Y) = oldPlane[U].chip(j  , Y) + delta[U];
                    m_fields[V].chip(k, Z).chip(j-1, Y) = oldPlane[U].chip(j-1, Y) + delta[V];
                    m_fields[W].chip(k, Z).chip(j  , Y) = oldPlane[U].chip(j  , Y) + delta[W];
                    m_fields[P].chip(k, Z).chip(j  , Y) = oldPlane[U].chip(j  , Y) + delta[P];

                    EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] += static_cast<array0D>( delta[field].abs().sum() )(0); } );

                }

                // Normalise residuals
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] /= ni*nj; } );
                if ( nIterations == 0 ) {
                    EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residualsInitialInv[field] = 1.0f / residuals[field]; } );
                }
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] *= residualsInitialInv[field];} );
                nIterations++;

                // Check residual tolerence
                if ( MetResidualTolerence( residuals, m_maxResiduals ) ) {
                    break;
                }

            }

            return solutionPlane;
        }


    private:
        ArrayAllocator<Fields, array3D> &m_fields;
        const FVCoefficients &m_fvCoeffs;
        intType m_maxIterations;
        EnumVector<Fields, floatType> m_maxResiduals;
        EnumVector<Fields, floatType> m_relaxation;
        LineSolver m_lineSolver;
};





void SweepSolve( ArrayAllocator<CFD::Fields, CFD::array3D>  &fields, 
                 const Mesh &mesh, 
                 const InputData &inputData) 
{
    using enum TransportCoefficients::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum Fields::ENUMDATA;

    // Extract from input data
    const InputData::PlaneSweepSettings  planeSweepSettings  = inputData.planeSweepSettings;
    const InputData::PlaneSolverSettings planeSolverSettings = inputData.planeSolverSettings;
    const InputData::LineSolverSettings  lineSolverSettings  = inputData.lineSolverSettings;

    // Initialise 
    ArrayAllocator<Fields, array3D> faceVelocities = InitialiseFaceVelocities( mesh, fields, inputData );
    ArrayAllocator<Fields, array3D> fieldsOld( fields );
    FVCoefficients fvCoeffs = InitialiseFVCoefficients( mesh, fields, inputData );
    
    // Plane solver
    PlaneSolver planeSolver(fields, fvCoeffs, planeSolverSettings, lineSolverSettings);

    // Counters and residuals
    intType nOuterIterations, nInnerIterations;
    EnumVector<Fields, floatType> residualsInner, residualsOuter, residualsOuterInitial, residualsInnerInitial;
    std::vector< EnumVector<Fields, floatType> > residualsHistory( planeSweepSettings.maxOuterIterations );

    // Useful data
    intType nk = mesh.nCells[Z];
    ArrayAllocator<Fields, array2D> newPlane( {U, V, P, W}, {mesh.nCells(X), mesh.nCells(Y)} );

    // For slicing ghost cells
    Eigen::array<intType, 2> offsets = { nGhost, nGhost };
    Eigen::array<intType, 2> extents = { mesh.nCells(X), mesh.nCells(Y) };

    // Outer iterations
    nOuterIterations = 0;
    while ( nOuterIterations < planeSweepSettings.maxOuterIterations ) 
    {
        
        // Inner iterations
        nInnerIterations = 0;
        while ( nInnerIterations < planeSweepSettings.maxInnerIterations ) {
        
            // Forward sweep
            for (intType k = 0; k != nk-1; k++) {

                newPlane = planeSolver.SolvePlane<t>(k);
                fields[U].chip(k  , Z).slice(offsets, extents) = newPlane[U];
                fields[V].chip(k  , Z).slice(offsets, extents) = newPlane[V];
                fields[W].chip(k+1, Z).slice(offsets, extents) = newPlane[W];
                fields[P].chip(k  , Z).slice(offsets, extents) = newPlane[P];
                // ***Update with relaxation
            }

            // Backward sweep
            for (intType k = nk-1; k != 1; k--) {

                newPlane = planeSolver.SolvePlane<b>(k);
                fields[U].chip(k  , Z).slice(offsets, extents) = newPlane[U];
                fields[V].chip(k  , Z).slice(offsets, extents) = newPlane[V];
                fields[W].chip(k-1, Z).slice(offsets, extents) = newPlane[W];
                fields[P].chip(k  , Z).slice(offsets, extents) = newPlane[P];
                // ***Update with relaxation
            }


            // Update residuals
            if ( nInnerIterations == 0 ) {
                UpdateResiduals(residualsInnerInitial, fields, fieldsOld);
            }
            UpdateResiduals(residualsInner, fields, fieldsOld);
            nInnerIterations++;

            // Update fields
            fieldsOld = fields;

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
    }


    // Strip unused data
    residualsHistory.resize( nOuterIterations );
 
}


}   // end namespace CFD