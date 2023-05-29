#include "Types.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"
#include "Solver.h"

#include "Utils.h"

#include <type_traits>
#include <iostream>

namespace CFD
{

namespace
{

// Calculated as the L1 norm of the difference between two arrays.
template< typename arrayType, typename enumStruct >
void L1ArrayDiff( EnumVector<enumStruct, floatType> &result,
                  const ArrayAllocator<enumStruct, arrayType> &array1,
                  const ArrayAllocator<enumStruct, arrayType> &array2)
{
    static_assert( std::is_same<arrayType, array1D>::value || 
                   std::is_same<arrayType, array2D>::value ||
                   std::is_same<arrayType, array3D>::value );

    EnumFor<enumStruct>( [&] (typename enumStruct::ENUMDATA enumName) { 
        auto fieldDiff = array1[enumName] - array2[enumName];  // auto lazily evaluates
        result[enumName] = static_cast<array0D>( fieldDiff.abs().mean() )(0);
    } );
}


// Turn the residual into a relative residual
void RelativeResidual( EnumVector<Fields, floatType> &residuals,
                       EnumVector<Fields, floatType> &residualsInitialInv,
                       const intType nIterations )
{
    if ( nIterations == 0 ) {
        EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residualsInitialInv[field] = 1.0f / residuals[field]; } );
    }
    EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] *= residualsInitialInv[field];} );
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
                cPleft    = TC::w;
                cPright   = TC::p;

                cMcoupled = TC::w;
                cMleft    = TC::p;
                cMright   = TC::e;

            } else if ( staggeredCoeff == TC::p ) {
                cPcoupled = TC::p;
                cPleft    = TC::w;
                cPright   = TC::e;

                cMcoupled = TC::p;
                cMleft    = TC::w;
                cMright   = TC::e;

            } else if ( staggeredCoeff == TC::e) {
                cPcoupled = TC::w;
                cPleft    = TC::p;
                cPright   = TC::e;

                cMcoupled = TC::e;
                cMleft    = TC::w;
                cMright   = TC::p;

            }

        } 

        constexpr void SetCompassV(TC staggeredCoeff)
        { 
            if        ( staggeredCoeff == TC::s ) {
                cPcoupled = TC::n;
                cPleft    = TC::s;
                cPright   = TC::p;

                cMcoupled = TC::s;
                cMleft    = TC::p;
                cMright   = TC::n;

            } else if ( staggeredCoeff == TC::p ) {
                cPcoupled = TC::p;
                cPleft    = TC::s;
                cPright   = TC::n;

                cMcoupled = TC::p;
                cMleft    = TC::s;
                cMright   = TC::n;

            } else if ( staggeredCoeff == TC::n) {
                cPcoupled = TC::s;
                cPleft   = TC::p;
                cPright  = TC::n;

                cMcoupled = TC::n;
                cMleft    = TC::s;
                cMright   = TC::p;

            }
        } 

        constexpr void SetCompassW(TC staggeredCoeff)
        { 
            if        ( staggeredCoeff == TC::b ) {
                cPcoupled = TC::t;
                cPleft    = TC::b;
                cPright   = TC::p;

                cMcoupled = TC::b;
                cMleft    = TC::p;
                cMright   = TC::t;

            } else if ( staggeredCoeff == TC::p ) {
                cPcoupled = TC::p;
                cPleft    = TC::b;
                cPright   = TC::t;

                cMcoupled = TC::p;
                cMleft    = TC::b;
                cMright   = TC::t;

            } else if ( staggeredCoeff == TC::t) {
                cPcoupled = TC::b;
                cPleft    = TC::p;
                cPright   = TC::t;

                cMcoupled = TC::t;
                cMleft    = TC::b;
                cMright   = TC::p;
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
consteval void AssertStaggerIndexing()
{
    using TC = TransportCoefficients::ENUMDATA;
    using F = Fields::ENUMDATA;

    if constexpr        ( field == F::U ) {

        static_assert( sI.cPcoupled == TC::w || sI.cPcoupled == TC::p || sI.cPcoupled == TC::e );
        static_assert( sI.cPleft    == TC::w || sI.cPleft    == TC::p || sI.cPleft    == TC::e );
        static_assert( sI.cPright   == TC::w || sI.cPright   == TC::p || sI.cPright   == TC::e );
        static_assert( sI.cMcoupled == TC::w || sI.cMcoupled == TC::p || sI.cMcoupled == TC::e );
        static_assert( sI.cMleft    == TC::w || sI.cMleft    == TC::p || sI.cMleft    == TC::e );
        static_assert( sI.cMright   == TC::w || sI.cMright   == TC::p || sI.cMright   == TC::e );

    } else if constexpr (field == F::V) {

        static_assert( sI.cPcoupled == TC::s || sI.cPcoupled == TC::p || sI.cPcoupled == TC::n );
        static_assert( sI.cPleft    == TC::s || sI.cPleft    == TC::p || sI.cPleft    == TC::n );
        static_assert( sI.cPright   == TC::s || sI.cPright   == TC::p || sI.cPright   == TC::n );
        static_assert( sI.cMcoupled == TC::s || sI.cMcoupled == TC::p || sI.cMcoupled == TC::n );
        static_assert( sI.cMleft    == TC::s || sI.cMleft    == TC::p || sI.cMleft    == TC::n );
        static_assert( sI.cMright   == TC::s || sI.cMright   == TC::p || sI.cMright   == TC::n );

    } else if constexpr (field == F::W) {

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
            floatType bP = m_fvCoeffs.Cont.B(i, j, k)
                            
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
                          - m_fvCoeffs.Cont.AP[bb](i, j, k) * m_fields[P]( G(i  , j  , k-2) );


            // This only needs to be updated at linearisation
            floatType K = m_fvCoeffs.Cont.AP[p](i, j, k)
                        - m_fvCoeffs.Cont.AU[sU.cMcoupled](i) * m_fvCoeffs.Umom.AP[sU.cPcoupled](iU) / m_fvCoeffs.Umom.AU[p](iU, jU, kU)
                        - m_fvCoeffs.Cont.AV[sV.cMcoupled](j) * m_fvCoeffs.Vmom.AP[sV.cPcoupled](jV) / m_fvCoeffs.Vmom.AV[p](iV, jV, kV)
                        - m_fvCoeffs.Cont.AW[sW.cMcoupled](k) * m_fvCoeffs.Wmom.AP[sW.cPcoupled](kW) / m_fvCoeffs.Wmom.AW[p](iW, jW, kW);


            // Update P from continuity
            m_fields[P]( G(i, j, k) ) = ( bP
                                        - m_fvCoeffs.Cont.AU[sU.cMcoupled](i) * bU
                                        - m_fvCoeffs.Cont.AV[sV.cMcoupled](j) * bV
                                        - m_fvCoeffs.Cont.AW[sW.cMcoupled](k) * bW
                                        ) / K;

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





class LineSolver
{
    using TC = TransportCoefficients::ENUMDATA; 

    public:

        LineSolver( ArrayAllocator<Fields, array3D> &fields,
                    FVCoefficients &fvCoeffs,
                    const InputData::LineSolverSettings &lineSolverSettings) :
        m_fields( fields ),
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
            static_assert( (Vstag == n) || (Vstag == s) || (Vstag == p), "Invalid V momentum staggering" );
            static_assert( (Wstag == t) || (Wstag == b) || (Wstag == p), "Invalid W momentum staggering" );

            // Temporary for storing new block update
            EnumVector<Fields, floatType> oldBlock;

            // Residuals 
            EnumVector<Fields, floatType> delta, residuals, residualsInitialInv;

            // Staggered indexing
            EnumVector<Fields, intType> iS( {0, 0, 0, 0} ), 
                                        jS( {j, j+CoeffIndex[Vstag], j, j} ), 
                                        kS( {k, k, k+CoeffIndex[Wstag], k} );


            // Lambda to update block in place, relax, and update residual
            auto UpdateAndRelax = [&]<TransportCoefficients::ENUMDATA Ustag>(intType i) {

                EnumFor<Fields>( [&] (Fields::ENUMDATA f) { iS[f] = i; } );      // Set iterating coefficient
                iS[U] += CoeffIndex[Ustag];                                      // U momentum is staggered    

                EnumFor<Fields>( [&] (Fields::ENUMDATA f) { oldBlock[f] = m_fields[f]( G(iS[f], jS[f], kS[f]) ); } );  // Set old block values

                m_blockSolver.UpdateBlock<Ustag, Vstag, Wstag>(i, j, k);

                EnumFor<Fields>( [&] (Fields::ENUMDATA f) 
                { 
                    auto &fieldBlock = m_fields[f]( G(iS[f], jS[f], kS[f]) );
                    delta[f] = m_relaxation[f] * ( fieldBlock - oldBlock[f] );   // Relaxed change in solution      
                    fieldBlock = oldBlock[f] + delta[f];                         // Apply relaxation
                    residuals[f] += abs( delta[f] );                             // Add to residual count
                } );
            };


            // Solver loop
            intType nIterations = 0;
            while ( nIterations < m_maxIterations ) 
            {
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] = 0.0f; } );

                // Update in place and relax
                for (intType i = 0; i != ni-1; i++) {   // Forward sweep
                    UpdateAndRelax.template operator()<e>(i);
                }

                for (intType i = ni-1; i != 0; i--) {   // Backward sweep
                    UpdateAndRelax.template operator()<w>(i);
                }


                // Normalise residuals
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] /= static_cast<floatType>( ni ); } );
                RelativeResidual( residuals, residualsInitialInv, nIterations );
                nIterations++;

                // Check residual
                if ( MetResidualTolerence( residuals, m_maxResiduals ) ) {
                    break;
                }

            }

        }

    private:

        ArrayAllocator<Fields, array3D> &m_fields;
        const intType m_maxIterations;
        const EnumVector<Fields, floatType> m_maxResiduals;
        const EnumVector<Fields, floatType> m_relaxation;
        BlockSolver m_blockSolver;
};





class PlaneSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    public:

        PlaneSolver( ArrayAllocator<Fields, array3D> &fields,
                     FVCoefficients &fvCoeffs,
                     const InputData::PlaneSolverSettings &planeSolverSettings,
                     const InputData::LineSolverSettings &lineSolverSettings) :
        m_fields( fields ),
        m_maxIterations( planeSolverSettings.maxIterations ),
        m_maxResiduals( planeSolverSettings.maxResiduals ),
        m_relaxation( planeSolverSettings.relaxation ),
        m_lineSolver( fields, fvCoeffs, lineSolverSettings )
        {}


        template<TC Wstag>
        void SolvePlane(const intType k)
        {
            using enum Fields::ENUMDATA;
            using enum Axis::ENUMDATA;
            using enum TransportCoefficients::ENUMDATA;
            intType ni = m_fields[U].dimension(X) - 2*nGhost;
            intType nj = m_fields[U].dimension(Y) - 2*nGhost;

            // Staggering must be valid
            static_assert( (Wstag == t) || (Wstag == b) || (Wstag == p), "Invalid W momentum staggering" );

            // Old plane solution
            EnumVector<Fields, array1D> oldLine( array1D( m_fields[U].dimension(X) ) );

            // Residuals
            EnumVector<Fields, floatType> residuals, residualsInitialInv;
            EnumVector<Fields, array1D> delta( array1D( m_fields[U].dimension(X) ) );

            // Staggered indexing
            EnumVector<Fields, intType> jS( {0, 0, 0, 0} ), 
                                        kS( {k, k, k+CoeffIndex[Wstag], k} );


            // Lambda to update line in place, relax, and update residual
            auto UpdateAndRelax = [&]<TransportCoefficients::ENUMDATA Vstag>(intType j) {

                EnumFor<Fields>( [&] (Fields::ENUMDATA f) { jS[f] = j; } );      // Set iterating coefficient
                jS[V] += CoeffIndex[Vstag];                                     // V momentum is staggered   

                EnumFor<Fields>( [&] (Fields::ENUMDATA f) { oldLine[f] = m_fields[f].chip(kS[f], Z).chip(jS[f], Y); } );      // Set old line

                m_lineSolver.SolveLine<Vstag, Wstag>(j, k);

                EnumFor<Fields> ( [&] (Fields::ENUMDATA f) {
                    auto fieldLine = m_fields[f].chip( G(kS[f]), Z).chip( G(jS[f]), Y);
                    delta[f] = delta[f].constant( m_relaxation[f] ) * ( fieldLine - oldLine[f] );   // Relaxed change in line
                    fieldLine = oldLine[f] + delta[f];                                              // Relax
                    residuals[f] += static_cast<array0D>( delta[f].abs().sum() )(0);                                                // Add to residual count
                } );
            };

            // Solver loop
            intType nIterations = 0;
            while ( nIterations < m_maxIterations ) 
            {
                // EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] = 0.0f; } );

                for (intType j = 0; j != nj-1; j++) {   // Forward sweep
                    UpdateAndRelax.template operator()<n>(j);
                }

                for (intType j = nj-1; j != 0; j--) {   // Backward sweep
                    UpdateAndRelax.template operator()<s>(j);
                }

                // Normalise residuals
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] /= static_cast<floatType>( ni*nj ); } );
                RelativeResidual( residuals, residualsInitialInv, nIterations );
                nIterations++;

                // Check residual tolerence
                if ( MetResidualTolerence( residuals, m_maxResiduals ) ) {
                    break;
                }
            }

        }


    private:
        ArrayAllocator<Fields, array3D> &m_fields;
        const intType m_maxIterations;
        const EnumVector<Fields, floatType> m_maxResiduals;
        const EnumVector<Fields, floatType> m_relaxation;
        LineSolver m_lineSolver;
};





void SweepSolve( ArrayAllocator<Fields, array3D> &fields, 
                 const Mesh &mesh, 
                 const InputData &inputData) 
{
    using enum Axis::ENUMDATA;
    using enum Fields::ENUMDATA;
    using TC = TransportCoefficients::ENUMDATA;

    // Extract from input data
    const InputData::PlaneSweepSettings  planeSweepSettings  = inputData.planeSweepSettings;
    const InputData::PlaneSolverSettings planeSolverSettings = inputData.planeSolverSettings;
    const InputData::LineSolverSettings  lineSolverSettings  = inputData.lineSolverSettings;

    // Initialise 
    ArrayAllocator<Fields, array3D> faceVelocities = InitialiseFaceVelocities( mesh, fields, inputData );
    ArrayAllocator<Fields, array3D> fieldsOld( fields );
    FVCoefficients fvCoeffs = InitialiseFVCoefficients( mesh, fields, faceVelocities, inputData );
    
    // Plane solver
    PlaneSolver planeSolver(fields, fvCoeffs, planeSolverSettings, lineSolverSettings);

    const EnumVector<Fields, floatType> &relaxation  = planeSweepSettings.relaxation;
    intType nOuterIterations, 
            nInnerIterations;
    EnumVector<Fields, floatType> residualsInner, 
                                  residualsOuter, 
                                  residualsOuterInitialInv, 
                                  residualsInnerInitialInv;
    std::vector< EnumVector<Fields, floatType> > residualsHistory;
    residualsHistory.reserve( static_cast<size_t>( planeSweepSettings.maxOuterIterations ) );
    
    EnumVector<Fields, array2D> delta( array2D( fields[U].dimension(X), fields[U].dimension(Y) ) ),
                                oldPlane( array2D( fields[U].dimension(X), fields[U].dimension(Y) ) );
    intType ni = mesh.nCells[X],
            nj = mesh.nCells[Y],
            nk = mesh.nCells[Z];
    EnumVector<Fields, intType> kS( {0, 0, 0, 0} );


    // Lambda to update plane in place, relax, and update residual
    auto UpdateAndRelax = [&]<TransportCoefficients::ENUMDATA Wstag>(intType k) {
        
        EnumFor<Fields>( [&] (Fields::ENUMDATA f) { kS[f] = k; } );       // Set iterating coefficient
        kS[W] += CoeffIndex[Wstag];                                       // W momentum is staggered

        EnumFor<Fields>( [&] (Fields::ENUMDATA f) { oldPlane[f] = fields[f].chip(kS[f], Z); } );    // Set old plane

        planeSolver.SolvePlane<Wstag>(k);   // Solve in place

        EnumFor<Fields> ( [&] (Fields::ENUMDATA f) {
            auto fieldPlane = fields[f].chip( G(kS[f]), Z);
            delta[f] = delta[f].constant( relaxation[f] ) * ( fieldPlane - oldPlane[f] );     // Relaxed change in plane
            fieldPlane = oldPlane[f] + delta[f];                                              // Relax
            residualsInner[f] += static_cast<array0D>( delta[f].abs().sum() )(0);             // Add to residual count
        } );
    };

    // Outer iterations
    nOuterIterations = 0;
    while ( nOuterIterations < planeSweepSettings.maxOuterIterations ) 
    {

        // Inner iterations
        nInnerIterations = 0;
        while ( nInnerIterations < planeSweepSettings.maxInnerIterations ) {
        
            EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residualsInner[field] = 0.0f; } );

            // Update plane
            for (intType k = 0; k != nk-1; k++) {   // Forward sweep
                UpdateAndRelax.template operator()<TC::t>(k);
            }

            for (intType k = nk-1; k != 0; k--) {   // Backward sweep
                UpdateAndRelax.template operator()<TC::b>(k);
            }

            UTIL::WriteArray("U_velocity.dbg", fields[U]);
            UTIL::WriteArray("V_velocity.dbg", fields[V]);
            UTIL::WriteArray("W_velocity.dbg", fields[W]);
            UTIL::WriteArray("Pressure.dbg"  , fields[P]);

            // Normalise residuals
            EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residualsInner[field] /= static_cast<floatType>( ni*nj*nk ); } );
            RelativeResidual( residualsInner, residualsInnerInitialInv, nInnerIterations );
            nInnerIterations++;

            std::cout << "Inner iteration: " << nInnerIterations << ", U residual: " << residualsInner[U] << "\n\n";

            // Check residual tolerence
            if ( MetResidualTolerence(residualsInner, planeSweepSettings.maxInnerResiduals) ) {
                std::cout << "*** INNER ITERATIONS CONVERGED ***" << "\n\n";
                break;
            }
        }
    
        

        // Update residuals
        L1ArrayDiff( residualsOuter, fields, fieldsOld );
        RelativeResidual( residualsOuter, residualsOuterInitialInv, nOuterIterations );
        residualsHistory.push_back( residualsOuter );
        nOuterIterations++;

        // Set old fields
        fieldsOld = fields;

        // Check residual tolerence
        // if ( MetResidualTolerence(residualsOuter, planeSweepSettings.maxOuterResiduals) ) {
        //     std::cout << "*** OUTER ITERATIONS CONVERGED ***" << "\n\n";
        //     break;
        // }

        // Update nonlinear coefficients
        UpdateFaceVelocities( faceVelocities, mesh, fields, inputData );
        UpdateFVCoefficients( fvCoeffs, mesh, fields, faceVelocities, inputData );

    }


    // Strip unused data
    residualsHistory.shrink_to_fit();
 
}


}   // end namespace CFD