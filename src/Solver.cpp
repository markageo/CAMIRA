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
        EnumFor<Fields>( [&] (Fields::ENUMDATA field) { 
            if ( residuals[field] != 0 ) {  // Division by zero
                residualsInitialInv[field] = 1.0f / residuals[field]; 
            } else {
                residualsInitialInv[field] = 1.0f;
            }
            
        } );
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
template< TransportCoefficients::ENUMDATA Ustag,
          TransportCoefficients::ENUMDATA Vstag,
          TransportCoefficients::ENUMDATA Wstag >
class BlockSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    // Staggering must be valid
    static_assert( (Ustag == TC::e) || (Ustag == TC::w) || (Ustag == TC::p), "Invalid U momentum staggering" );
    static_assert( (Vstag == TC::n) || (Vstag == TC::s) || (Vstag == TC::p), "Invalid V momentum staggering" );
    static_assert( (Wstag == TC::t) || (Wstag == TC::b) || (Wstag == TC::p), "Invalid W momentum staggering" );

    public:

        BlockSolver(ArrayAllocator<Fields, array3D> &fields, const FVCoefficients &fvCoeffs) :
            m_fields( fields ), 
            m_fvCoeffs( fvCoeffs ) 
        {};


        // Core function which updates the local coupled system. Templated by staggering direction.
        void UpdateBlock(const intType i, const intType j, const intType k)
        {
            using enum Fields::ENUMDATA;
            using enum TransportCoefficients::ENUMDATA;

            // Indexing variables to take care of staggering - TODO: these should be made member variables
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
            K  = 1.0f / K;


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




template< TransportCoefficients::ENUMDATA Vstag,
          TransportCoefficients::ENUMDATA Wstag >
class LineSolver
{
    using TC = TransportCoefficients::ENUMDATA; 

    // Staggering must be valid
    static_assert( (Vstag == TC::n) || (Vstag == TC::s) || (Vstag == TC::p), "Invalid V momentum staggering" );
    static_assert( (Wstag == TC::t) || (Wstag == TC::b) || (Wstag == TC::p), "Invalid W momentum staggering" );

    public:

        LineSolver( ArrayAllocator<Fields, array3D> &fields,
                    FVCoefficients &fvCoeffs,
                    const InputData::LineSolverSettings &lineSolverSettings) :
        m_fields( fields ),
        m_maxIterations( lineSolverSettings.maxIterations ),
        m_maxResiduals( lineSolverSettings.maxResiduals ),
        m_relaxation( lineSolverSettings.relaxation ),
        m_blockSolverEast( fields, fvCoeffs ),
        m_blockSolverWest( fields, fvCoeffs)
        {}


        void SolveLine(const intType j, const intType k)
        {
            using enum TransportCoefficients::ENUMDATA;
            using enum Fields::ENUMDATA;
            using enum Axis::ENUMDATA;
            intType ni = m_fields[U].dimension(X) - 2*nGhost;

            // Temporary for storing new block update
            EnumVector<Fields, floatType> oldBlock;

            // Residuals 
            EnumVector<Fields, floatType> delta, residuals, residualsInitialInv;

            // Staggered indexing
            EnumVector<Fields, intType> iS( {0, 0, 0, 0} ), 
                                        jS( {j, j+CoeffIndex[Vstag], j, j} ), 
                                        kS( {k, k, k+CoeffIndex[Wstag], k} );


            // Lambda to update block in place, relax, and update residual
            auto UpdateAndRelax = [&]<TC Ustag>(BlockSolver< Ustag, Vstag, Wstag > &blockSolver, intType i ) {

                EnumFor<Fields>( [&] (Fields::ENUMDATA f) { iS[f] = i; } );      // Set iterating coefficient
                iS[U] += CoeffIndex[Ustag];                                      // U momentum is staggered    

                EnumFor<Fields>( [&] (Fields::ENUMDATA f) { oldBlock[f] = m_fields[f]( G(iS[f], jS[f], kS[f]) ); } );  // Set old block values

                blockSolver.UpdateBlock(i, j, k);

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
                    UpdateAndRelax( m_blockSolverEast, i );
                }

                for (intType i = ni-1; i != 0; i--) {   // Backward sweep
                    UpdateAndRelax( m_blockSolverWest, i );
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

        BlockSolver< TC::e, Vstag, Wstag > m_blockSolverEast;
        BlockSolver< TC::w, Vstag, Wstag > m_blockSolverWest; 
};





template< TransportCoefficients::ENUMDATA Wstag >
class PlaneSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    // Staggering must be valid
    static_assert( (Wstag == TC::t) || (Wstag == TC::b) || (Wstag == TC::p), "Invalid W momentum staggering" );

    public:

        PlaneSolver( ArrayAllocator<Fields, array3D> &fields,
                     FVCoefficients &fvCoeffs,
                     const InputData::PlaneSolverSettings &planeSolverSettings,
                     const InputData::LineSolverSettings &lineSolverSettings) :
        m_fields( fields ),
        m_maxIterations( planeSolverSettings.maxIterations ),
        m_maxResiduals( planeSolverSettings.maxResiduals ),
        m_relaxation( planeSolverSettings.relaxation ),
        m_lineSolverNorth( fields, fvCoeffs, lineSolverSettings ),
        m_lineSolverSouth( fields, fvCoeffs, lineSolverSettings )
        {}


        void SolvePlane(const intType k)
        {
            using enum Fields::ENUMDATA;
            using enum Axis::ENUMDATA;
            using enum TransportCoefficients::ENUMDATA;
            intType ni = m_fields[U].dimension(X) - 2*nGhost;
            intType nj = m_fields[U].dimension(Y) - 2*nGhost;

            // Old plane solution
            EnumVector<Fields, array1D> oldLine( array1D( m_fields[U].dimension(X) ) );

            // Residuals
            EnumVector<Fields, floatType> residuals, residualsInitialInv;
            EnumVector<Fields, array1D> delta( array1D( m_fields[U].dimension(X) ) );

            // Staggered indexing
            EnumVector<Fields, intType> jS( {0, 0, 0, 0} ), 
                                        kS( {k, k, k+CoeffIndex[Wstag], k} );


            // Lambda to update line in place, relax, and update residual
            auto UpdateAndRelax = [&]<TC Vstag>( LineSolver< Vstag, Wstag > &lineSolver, intType j) {

                EnumFor<Fields>( [&] (Fields::ENUMDATA f) { jS[f] = j; } );      // Set iterating coefficient
                jS[V] += CoeffIndex[Vstag];                                     // V momentum is staggered   

                EnumFor<Fields>( [&] (Fields::ENUMDATA f) { oldLine[f] = m_fields[f].chip( G(kS[f]), Z).chip( G(jS[f]), Y); } );      // Set old line

                lineSolver.SolveLine(j, k);

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
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] = 0.0f; } );

                for (intType j = 0; j != nj-1; j++) {   // Forward sweep
                    UpdateAndRelax( m_lineSolverNorth, j );
                }

                for (intType j = nj-1; j != 0; j--) {   // Backward sweep
                    UpdateAndRelax( m_lineSolverSouth, j );
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

        LineSolver< TC::n, Wstag > m_lineSolverNorth;
        LineSolver< TC::s, Wstag > m_lineSolverSouth;
};





class LinearSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    public:

        LinearSolver( ArrayAllocator<Fields, array3D> &fields,
                      FVCoefficients &fvCoeffs,
                      const InputData::LinearSolverSettings &linearSolverSettings,
                      const InputData::PlaneSolverSettings &planeSolverSettings,
                      const InputData::LineSolverSettings &lineSolverSettings ) :
        m_fields( fields ),
        m_maxIterations( linearSolverSettings.maxIterations ),
        m_maxResiduals( linearSolverSettings.maxResiduals ),
        m_relaxation( linearSolverSettings.relaxation ),
        m_planeSolverTop( fields, fvCoeffs, planeSolverSettings, lineSolverSettings ),
        m_planeSolverBottom( fields, fvCoeffs, planeSolverSettings, lineSolverSettings )
        {}


        void Solve()
        {
            using enum Fields::ENUMDATA;
            using enum Axis::ENUMDATA;
            using enum TransportCoefficients::ENUMDATA;
            intType ni = m_fields[U].dimension(X) - 2*nGhost;
            intType nj = m_fields[U].dimension(Y) - 2*nGhost;
            intType nk = m_fields[U].dimension(Z) - 2*nGhost;

            EnumVector<Fields, array2D> delta( array2D( m_fields[U].dimension(X), m_fields[U].dimension(Y) ) ),
                                        oldPlane( array2D( m_fields[U].dimension(X), m_fields[U].dimension(Y) ) );

            EnumVector<Fields, floatType> residuals, residualsInitialInv;

            EnumVector<Fields, intType> kS( {0, 0, 0, 0} );


            // Lambda to update plane in place, relax, and update residual
            auto UpdateAndRelax = [&]<TC Wstag>( PlaneSolver< Wstag > &planeSolver, intType k) {
                
                EnumFor<Fields>( [&] (Fields::ENUMDATA f) { kS[f] = k; } );       // Set iterating coefficient
                kS[W] += CoeffIndex[Wstag];                                       // W momentum is staggered

                EnumFor<Fields>( [&] (Fields::ENUMDATA f) { oldPlane[f] = m_fields[f].chip( G(kS[f]), Z); } );    // Set old plane

                planeSolver.SolvePlane(k);

                EnumFor<Fields> ( [&] (Fields::ENUMDATA f) {
                    auto fieldPlane = m_fields[f].chip( G(kS[f]), Z);
                    delta[f] = delta[f].constant( m_relaxation[f] ) * ( fieldPlane - oldPlane[f] );     // Relaxed change in plane
                    fieldPlane = oldPlane[f] + delta[f];                                              // Relax
                    residuals[f] += static_cast<array0D>( delta[f].abs().sum() )(0);                  // Add to residual count
                } );
            };


            intType nIterations = 0;
            while ( nIterations < m_maxIterations ) {
            
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] = 0.0f; } );

                // Update plane
                for (intType k = 0; k != nk-1; k++) {   // Forward sweep
                    UpdateAndRelax( m_planeSolverTop, k );
                }

                for (intType k = nk-1; k != 0; k--) {   // Backward sweep
                    UpdateAndRelax( m_planeSolverBottom, k );
                }

                // Normalise residuals
                EnumFor<Fields>( [&] (Fields::ENUMDATA field) { residuals[field] /= static_cast<floatType>( ni*nj*nk ); } );
                RelativeResidual( residuals, residualsInitialInv, nIterations );
                nIterations++;

                // Check residual tolerence
                if ( MetResidualTolerence( residuals, m_maxResiduals ) ) {
                    std::cout << "*** INNER ITERATIONS CONVERGED ***" << "\n\n";
                    break;
                }
            }

        }


    private:
        ArrayAllocator<Fields, array3D> &m_fields;
        const intType m_maxIterations;
        const EnumVector<Fields, floatType> m_maxResiduals;
        const EnumVector<Fields, floatType> m_relaxation;

        PlaneSolver< TC::t > m_planeSolverTop;
        PlaneSolver< TC::b > m_planeSolverBottom;

};




void SweepSolve( ArrayAllocator<Fields, array3D> &fields, 
                 const Mesh &mesh, 
                 const InputData &inputData) 
{
    using enum Axis::ENUMDATA;
    using enum Fields::ENUMDATA;

    // Extract from input data
    const InputData::LinearSolverSettings linearSolverSettings = inputData.linearSolverSettings;
    const InputData::PlaneSolverSettings  planeSolverSettings  = inputData.planeSolverSettings;
    const InputData::LineSolverSettings   lineSolverSettings   = inputData.lineSolverSettings;

    const intType maxOuterIterations = inputData.schemes.maxOuterIterations;
    const EnumVector<Fields, floatType> maxOuterResiduals = inputData.schemes.maxOuterResiduals;

    // Initialise 
    ArrayAllocator<Fields, array3D> faceVelocities = InitialiseFaceVelocities( mesh, fields, inputData );
    ArrayAllocator<Fields, array3D> fieldsOld( fields );
    FVCoefficients fvCoeffs = InitialiseFVCoefficients( mesh, fields, faceVelocities, inputData );
    
    intType nOuterIterations;
    EnumVector<Fields, floatType> residualsOuter, residualsOuterInitialInv;
    std::vector< EnumVector<Fields, floatType> > residualsHistory;
    residualsHistory.reserve( static_cast<size_t>( inputData.schemes.maxOuterIterations ) );
    
    // Instantiate linear solver
    LinearSolver linearSolver( fields, fvCoeffs, linearSolverSettings, planeSolverSettings, lineSolverSettings );


    // Outer iterations
    nOuterIterations = 0;
    while ( nOuterIterations < maxOuterIterations ) 
    {

        // Solve linearised equations    
        linearSolver.Solve();

        // Update residuals
        L1ArrayDiff( residualsOuter, fields, fieldsOld );
        RelativeResidual( residualsOuter, residualsOuterInitialInv, nOuterIterations );
        residualsHistory.push_back( residualsOuter );
        nOuterIterations++;

        // Set old fields
        fieldsOld = fields;

        std::cout << "Outer iteration: " << nOuterIterations << ", U outer residual: " << residualsOuter[U]
                                                             << ", V outer residual: " << residualsOuter[V]
                                                             << ", W outer residual: " << residualsOuter[W]
                                                             << "\n\n";

        // Check residual tolerence
        // if ( MetResidualTolerence(residualsOuter, maxOuterResiduals) ) {
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