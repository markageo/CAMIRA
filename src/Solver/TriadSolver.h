#ifndef TRIAD_SOLVER
#define TRIAD_SOLVER

#include "StaggerIndexing.h"

#include "../Types.h"
#include "../Macros.h"
#include "../IO/InputProcessing.h"
#include "../Tools/SweepTransformations.h"
#include "../Tools/FVTools.h"
#include "../Tools/FVLookups.h"
#include "../FiniteVolume/FiniteVolume.h"

namespace CFD
{

template < TransportCoefficients::ENUMDATA Ustag,
           TransportCoefficients::ENUMDATA Vstag,
           TransportCoefficients::ENUMDATA Wstag, 
           MomentumInterpolation MI, 
           Linearisation LI>
class TriadSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    // Staggering must be valid
    static_assert( (Ustag == TC::e) || (Ustag == TC::w) || (Ustag == TC::p), "Invalid U momentum staggering" );
    static_assert( (Vstag == TC::n) || (Vstag == TC::s) || (Vstag == TC::p), "Invalid V momentum staggering" );
    static_assert( (Wstag == TC::t) || (Wstag == TC::b) || (Wstag == TC::p), "Invalid W momentum staggering" );

    // Aliases for the staggering offsets
    using sCU = typename StaggerIndexing< Axis::X, Ustag >::ContinuityVelocity;
    using sUP = typename StaggerIndexing< Axis::X, Ustag >::MomentumPressure;

    using sCV = typename StaggerIndexing< Axis::Y, Vstag >::ContinuityVelocity;
    using sVP = typename StaggerIndexing< Axis::Y, Vstag >::MomentumPressure;

    using sCW = typename StaggerIndexing< Axis::Z, Wstag >::ContinuityVelocity;
    using sWP = typename StaggerIndexing< Axis::Z, Wstag >::MomentumPressure;

public:
    TriadSolver( FieldData<array3D> &fields,
                 const FieldData<array3D> &fieldsOld,
                 const FVCoefficients &fvCoeffs ) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_fvCoeffs( fvCoeffs ),
                    m_ni( fvCoeffs.nCells(0) ),
                    m_nj( fvCoeffs.nCells(1) ),
                    m_nk (fvCoeffs.nCells(2) ),
                    m_K( array3D(m_ni, m_nj, m_nk).setZero() )
    { UpdateGlobalConstants(); };



    // Core function which updates the local coupled system. Templated by staggering direction.
    // This makes use of precomputed line constants
    __attribute__((always_inline)) 
    inline void UpdateTriad( const intType i, 
                             const intType j, 
                             const intType k,
                             const FieldData<array1D> &lineConstants )
    {
        using namespace FVT;
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        // For indexing the staggered cells
        intType iU{ i + sCU::iCoupled }, jU{ j                 }, kU{ k                 }; // U momentum
        intType iV{ i                 }, jV{ j + sCV::iCoupled }, kV{ k                 }; // V momentum
        intType iW{ i                 }, jW{ j                 }, kW{ k + sCW::iCoupled }; // W momentum

        // With ghost cells, this is faster than using the G() function inline every time
        intType igU{ G(iU) }, jgU{ G(jU) }, kgU{ G(kU) };
        intType igV{ G(iV) }, jgV{ G(jV) }, kgV{ G(kV) };
        intType igW{ G(iW) }, jgW{ G(jW) }, kgW{ G(kW) };
        intType   ig{ G(i) },   jg{ G(j) },   kg{ G(k) };


        // Precompute momentum RHS divided by AP coefficients
        // U momentum
        floatType newtonStencilX = 0.0f;
        if constexpr ( LI == Linearisation::Newton ) {
            newtonStencilX = - m_fvCoeffs.Mom[X].AU[Y][sCV::cCoupled]( i, j, k ) * m_fields.U[Y]( ig, jg+sCV::iCoupled, kg )

                             - m_fvCoeffs.Mom[X].AU[Z][sCW::cCoupled]( i, j, k ) * m_fields.U[Z]( ig, jg, kg+sCW::iRight );
        }
        floatType bU = ( lineConstants.U[X](iU)  

                       - m_fvCoeffs.Mom[X].AU[X][e](iU, jU, kU) * m_fields.U[X]( igU+1, jgU  , kgU  )
                       - m_fvCoeffs.Mom[X].AU[X][w](iU, jU, kU) * m_fields.U[X]( igU-1, jgU  , kgU  )

                       - m_fvCoeffs.Mom[X].AP[sUP::cLeft ](iU) * m_fields.P( igU + sUP::iLeft , jgU, kgU)
                       - m_fvCoeffs.Mom[X].AP[sUP::cRight](iU) * m_fields.P( igU + sUP::iRight, jgU, kgU) 

                       + newtonStencilX

                       ) * m_fvCoeffs.Mom[X].diagCoeffInv(iU, jU, kU);


        // V momentum
        floatType newtonStencilY = 0.0f;
        floatType bV = ( lineConstants.U[Y](iV)

                       - m_fvCoeffs.Mom[Y].AU[Y][e](iV, jV, kV) * m_fields.U[Y]( igV+1, jgV  , kgV  ) 
                       - m_fvCoeffs.Mom[Y].AU[Y][w](iV, jV, kV) * m_fields.U[Y]( igV-1, jgV  , kgV  ) 

                       + newtonStencilY

                       ) * m_fvCoeffs.Mom[Y].diagCoeffInv(iV, jV, kV);


        // W momentum
        floatType newtonStencilZ = 0.0f;
        floatType bW = ( lineConstants.U[Z](iW)

                       - m_fvCoeffs.Mom[Z].AU[Z][e](iW, jW, kW) * m_fields.U[Z]( igW+1, jgW  , kgW  ) 
                       - m_fvCoeffs.Mom[Z].AU[Z][w](iW, jW, kW) * m_fields.U[Z]( igW-1, jgW  , kgW  ) 

                       + newtonStencilZ

                       ) * m_fvCoeffs.Mom[Z].diagCoeffInv(iW, jW, kW);


        // Continuity for pressure
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - m_fvCoeffs.Cont.AP[ee](i, j, k) * m_fields.P( ig+2, jg  , kg  ) 
                                  - m_fvCoeffs.Cont.AP[ww](i, j, k) * m_fields.P( ig-2, jg  , kg  );
        }
        floatType bP = lineConstants.P(i)

                     - m_fvCoeffs.Cont.AU[X][sCU::cLeft ](i) * m_fields.U[X]( ig + sCU::iLeft , jg, kg)
                     - m_fvCoeffs.Cont.AU[X][sCU::cRight](i) * m_fields.U[X]( ig + sCU::iRight, jg, kg)

                     - m_fvCoeffs.Cont.AP[e](i, j, k) * m_fields.P( ig+1, jg  , kg  ) 
                     - m_fvCoeffs.Cont.AP[w](i, j, k) * m_fields.P( ig-1, jg  , kg  ) 

                     + pressureWideStencil;


        // Update P from continuity
        m_fields.P( ig, jg, kg ) = ( 1 - m_fvCoeffs.Cont.relaxation ) * m_fieldsOld.P( ig, jg, kg )
                                 + m_fvCoeffs.Cont.relaxation * 
                                   ( bP 
                                   - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](i) * bU 
                                   - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](j) * bV 
                                   - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](k) * bW 
                                   ) * m_K(i, j, k);


        // Update U from momentum
        m_fields.U[X]( igU, jgU, kgU ) = ( 1 - m_fvCoeffs.Mom[X].relaxation ) * m_fieldsOld.U[X]( igU, jgU, kgU )
                                       + m_fvCoeffs.Mom[X].relaxation * ( bU - m_fvCoeffs.Mom[X].AP[sUP::cCoupled](iU) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[X].diagCoeffInv(iU, jU, kU) );

        // Update V from momentum
        m_fields.U[Y]( igV, jgV, kgV ) = ( 1 - m_fvCoeffs.Mom[Y].relaxation ) * m_fieldsOld.U[Y]( igV, jgV, kgV )
                                       + m_fvCoeffs.Mom[Y].relaxation * ( bV - m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jV) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Y].diagCoeffInv(iV, jV, kV) );

        // Update W from momentum
        m_fields.U[Z]( igW, jgW, kgW ) = ( 1 - m_fvCoeffs.Mom[Z].relaxation ) * m_fieldsOld.U[Z]( igW, jgW, kgW ) 
                                       + m_fvCoeffs.Mom[Z].relaxation * ( bW - m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kW) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Z].diagCoeffInv(iW, jW, kW) );

    }


    // Constants which are global to the linear solver
    void UpdateGlobalConstants()
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Axis::ENUMDATA;

        // Staggered indexing for fields
        intType iU, jU, kU,
                iV, jV, kV,
                iW, jW, kW;

        // Starting and ending indices, since K cannot be calculated on some boundaries due to the staggering
        intType iStart = 1 + sCU::iLeft,
                iLength = m_ni - 1 + sCU::iRight,

                jStart = 1 + sCV::iLeft,
                jLength = m_nj - 1 + sCV::iRight,

                kStart = 1 + sCW::iLeft,
                kLength = m_nk - 1 + sCW::iRight;

        for (intType k = kStart; k != kLength; k++) {

            kU = k;
            kV = k;
            kW = k + sCW::iCoupled;

            for (intType j = jStart; j != jLength; j++) {

                jU = j;
                jV = j + sCV::iCoupled;
                jW = j;

                CFD_PRAGMA_VECTORIZE
                for (intType i = iStart; i != iLength; i++) {

                    iU = i + sCU::iCoupled;
                    iV = i;
                    iW = i;

                    m_K(i, j, k) = m_fvCoeffs.Cont.AP[p](i, j, k) 
                                 - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](i) * m_fvCoeffs.Mom[X].AP[sUP::cCoupled](iU) * m_fvCoeffs.Mom[X].diagCoeffInv(iU, jU, kU) 
                                 - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](j) * m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jV) * m_fvCoeffs.Mom[Y].diagCoeffInv(iV, jV, kV) 
                                 - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](k) * m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kW) * m_fvCoeffs.Mom[Z].diagCoeffInv(iW, jW, kW);
                    m_K(i, j, k) = 1.0f / m_K(i, j, k);
                }
            }
        }
    }


private:
    FieldData<array3D> &m_fields;
    const FieldData<array3D> &m_fieldsOld;
    const FVCoefficients &m_fvCoeffs;
    const intType m_ni, m_nj, m_nk;
    array3D m_K;

};

}   // end namespace CFD    


#endif // TRIAD_SOLVER