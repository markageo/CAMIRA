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

#include <iostream>

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
    TriadSolver( FieldData<Tensor3D> &fields,
                 const FieldData<Tensor3D> &fieldsOld,
                 const Tensor3D &mask,
                 const FVCoefficients &fvCoeffs ) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_mask( mask ),
                    m_fvCoeffs( fvCoeffs ),
                    m_ni( fvCoeffs.nCells(0) ),
                    m_nj( fvCoeffs.nCells(1) ),
                    m_nk (fvCoeffs.nCells(2) ),
                    m_K( Tensor3D(m_ni, m_nj, m_nk).setZero() )
    { UpdateGlobalConstants(); };



    // Core function which updates the local coupled system. Templated by staggering direction.
    // This makes use of precomputed line constants
    __attribute__((always_inline)) 
    inline void UpdateTriad( const intType i, 
                             const intType j, 
                             const intType k,
                             const FieldData<Tensor1D> &lineConstants )
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
            newtonStencilX = - m_fvCoeffs.Mom[X].AU[Y][sCV::cCoupled]( iU, jU, kU ) * m_fields.U[Y]( igU, jgU+sCV::iCoupled, kgU )

                             - m_fvCoeffs.Mom[X].AU[Z][sCW::cCoupled]( iU, jU, kU ) * m_fields.U[Z]( igU, jgU, kgU+sCW::iCoupled );
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


        // Only update the molecule if none of the cells are within the immersed boundary
        floatType masterMask = m_mask(iU, jU, kU) * m_mask(iV, jV, kV) * m_mask(iW, jW, kW) * m_mask(i, j, k);

        // Update P from continuity
        floatType newP = ( 1 - m_fvCoeffs.Cont.relaxation ) * m_fields.P( ig, jg, kg )
                                 + m_fvCoeffs.Cont.relaxation * 
                                   ( bP 
                                   - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](i) * bU 
                                   - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](j) * bV 
                                   - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](k) * bW 
                                   ) * m_K(i, j, k);

        // Update U from momentum
        floatType newU = ( 1 - m_fvCoeffs.Mom[X].relaxation) * m_fields.U[X]( igU, jgU, kgU )
                                       + m_fvCoeffs.Mom[X].relaxation * ( bU - m_fvCoeffs.Mom[X].AP[sUP::cCoupled](iU) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[X].diagCoeffInv(iU, jU, kU) );

        // Update V from momentum
        floatType newV = ( 1 - m_fvCoeffs.Mom[Y].relaxation ) * m_fields.U[Y]( igV, jgV, kgV )
                                       + m_fvCoeffs.Mom[Y].relaxation * ( bV - m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jV) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Y].diagCoeffInv(iV, jV, kV) );

        // Update W from momentum
        floatType newW = ( 1 - m_fvCoeffs.Mom[Z].relaxation) * m_fields.U[Z]( igW, jgW, kgW ) 
                                       + m_fvCoeffs.Mom[Z].relaxation * ( bW - m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kW) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Z].diagCoeffInv(iW, jW, kW) );


        // Updating like this means that the old pressure value is used in the momentum update. 
        m_fields.P( ig, jg, kg )       = (1.0f - masterMask) * m_fields.P( ig, jg, kg )        +  masterMask * newP;
        m_fields.U[X]( igU, jgU, kgU ) = (1.0f - masterMask) * m_fields.U[X]( igU, jgU, kgU )  +  masterMask * newU;
        m_fields.U[Y]( igV, jgV, kgV ) = (1.0f - masterMask) * m_fields.U[Y]( igV, jgV, kgV )  +  masterMask * newV;
        m_fields.U[Z]( igW, jgW, kgW ) = (1.0f - masterMask) * m_fields.U[Z]( igW, jgW, kgW )  +  masterMask * newW;
    }



    // Core function which updates the local coupled system. Templated by staggering direction.
    // This function evaluates the full stencil
    __attribute__((always_inline)) 
    inline void UpdateTriad( const intType i, 
                             const intType j, 
                             const intType k)
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

        // Just for now
        bool xOutOfBounds = (iU < 0) || (iU > m_ni-1),
             yOutOfBounds = (jV < 0) || (jV > m_nj-1),
             zOutOfBounds = (kW < 0) || (kW > m_nk-1);


        // Precompute momentum RHS divided by AP coefficients
        // U momentum
        floatType bU = 0.0f;
        if ( !xOutOfBounds ) {
            floatType newtonStencilX = 0.0f;
            if constexpr ( LI == Linearisation::Newton ) {
                newtonStencilX = - m_fvCoeffs.Mom[X].AU[Y][n]( iU, jU, kU ) * m_fields.U[Y]( igU  , jgU+1, kgU  )
                                 - m_fvCoeffs.Mom[X].AU[Y][p]( iU, jU, kU ) * m_fields.U[Y]( igU  , jgU  , kgU  )
                                 - m_fvCoeffs.Mom[X].AU[Y][s]( iU, jU, kU ) * m_fields.U[Y]( igU  , jgU-1, kgU  )

                                 - m_fvCoeffs.Mom[X].AU[Z][t]( iU, jU, kU ) * m_fields.U[Z]( igU  , jgU  , kgU+1)
                                 - m_fvCoeffs.Mom[X].AU[Z][p]( iU, jU, kU ) * m_fields.U[Z]( igU  , jgU  , kgU  )
                                 - m_fvCoeffs.Mom[X].AU[Z][b]( iU, jU, kU ) * m_fields.U[Z]( igU  , jgU  , kgU-1);
            }
            bU = ( m_fvCoeffs.Mom[X].F(iU, jU, kU)
                   - m_fvCoeffs.Mom[X].B(iU, jU, kU)  

                        - m_fvCoeffs.Mom[X].AU[X][n](iU, jU, kU) * m_fields.U[X]( igU  , jgU+1, kgU  )
                        - m_fvCoeffs.Mom[X].AU[X][e](iU, jU, kU) * m_fields.U[X]( igU+1, jgU  , kgU  )
                        - m_fvCoeffs.Mom[X].AU[X][s](iU, jU, kU) * m_fields.U[X]( igU  , jgU-1, kgU  )
                        - m_fvCoeffs.Mom[X].AU[X][w](iU, jU, kU) * m_fields.U[X]( igU-1, jgU  , kgU  )
                        - m_fvCoeffs.Mom[X].AU[X][t](iU, jU, kU) * m_fields.U[X]( igU  , jgU  , kgU+1) 
                        - m_fvCoeffs.Mom[X].AU[X][b](iU, jU, kU) * m_fields.U[X]( igU  , jgU  , kgU-1)

                        - m_fvCoeffs.Mom[X].AP[sUP::cLeft ](iU) * m_fields.P( igU + sUP::iLeft , jgU, kgU)
                        - m_fvCoeffs.Mom[X].AP[sUP::cRight](iU) * m_fields.P( igU + sUP::iRight, jgU, kgU) 

                        + newtonStencilX

                        ) * m_fvCoeffs.Mom[X].diagCoeffInv(iU, jU, kU);
        }
        


        // V momentum
        floatType bV = 0.0f;
        if ( !yOutOfBounds ) {
            floatType newtonStencilY = 0.0f;
            if constexpr ( LI == Linearisation::Newton ) {
                newtonStencilY = - m_fvCoeffs.Mom[Y].AU[X][e]( iV, jV, kV ) * m_fields.U[X]( igV+1, jgV  , kgV  )
                                - m_fvCoeffs.Mom[Y].AU[X][p]( iV, jV, kV ) * m_fields.U[X]( igV  , jgV  , kgV  )
                                - m_fvCoeffs.Mom[Y].AU[X][w]( iV, jU, kU ) * m_fields.U[X]( igV-1, jgV  , kgV )

                                - m_fvCoeffs.Mom[Y].AU[Z][t]( iV, jV, kV ) * m_fields.U[Z]( igV  , jgV  , kgV+1)
                                - m_fvCoeffs.Mom[Y].AU[Z][p]( iV, jV, kV ) * m_fields.U[Z]( igV  , jgV  , kgV  )
                                - m_fvCoeffs.Mom[Y].AU[Z][b]( iV, jV, kV ) * m_fields.U[Z]( igV  , jgV  , kgV-1);
            }
            bV = ( m_fvCoeffs.Mom[Y].F(iV, jV, kV)
                 - m_fvCoeffs.Mom[Y].B(iV, jV, kV)

                       - m_fvCoeffs.Mom[Y].AU[Y][n](iV, jV, kV) * m_fields.U[Y]( igV  , jgV+1, kgV  ) 
                       - m_fvCoeffs.Mom[Y].AU[Y][e](iV, jV, kV) * m_fields.U[Y]( igV+1, jgV  , kgV  ) 
                       - m_fvCoeffs.Mom[Y].AU[Y][s](iV, jV, kV) * m_fields.U[Y]( igV  , jgV-1, kgV  ) 
                       - m_fvCoeffs.Mom[Y].AU[Y][w](iV, jV, kV) * m_fields.U[Y]( igV-1, jgV  , kgV  ) 
                       - m_fvCoeffs.Mom[Y].AU[Y][t](iV, jV, kV) * m_fields.U[Y]( igV  , jgV  , kgV+1) 
                       - m_fvCoeffs.Mom[Y].AU[Y][b](iV, jV, kV) * m_fields.U[Y]( igV  , jgV  , kgV-1)

                       - m_fvCoeffs.Mom[Y].AP[sVP::cLeft ](jV) * m_fields.P( igV, jgV + sVP::iLeft , kgV)
                       - m_fvCoeffs.Mom[Y].AP[sVP::cRight](jV) * m_fields.P( igV, jgV + sVP::iRight, kgV)

                       + newtonStencilY

                       ) * m_fvCoeffs.Mom[Y].diagCoeffInv(iV, jV, kV);
        }
        


        // W momentum
        floatType bW = 0.0f;
        if ( !zOutOfBounds ) {
            floatType newtonStencilZ = 0.0f;
            if constexpr ( LI == Linearisation::Newton ) {
                newtonStencilZ =  - m_fvCoeffs.Mom[Z].AU[X][e](iW, jW, kW) * m_fields.U[X]( igW+1, jgW  , kgW  )
                                - m_fvCoeffs.Mom[Z].AU[X][p](iW, jW, kW) * m_fields.U[X]( igW  , jgW  , kgW  )
                                - m_fvCoeffs.Mom[Z].AU[X][w](iW, jW, kW) * m_fields.U[X]( igW-1, jgW  , kgW  )

                                - m_fvCoeffs.Mom[Z].AU[Y][n](iW, jW, kW) * m_fields.U[Y]( igW  , jgW+1, kgW  )
                                - m_fvCoeffs.Mom[Z].AU[Y][p](iW, jW, kW) * m_fields.U[Y]( igW  , jgW  , kgW  )
                                - m_fvCoeffs.Mom[Z].AU[Y][s](iW, jW, kW) * m_fields.U[Y]( igW  , jgW-1, kgW  );
            }
            bW = ( m_fvCoeffs.Mom[Z].F(iW, jW, kW)
                  - m_fvCoeffs.Mom[Z].B(iW, jW, kW)
                            
                       - m_fvCoeffs.Mom[Z].AU[Z][n](iW, jW, kW) * m_fields.U[Z]( igW  , jgW+1, kgW  ) 
                       - m_fvCoeffs.Mom[Z].AU[Z][e](iW, jW, kW) * m_fields.U[Z]( igW+1, jgW  , kgW  ) 
                       - m_fvCoeffs.Mom[Z].AU[Z][s](iW, jW, kW) * m_fields.U[Z]( igW  , jgW-1, kgW  ) 
                       - m_fvCoeffs.Mom[Z].AU[Z][w](iW, jW, kW) * m_fields.U[Z]( igW-1, jgW  , kgW  ) 
                       - m_fvCoeffs.Mom[Z].AU[Z][t](iW, jW, kW) * m_fields.U[Z]( igW  , jgW  , kgW+1) 
                       - m_fvCoeffs.Mom[Z].AU[Z][b](iW, jW, kW) * m_fields.U[Z]( igW  , jgW  , kgW-1)

                       - m_fvCoeffs.Mom[Z].AP[sWP::cLeft ](kW) * m_fields.P( igW, jgW, kgW + sWP::iLeft ) 
                       - m_fvCoeffs.Mom[Z].AP[sWP::cRight](kW) * m_fields.P( igW, jgW, kgW + sWP::iRight)

                       + newtonStencilZ

                       ) * m_fvCoeffs.Mom[Z].diagCoeffInv(iW, jW, kW);
        }


        // Continuity for pressure
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - m_fvCoeffs.Cont.AP[nn](i, j, k) * m_fields.P( ig  , jg+2, kg  )
                                  - m_fvCoeffs.Cont.AP[ee](i, j, k) * m_fields.P( ig+2, jg  , kg  ) 
                                  - m_fvCoeffs.Cont.AP[ss](i, j, k) * m_fields.P( ig  , jg-2, kg  ) 
                                  - m_fvCoeffs.Cont.AP[ww](i, j, k) * m_fields.P( ig-2, jg  , kg  ) 
                                  - m_fvCoeffs.Cont.AP[tt](i, j, k) * m_fields.P( ig  , jg  , kg+2) 
                                  - m_fvCoeffs.Cont.AP[bb](i, j, k) * m_fields.P( ig  , jg  , kg-2);
        }
        floatType bP = m_fvCoeffs.Cont.F(i, j, k)
                     - m_fvCoeffs.Cont.B(i, j, k)

                     - m_fvCoeffs.Cont.AU[X][sCU::cLeft ](i) * m_fields.U[X]( ig + sCU::iLeft , jg, kg)
                     - m_fvCoeffs.Cont.AU[X][sCU::cRight](i) * m_fields.U[X]( ig + sCU::iRight, jg, kg)

                     - m_fvCoeffs.Cont.AU[Y][sCV::cLeft ](j) * m_fields.U[Y]( ig, jg + sCV::iLeft , kg)
                     - m_fvCoeffs.Cont.AU[Y][sCV::cRight](j) * m_fields.U[Y]( ig, jg + sCV::iRight, kg)

                     - m_fvCoeffs.Cont.AU[Z][sCW::cLeft ](k) * m_fields.U[Z]( ig, jg, kg + sCW::iLeft )
                     - m_fvCoeffs.Cont.AU[Z][sCW::cRight](k) * m_fields.U[Z]( ig, jg, kg + sCW::iRight)

                     - m_fvCoeffs.Cont.AP[n](i, j, k) * m_fields.P( ig  , jg+1, kg  ) 
                     - m_fvCoeffs.Cont.AP[e](i, j, k) * m_fields.P( ig+1, jg  , kg  ) 
                     - m_fvCoeffs.Cont.AP[s](i, j, k) * m_fields.P( ig  , jg-1, kg  ) 
                     - m_fvCoeffs.Cont.AP[w](i, j, k) * m_fields.P( ig-1, jg  , kg  ) 
                     - m_fvCoeffs.Cont.AP[t](i, j, k) * m_fields.P( ig  , jg  , kg+1) 
                     - m_fvCoeffs.Cont.AP[b](i, j, k) * m_fields.P( ig  , jg  , kg-1)

                     + pressureWideStencil;


        // Only update the molecule if none of the cells are within the immersed boundary
        floatType masterMask = 1.0f;
        // floatType masterMask = m_mask(iU, jU, kU) * m_mask(iV, jV, kV) * m_mask(iW, jW, kW) * m_mask(i, j, k);

        // This should be precomputed
        floatType term1 = ( xOutOfBounds ) ? 0.0f : - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](i) * m_fvCoeffs.Mom[X].AP[sUP::cCoupled](iU) * m_fvCoeffs.Mom[X].diagCoeffInv(iU, jU, kU),
                  term2 = ( yOutOfBounds ) ? 0.0f : - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](j) * m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jV) * m_fvCoeffs.Mom[Y].diagCoeffInv(iV, jV, kV),
                  term3 = ( zOutOfBounds ) ? 0.0f : - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](k) * m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kW) * m_fvCoeffs.Mom[Z].diagCoeffInv(iW, jW, kW);
        floatType K = m_fvCoeffs.Cont.AP[p](i, j, k) 
                    + term1 + term2 + term3; 
        K = 1.0f / K;

        // Update P from continuity
        floatType newP = ( 1 - m_fvCoeffs.Cont.relaxation ) * m_fields.P( ig, jg, kg )
                                 + m_fvCoeffs.Cont.relaxation * 
                                   ( bP 
                                   - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](i) * bU 
                                   - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](j) * bV 
                                   - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](k) * bW 
                                   ) * K;

        // if ( xOutOfBounds || yOutOfBounds || zOutOfBounds ) {
        //     std::cout << bU << ", " << bV << ", " << bW << ", " << K << ":   " << newP << "\n";
        // }

        m_fields.P( ig, jg, kg )       = (1.0f - masterMask) * m_fields.P( ig, jg, kg )        +  masterMask * newP;

        // Update U from momentum
        floatType newU = m_fields.U[X]( igU, jgU, kgU );
        if ( !xOutOfBounds ) {
            newU = ( 1 - m_fvCoeffs.Mom[X].relaxation) * m_fields.U[X]( igU, jgU, kgU )
                                    + m_fvCoeffs.Mom[X].relaxation * ( bU - m_fvCoeffs.Mom[X].AP[sUP::cCoupled](iU) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[X].diagCoeffInv(iU, jU, kU) );
        } 
        
        // Update V from momentum
        floatType newV = m_fields.U[Y]( igV, jgV, kgV );
        if ( !yOutOfBounds ) {
            newV = ( 1 - m_fvCoeffs.Mom[Y].relaxation ) * m_fields.U[Y]( igV, jgV, kgV )
                                    + m_fvCoeffs.Mom[Y].relaxation * ( bV - m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jV) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Y].diagCoeffInv(iV, jV, kV) );
        }
        

        // Update W from momentum
        floatType newW = m_fields.U[Z]( igW, jgW, kgW );
        if ( !zOutOfBounds ) {
            newW = ( 1 - m_fvCoeffs.Mom[Z].relaxation) * m_fields.U[Z]( igW, jgW, kgW ) 
                                    + m_fvCoeffs.Mom[Z].relaxation * ( bW - m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kW) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Z].diagCoeffInv(iW, jW, kW) );
        }


        // Updating like this means that the old pressure value is used in the momentum update. 
        m_fields.U[X]( igU, jgU, kgU ) = (1.0f - masterMask) * m_fields.U[X]( igU, jgU, kgU )  +  masterMask * newU;
        m_fields.U[Y]( igV, jgV, kgV ) = (1.0f - masterMask) * m_fields.U[Y]( igV, jgV, kgV )  +  masterMask * newV;
        m_fields.U[Z]( igW, jgW, kgW ) = (1.0f - masterMask) * m_fields.U[Z]( igW, jgW, kgW )  +  masterMask * newW;

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
    FieldData<Tensor3D> &m_fields;
    const FieldData<Tensor3D> &m_fieldsOld;
    const Tensor3D &m_mask;
    const FVCoefficients &m_fvCoeffs;
    const intType m_ni, m_nj, m_nk;
    Tensor3D m_K;

};

}   // end namespace CFD    


#endif // TRIAD_SOLVER