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
                 const Tensor3D &mask,
                 const FVCoefficients &fvCoeffs ) : 
                    m_fields( fields ),
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

        // For indexing the staggered cells (with ghost)
        const intType igU{ G( i + sCU::iCoupled ) }, jgU{ G( j                 ) }, kgU{ G( k                 ) },
                      igV{ G( i                 ) }, jgV{ G( j + sCV::iCoupled ) }, kgV{ G( k                 ) },
                      igW{ G( i                 ) }, jgW{ G( j                 ) }, kgW{ G( k + sCW::iCoupled ) },
                       ig{ G( i                 ) },  jg{ G( j                 ) },  kg{ G( k                 ) };


        // Precompute momentum RHS divided by AP coefficients
        // U momentum
        floatType newtonStencilX = 0.0f;
        if constexpr ( LI == Linearisation::Newton ) {
            newtonStencilX = - m_fvCoeffs.Mom[X].AU[Y][sCV::cCoupled]( igU, jgU, kgU ) * m_fields.U[Y]( igU, jgU+sCV::iCoupled, kgU )

                             - m_fvCoeffs.Mom[X].AU[Z][sCW::cCoupled]( igU, jgU, kgU ) * m_fields.U[Z]( igU, jgU, kgU+sCW::iCoupled );
        }
        const floatType bU = ( lineConstants.U[X](igU)  

                             - m_fvCoeffs.Mom[X].AU[X][e](igU, jgU, kgU) * m_fields.U[X]( igU+1, jgU  , kgU  )
                             - m_fvCoeffs.Mom[X].AU[X][w](igU, jgU, kgU) * m_fields.U[X]( igU-1, jgU  , kgU  )

                             - m_fvCoeffs.Mom[X].AP[sUP::cLeft ](igU) * m_fields.P( igU + sUP::iLeft , jgU, kgU)
                             - m_fvCoeffs.Mom[X].AP[sUP::cRight](igU) * m_fields.P( igU + sUP::iRight, jgU, kgU) 

                             + newtonStencilX

                             ) * m_fvCoeffs.Mom[X].diagCoeffInv(igU, jgU, kgU);


        // V momentum
        floatType newtonStencilY = 0.0f;
        const floatType bV = ( lineConstants.U[Y](igV)

                             - m_fvCoeffs.Mom[Y].AU[Y][e](igV, jgV, kgV) * m_fields.U[Y]( igV+1, jgV  , kgV  ) 
                             - m_fvCoeffs.Mom[Y].AU[Y][w](igV, jgV, kgV) * m_fields.U[Y]( igV-1, jgV  , kgV  ) 

                             + newtonStencilY

                             ) * m_fvCoeffs.Mom[Y].diagCoeffInv(igV, jgV, kgV);


        // W momentum
        floatType newtonStencilZ = 0.0f;
        const floatType bW = ( lineConstants.U[Z](igW)

                             - m_fvCoeffs.Mom[Z].AU[Z][e](igW, jgW, kgW) * m_fields.U[Z]( igW+1, jgW  , kgW  ) 
                             - m_fvCoeffs.Mom[Z].AU[Z][w](igW, jgW, kgW) * m_fields.U[Z]( igW-1, jgW  , kgW  ) 

                             + newtonStencilZ

                             ) * m_fvCoeffs.Mom[Z].diagCoeffInv(igW, jgW, kgW);


        // Continuity for pressure
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - m_fvCoeffs.Cont.AP[ee](ig, jg, kg) * m_fields.P( ig+2, jg  , kg  ) 
                                  - m_fvCoeffs.Cont.AP[ww](ig, jg, kg) * m_fields.P( ig-2, jg  , kg  );
        }
        const floatType bP = lineConstants.P(ig)

                             - m_fvCoeffs.Cont.AU[X][sCU::cLeft ](ig) * m_fields.U[X]( ig + sCU::iLeft , jg, kg)
                             - m_fvCoeffs.Cont.AU[X][sCU::cRight](ig) * m_fields.U[X]( ig + sCU::iRight, jg, kg)

                             - m_fvCoeffs.Cont.AP[e](ig, jg, kg) * m_fields.P( ig+1, jg  , kg  ) 
                             - m_fvCoeffs.Cont.AP[w](ig, jg, kg) * m_fields.P( ig-1, jg  , kg  ) 

                             + pressureWideStencil;

        UpdateMolecule( i, j, k, bU, bV, bW, bP );
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

        // For indexing the staggered cells (with ghost)
        const intType igU{ G( i + sCU::iCoupled ) }, jgU{ G( j                 ) }, kgU{ G( k                 ) },
                      igV{ G( i                 ) }, jgV{ G( j + sCV::iCoupled ) }, kgV{ G( k                 ) },
                      igW{ G( i                 ) }, jgW{ G( j                 ) }, kgW{ G( k + sCW::iCoupled ) },
                       ig{ G( i                 ) },  jg{ G( j                 ) },  kg{ G( k                 ) };

        // Precompute momentum RHS divided by AP coefficients
        // U momentum
        floatType newtonStencilX = 0.0f;
        if constexpr ( LI == Linearisation::Newton ) {
            newtonStencilX = - m_fvCoeffs.Mom[X].AU[Y][n]( igU, jgU, kgU ) * m_fields.U[Y]( igU  , jgU+1, kgU  )
                             - m_fvCoeffs.Mom[X].AU[Y][p]( igU, jgU, kgU ) * m_fields.U[Y]( igU  , jgU  , kgU  )
                             - m_fvCoeffs.Mom[X].AU[Y][s]( igU, jgU, kgU ) * m_fields.U[Y]( igU  , jgU-1, kgU  )

                             - m_fvCoeffs.Mom[X].AU[Z][t]( igU, jgU, kgU ) * m_fields.U[Z]( igU  , jgU  , kgU+1)
                             - m_fvCoeffs.Mom[X].AU[Z][p]( igU, jgU, kgU ) * m_fields.U[Z]( igU  , jgU  , kgU  )
                             - m_fvCoeffs.Mom[X].AU[Z][b]( igU, jgU, kgU ) * m_fields.U[Z]( igU  , jgU  , kgU-1);
        }
        const floatType bU = ( m_fvCoeffs.Mom[X].F(igU, jgU, kgU)
                             - m_fvCoeffs.Mom[X].B(igU, jgU, kgU)  

                             - m_fvCoeffs.Mom[X].AU[X][n](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU+1, kgU  )
                             - m_fvCoeffs.Mom[X].AU[X][e](igU, jgU, kgU) * m_fields.U[X]( igU+1, jgU  , kgU  )
                             - m_fvCoeffs.Mom[X].AU[X][s](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU-1, kgU  )
                             - m_fvCoeffs.Mom[X].AU[X][w](igU, jgU, kgU) * m_fields.U[X]( igU-1, jgU  , kgU  )
                             - m_fvCoeffs.Mom[X].AU[X][t](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU  , kgU+1) 
                             - m_fvCoeffs.Mom[X].AU[X][b](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU  , kgU-1)

                             - m_fvCoeffs.Mom[X].AP[sUP::cLeft ](igU) * m_fields.P( igU + sUP::iLeft , jgU, kgU)
                             - m_fvCoeffs.Mom[X].AP[sUP::cRight](igU) * m_fields.P( igU + sUP::iRight, jgU, kgU) 

                             + newtonStencilX

                             ) * m_fvCoeffs.Mom[X].diagCoeffInv(igU, jgU, kgU);
    


        // V momentum
        floatType newtonStencilY = 0.0f;
        if constexpr ( LI == Linearisation::Newton ) {
            newtonStencilY = - m_fvCoeffs.Mom[Y].AU[X][e]( igV, jgV, kgV ) * m_fields.U[X]( igV+1, jgV  , kgV  )
                             - m_fvCoeffs.Mom[Y].AU[X][p]( igV, jgV, kgV ) * m_fields.U[X]( igV  , jgV  , kgV  )
                             - m_fvCoeffs.Mom[Y].AU[X][w]( igV, jgV, kgV ) * m_fields.U[X]( igV-1, jgV  , kgV )

                             - m_fvCoeffs.Mom[Y].AU[Z][t]( igV, jgV, kgV ) * m_fields.U[Z]( igV  , jgV  , kgV+1)
                             - m_fvCoeffs.Mom[Y].AU[Z][p]( igV, jgV, kgV ) * m_fields.U[Z]( igV  , jgV  , kgV  )
                             - m_fvCoeffs.Mom[Y].AU[Z][b]( igV, jgV, kgV ) * m_fields.U[Z]( igV  , jgV  , kgV-1);
        }
        const floatType bV = ( m_fvCoeffs.Mom[Y].F(igV, jgV, kgV)
                             - m_fvCoeffs.Mom[Y].B(igV, jgV, kgV)

                             - m_fvCoeffs.Mom[Y].AU[Y][n](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV+1, kgV  ) 
                             - m_fvCoeffs.Mom[Y].AU[Y][e](igV, jgV, kgV) * m_fields.U[Y]( igV+1, jgV  , kgV  ) 
                             - m_fvCoeffs.Mom[Y].AU[Y][s](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV-1, kgV  ) 
                             - m_fvCoeffs.Mom[Y].AU[Y][w](igV, jgV, kgV) * m_fields.U[Y]( igV-1, jgV  , kgV  ) 
                             - m_fvCoeffs.Mom[Y].AU[Y][t](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV  , kgV+1) 
                             - m_fvCoeffs.Mom[Y].AU[Y][b](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV  , kgV-1)

                             - m_fvCoeffs.Mom[Y].AP[sVP::cLeft ](jgV) * m_fields.P( igV, jgV + sVP::iLeft , kgV)
                             - m_fvCoeffs.Mom[Y].AP[sVP::cRight](jgV) * m_fields.P( igV, jgV + sVP::iRight, kgV)

                             + newtonStencilY

                             ) * m_fvCoeffs.Mom[Y].diagCoeffInv(igV, jgV, kgV);
        


        // W momentum
        floatType newtonStencilZ = 0.0f;
        if constexpr ( LI == Linearisation::Newton ) {
            newtonStencilZ =  - m_fvCoeffs.Mom[Z].AU[X][e](igW, jgW, kgW) * m_fields.U[X]( igW+1, jgW  , kgW  )
                              - m_fvCoeffs.Mom[Z].AU[X][p](igW, jgW, kgW) * m_fields.U[X]( igW  , jgW  , kgW  )
                              - m_fvCoeffs.Mom[Z].AU[X][w](igW, jgW, kgW) * m_fields.U[X]( igW-1, jgW  , kgW  )

                              - m_fvCoeffs.Mom[Z].AU[Y][n](igW, jgW, kgW) * m_fields.U[Y]( igW  , jgW+1, kgW  )
                              - m_fvCoeffs.Mom[Z].AU[Y][p](igW, jgW, kgW) * m_fields.U[Y]( igW  , jgW  , kgW  )
                              - m_fvCoeffs.Mom[Z].AU[Y][s](igW, jgW, kgW) * m_fields.U[Y]( igW  , jgW-1, kgW  );
        }
        const floatType bW = ( m_fvCoeffs.Mom[Z].F(igW, jgW, kgW)
                             - m_fvCoeffs.Mom[Z].B(igW, jgW, kgW)
                            
                             - m_fvCoeffs.Mom[Z].AU[Z][n](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW+1, kgW  ) 
                             - m_fvCoeffs.Mom[Z].AU[Z][e](igW, jgW, kgW) * m_fields.U[Z]( igW+1, jgW  , kgW  ) 
                             - m_fvCoeffs.Mom[Z].AU[Z][s](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW-1, kgW  ) 
                             - m_fvCoeffs.Mom[Z].AU[Z][w](igW, jgW, kgW) * m_fields.U[Z]( igW-1, jgW  , kgW  ) 
                             - m_fvCoeffs.Mom[Z].AU[Z][t](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW  , kgW+1) 
                             - m_fvCoeffs.Mom[Z].AU[Z][b](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW  , kgW-1)

                             - m_fvCoeffs.Mom[Z].AP[sWP::cLeft ](kgW) * m_fields.P( igW, jgW, kgW + sWP::iLeft ) 
                             - m_fvCoeffs.Mom[Z].AP[sWP::cRight](kgW) * m_fields.P( igW, jgW, kgW + sWP::iRight)

                             + newtonStencilZ

                             ) * m_fvCoeffs.Mom[Z].diagCoeffInv(igW, jgW, kgW);
    


        // Continuity for pressure
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - m_fvCoeffs.Cont.AP[nn](ig, jg, kg) * m_fields.P( ig  , jg+2, kg  )
                                  - m_fvCoeffs.Cont.AP[ee](ig, jg, kg) * m_fields.P( ig+2, jg  , kg  ) 
                                  - m_fvCoeffs.Cont.AP[ss](ig, jg, kg) * m_fields.P( ig  , jg-2, kg  ) 
                                  - m_fvCoeffs.Cont.AP[ww](ig, jg, kg) * m_fields.P( ig-2, jg  , kg  ) 
                                  - m_fvCoeffs.Cont.AP[tt](ig, jg, kg) * m_fields.P( ig  , jg  , kg+2) 
                                  - m_fvCoeffs.Cont.AP[bb](ig, jg, kg) * m_fields.P( ig  , jg  , kg-2);
        }
        const floatType bP = m_fvCoeffs.Cont.F(ig, jg, kg)
                             - m_fvCoeffs.Cont.B(ig, jg, kg)

                             - m_fvCoeffs.Cont.AU[X][sCU::cLeft ](ig) * m_fields.U[X]( ig + sCU::iLeft , jg, kg)
                             - m_fvCoeffs.Cont.AU[X][sCU::cRight](ig) * m_fields.U[X]( ig + sCU::iRight, jg, kg)

                             - m_fvCoeffs.Cont.AU[Y][sCV::cLeft ](jg) * m_fields.U[Y]( ig, jg + sCV::iLeft , kg)
                             - m_fvCoeffs.Cont.AU[Y][sCV::cRight](jg) * m_fields.U[Y]( ig, jg + sCV::iRight, kg)

                             - m_fvCoeffs.Cont.AU[Z][sCW::cLeft ](kg) * m_fields.U[Z]( ig, jg, kg + sCW::iLeft )
                             - m_fvCoeffs.Cont.AU[Z][sCW::cRight](kg) * m_fields.U[Z]( ig, jg, kg + sCW::iRight)

                             - m_fvCoeffs.Cont.AP[n](ig, jg, kg) * m_fields.P( ig  , jg+1, kg  ) 
                             - m_fvCoeffs.Cont.AP[e](ig, jg, kg) * m_fields.P( ig+1, jg  , kg  ) 
                             - m_fvCoeffs.Cont.AP[s](ig, jg, kg) * m_fields.P( ig  , jg-1, kg  ) 
                             - m_fvCoeffs.Cont.AP[w](ig, jg, kg) * m_fields.P( ig-1, jg  , kg  ) 
                             - m_fvCoeffs.Cont.AP[t](ig, jg, kg) * m_fields.P( ig  , jg  , kg+1) 
                             - m_fvCoeffs.Cont.AP[b](ig, jg, kg) * m_fields.P( ig  , jg  , kg-1)

                             + pressureWideStencil;

        UpdateMolecule( i, j, k, bU, bV, bW, bP );

    }


    // Perform molecule update step with masking
    __attribute__((always_inline)) 
    inline void UpdateMolecule( const intType i, 
                                const intType j, 
                                const intType k,
                                const floatType bU, 
                                const floatType bV, 
                                const floatType bW, 
                                const floatType bP )
    {
        using namespace FVT;
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        // For indexing the staggered cells (with ghost)
        const intType igU{ G( i + sCU::iCoupled ) }, jgU{ G( j                 ) }, kgU{ G( k                 ) },
                      igV{ G( i                 ) }, jgV{ G( j + sCV::iCoupled ) }, kgV{ G( k                 ) },
                      igW{ G( i                 ) }, jgW{ G( j                 ) }, kgW{ G( k + sCW::iCoupled ) },
                       ig{ G( i                 ) },  jg{ G( j                 ) },  kg{ G( k                 ) };

        // Mask for each staggered cell
        const floatType maskU = m_mask(igU, jgU, kgU),
                        maskV = m_mask(igV, jgV, kgV),
                        maskW = m_mask(igW, jgW, kgW),
                        maskP = m_mask(ig, jg, kg);

        // Update P from continuity
        const floatType newP = ( 1 - m_fvCoeffs.Cont.relaxation ) * m_fields.P( ig, jg, kg )
                                 + m_fvCoeffs.Cont.relaxation * 
                                   ( bP 
                                   - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](ig) * bU * maskU
                                   - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](jg) * bV * maskV
                                   - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](kg) * bW * maskW
                                   ) * m_K(i, j, k);

        // Update U from momentum
        const floatType newU = ( 1.0f - m_fvCoeffs.Mom[X].relaxation) * m_fields.U[X]( igU, jgU, kgU )
                                + m_fvCoeffs.Mom[X].relaxation * ( bU - m_fvCoeffs.Mom[X].AP[sUP::cCoupled](igU) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[X].diagCoeffInv(igU, jgU, kgU) );
        
        // Update V from momentum
        const floatType newV = ( 1.0f - m_fvCoeffs.Mom[Y].relaxation ) * m_fields.U[Y]( igV, jgV, kgV )
                                + m_fvCoeffs.Mom[Y].relaxation * ( bV - m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jgV) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Y].diagCoeffInv(igV, jgV, kgV) );
        
        // Update W from momentum
        const floatType newW = ( 1.0f - m_fvCoeffs.Mom[Z].relaxation) * m_fields.U[Z]( igW, jgW, kgW ) 
                                + m_fvCoeffs.Mom[Z].relaxation * ( bW - m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kgW) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Z].diagCoeffInv(igW, jgW, kgW) );

        // Pressure  update
        m_fields.P( ig, jg, kg )       = (1.0f - maskP ) * m_fields.P( ig, jg, kg )        +  maskP * newP;

        // Momentum update
        m_fields.U[X]( igU, jgU, kgU ) = (1.0f - maskP * maskU ) * m_fields.U[X]( igU, jgU, kgU )  +  maskP * maskU * newU;
        m_fields.U[Y]( igV, jgV, kgV ) = (1.0f - maskP * maskV ) * m_fields.U[Y]( igV, jgV, kgV )  +  maskP * maskV * newV;
        m_fields.U[Z]( igW, jgW, kgW ) = (1.0f - maskP * maskW ) * m_fields.U[Z]( igW, jgW, kgW )  +  maskP * maskW * newW;

    }



    // Constants which are global to the linear solver
    void UpdateGlobalConstants()
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Axis::ENUMDATA;
        using FVT::G;

        // Staggered indexing for fields
        intType igU, jgU, kgU,
                igV, jgV, kgV,
                igW, jgW, kgW;

        intType iStart = 0,
                iLength = m_ni,

                jStart = 0,
                jLength = m_nj,

                kStart = 0,
                kLength = m_nk;

        for (intType k = kStart; k != kLength; k++) {

            kgU = G(k);
            kgV = G(k);
            kgW = G(k + sCW::iCoupled);

            for (intType j = jStart; j != jLength; j++) {

                jgU = G(j);
                jgV = G(j + sCV::iCoupled);
                jgW = G(j);

                CFD_PRAGMA_VECTORIZE
                for (intType i = iStart; i != iLength; i++) {

                    igU = G(i + sCU::iCoupled);
                    igV = G(i);
                    igW = G(i);

                    const floatType maskU = m_mask(igU, jgU, kgU),
                                    maskV = m_mask(igV, jgV, kgV),
                                    maskW = m_mask(igW, jgW, kgW);

                    m_K(i, j, k) = m_fvCoeffs.Cont.AP[p](G(i, j, k))
                                 - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](G(i)) * m_fvCoeffs.Mom[X].AP[sUP::cCoupled](igU) * m_fvCoeffs.Mom[X].diagCoeffInv(igU, jgU, kgU) * maskU
                                 - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](G(j)) * m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jgV) * m_fvCoeffs.Mom[Y].diagCoeffInv(igV, jgV, kgV) * maskV
                                 - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](G(k)) * m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kgW) * m_fvCoeffs.Mom[Z].diagCoeffInv(igW, jgW, kgW) * maskW;
                    if ( m_K(i, j, k) == 0.0f ) {
                        m_K(i, j, k) = 0.0f;
                    } else {
                        m_K(i, j, k) = 1.0f / m_K(i, j, k);
                    }
                    
                }
            }
        }
    }


private:
    FieldData<Tensor3D> &m_fields;
    const Tensor3D &m_mask;
    const FVCoefficients &m_fvCoeffs;
    const intType m_ni, m_nj, m_nk;
    Tensor3D m_K;

};

}   // end namespace CFD    


#endif // TRIAD_SOLVER