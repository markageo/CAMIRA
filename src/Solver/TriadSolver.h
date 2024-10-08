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
        // floatType masterMask = m_mask(iU, jU, kU) * m_mask(iV, jV, kV) * m_mask(iW, jW, kW) * m_mask(i, j, k);
        floatType masterMask = 1.0f;

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
        floatType bU = ( m_fvCoeffs.Mom[X].F(igU, jgU, kgU)
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
        floatType bV = ( m_fvCoeffs.Mom[Y].F(igV, jgV, kgV)
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
        floatType bW = ( m_fvCoeffs.Mom[Z].F(igW, jgW, kgW)
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
        floatType bP = m_fvCoeffs.Cont.F(ig, jg, kg)
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


        // Only update the molecule if none of the cells are within the immersed boundary
        floatType masterMask = m_mask(igU, jgU, kgU) * m_mask(igV, jgV, kgV) * m_mask(igW, jgW, kgW) * m_mask(ig, jg, kg); 

        // Update P from continuity
        floatType newP = ( 1 - m_fvCoeffs.Cont.relaxation ) * m_fields.P( ig, jg, kg )
                                 + m_fvCoeffs.Cont.relaxation * 
                                   ( bP 
                                   - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](ig) * bU 
                                   - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](jg) * bV 
                                   - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](kg) * bW 
                                   ) * m_K(i, j, k);


        m_fields.P( ig, jg, kg )       = (1.0f - masterMask) * m_fields.P( ig, jg, kg )        +  masterMask * newP;

        // Update U from momentum
        floatType newU = ( 1 - m_fvCoeffs.Mom[X].relaxation) * m_fields.U[X]( igU, jgU, kgU )
                         + m_fvCoeffs.Mom[X].relaxation * ( bU - m_fvCoeffs.Mom[X].AP[sUP::cCoupled](igU) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[X].diagCoeffInv(igU, jgU, kgU) );
        
        // Update V from momentum
        floatType newV = ( 1 - m_fvCoeffs.Mom[Y].relaxation ) * m_fields.U[Y]( igV, jgV, kgV )
                         + m_fvCoeffs.Mom[Y].relaxation * ( bV - m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jgV) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Y].diagCoeffInv(igV, jgV, kgV) );
        

        // Update W from momentum
        floatType newW = ( 1 - m_fvCoeffs.Mom[Z].relaxation) * m_fields.U[Z]( igW, jgW, kgW ) 
                         + m_fvCoeffs.Mom[Z].relaxation * ( bW - m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kgW) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Z].diagCoeffInv(igW, jgW, kgW) );


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
        using FVT::G;

        // Staggered indexing for fields
        intType igU, jgU, kgU,
                igV, jgV, kgV,
                igW, jgW, kgW;

        // // Starting and ending indices, since K cannot be calculated on some boundaries due to the staggering
        // intType iStart = 1 + sCU::iLeft,
        //         iLength = m_ni - 1 + sCU::iRight,

        //         jStart = 1 + sCV::iLeft,
        //         jLength = m_nj - 1 + sCV::iRight,

        //         kStart = 1 + sCW::iLeft,
        //         kLength = m_nk - 1 + sCW::iRight;

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

                    m_K(i, j, k) = m_fvCoeffs.Cont.AP[p](G(i, j, k)) 
                                 - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](G(i)) * m_fvCoeffs.Mom[X].AP[sUP::cCoupled](igU) * m_fvCoeffs.Mom[X].diagCoeffInv(igU, jgU, kgU) 
                                 - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](G(j)) * m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jgV) * m_fvCoeffs.Mom[Y].diagCoeffInv(igV, jgV, kgV) 
                                 - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](G(k)) * m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kgW) * m_fvCoeffs.Mom[Z].diagCoeffInv(igW, jgW, kgW);
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