#ifndef STENCIL_CONSTANTS
#define STENCIL_CONSTANTS

#include "StaggerIndexing.h"

#include "../Types.h"
#include "../Macros.h"
#include "../Tools/FVTools.h"
#include "../FiniteVolume/FiniteVolume.h"

namespace CFD
{

// Precalculate parts of stencil that are constant along a plane
template <TransportCoefficients::ENUMDATA Wstag, 
          MomentumInterpolation MI,
          Linearisation LI>
FieldData<Tensor2D> CalculatePlaneConstants( const intType k,
                                             const FVCoefficients &fvCoeffs,
                                             const FieldData<Tensor3D> &fields )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    using sCW = typename StaggerIndexing< Axis::Z, Wstag >::ContinuityVelocity;
    using sWP = typename StaggerIndexing< Axis::Z, Wstag >::MomentumPressure;

    intType ni = fvCoeffs.nCells[0],
            nj = fvCoeffs.nCells[1];

    // Staggered indices, U and V momentum is not staggered with respect to a plane
    intType kW{ k + sCW::iCoupled }; // W momentum

    // Ghost cells
    intType kgW{ G(kW) };
    intType  kg{ G(k) };

    FieldData<Tensor2D> planeConstants( Tensor2D( fvCoeffs.nCells(Axis::X) + 2*nGhost, fvCoeffs.nCells(Axis::Y) + 2*nGhost ).setZero() );

    for ( intType j = 0; j != nj; j++ ) {

        CFD_PRAGMA_VECTORIZE
        for ( intType i = 0; i != ni; i++ ) {
            intType ig{ G(i) }, jg{ G(j) };

            // U momentum
            floatType newtonStencilX = 0.0f;
            if constexpr ( LI == Linearisation::Newton ) {
                newtonStencilX = - fvCoeffs.Mom[X].AU[Z][sCW::cLeft ](ig, jg, kg) * fields.U[Z]( ig, jg, kg+sCW::iLeft  )
                                 - fvCoeffs.Mom[X].AU[Z][sCW::cRight](ig, jg, kg) * fields.U[Z]( ig, jg, kg+sCW::iRight );
            }
            planeConstants.U[X](ig, jg) = fvCoeffs.Mom[X].F(ig, jg, kg)
                                        - fvCoeffs.Mom[X].B(ig, jg, kg)

                                        - fvCoeffs.Mom[X].AU[X][t](ig, jg, kg) * fields.U[X]( ig  , jg  , kg+1) 
                                        - fvCoeffs.Mom[X].AU[X][b](ig, jg, kg) * fields.U[X]( ig  , jg  , kg-1)
                                        
                                        + newtonStencilX;


            // V momentum
            floatType newtonStencilY = 0.0f;
            if constexpr ( LI == Linearisation::Newton ) {
                newtonStencilY = - fvCoeffs.Mom[Y].AU[Z][sCW::cLeft ](ig, jg, kg) * fields.U[Z]( ig, jg, kg+sCW::iLeft  )
                                 - fvCoeffs.Mom[Y].AU[Z][sCW::cRight](ig, jg, kg) * fields.U[Z]( ig, jg, kg+sCW::iRight );
            }
            planeConstants.U[Y](ig, jg) = fvCoeffs.Mom[Y].F(ig, jg, kg)
                                        - fvCoeffs.Mom[Y].B(ig, jg, kg)

                                        - fvCoeffs.Mom[Y].AU[Y][t](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg+1) 
                                        - fvCoeffs.Mom[Y].AU[Y][b](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg-1)

                                        + newtonStencilY;
                                        

            // W momentum 
            floatType newtonStencilZ = 0.0f;
            if constexpr ( LI == Linearisation::Newton ) {
                newtonStencilZ = - fvCoeffs.Mom[Z].AU[X][e](ig, jg, kgW) * fields.U[X]( ig+1, jg  , kgW  )
                                 - fvCoeffs.Mom[Z].AU[X][p](ig, jg, kgW) * fields.U[X]( ig  , jg  , kgW  )
                                 - fvCoeffs.Mom[Z].AU[X][w](ig, jg, kgW) * fields.U[X]( ig-1, jg  , kgW  )

                                 - fvCoeffs.Mom[Z].AU[Y][n](ig, jg, kgW) * fields.U[Y]( ig  , jg+1, kgW  )
                                 - fvCoeffs.Mom[Z].AU[Y][p](ig, jg, kgW) * fields.U[Y]( ig  , jg  , kgW  )
                                 - fvCoeffs.Mom[Z].AU[Y][s](ig, jg, kgW) * fields.U[Y]( ig  , jg-1, kgW  );
            }
            planeConstants.U[Z](ig, jg) = fvCoeffs.Mom[Z].F(ig, jg, kgW)
                                        - fvCoeffs.Mom[Z].B(ig, jg, kgW)
                            
                                        - fvCoeffs.Mom[Z].AU[Z][t](ig, jg, kgW) * fields.U[Z]( ig  , jg  , kgW+1) 
                                        - fvCoeffs.Mom[Z].AU[Z][b](ig, jg, kgW) * fields.U[Z]( ig  , jg  , kgW-1)

                                        - fvCoeffs.Mom[Z].AP[sWP::cLeft ](kgW) * fields.P( ig, jg, kgW + sWP::iLeft ) 
                                        - fvCoeffs.Mom[Z].AP[sWP::cRight](kgW) * fields.P( ig, jg, kgW + sWP::iRight)

                                        + newtonStencilZ;


            // Continuity equation
            floatType pressureWideStencil = 0.0f;
            if constexpr ( MI == MomentumInterpolation::Implicit ) {
                pressureWideStencil = - fvCoeffs.Cont.AP[tt](ig, jg, kg) * fields.P( ig  , jg  , kg+2) 
                                      - fvCoeffs.Cont.AP[bb](ig, jg, kg) * fields.P( ig  , jg  , kg-2);
            }
            planeConstants.P(ig, jg) = fvCoeffs.Cont.F(ig, jg, kg)
                                     - fvCoeffs.Cont.B(ig, jg, kg)

                                     - fvCoeffs.Cont.AU[Z][sCW::cLeft ](kg) * fields.U[Z]( ig, jg, kg + sCW::iLeft )
                                     - fvCoeffs.Cont.AU[Z][sCW::cRight](kg) * fields.U[Z]( ig, jg, kg + sCW::iRight)

                                     - fvCoeffs.Cont.AP[t](ig, jg, kg) * fields.P( ig  , jg  , kg+1) 
                                     - fvCoeffs.Cont.AP[b](ig, jg, kg) * fields.P( ig  , jg  , kg-1)
    
                                     + pressureWideStencil;
            
        }
    }

    return planeConstants;
}



// Precalculate parts of stencil that are constant along a line
template < TransportCoefficients::ENUMDATA Vstag,
           TransportCoefficients::ENUMDATA Wstag,
           MomentumInterpolation MI,
           Linearisation LI >
FieldData<Tensor1D> CalculateLineConstants( const intType j, 
                                            const intType k, 
                                            const FieldData<Tensor2D> &planeConstants,
                                            const FVCoefficients &fvCoeffs,
                                            const FieldData<Tensor3D> &fields )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    using sCV = typename StaggerIndexing< Axis::Y, Vstag >::ContinuityVelocity;
    using sVP = typename StaggerIndexing< Axis::Y, Vstag >::MomentumPressure;
    using sCW = typename StaggerIndexing< Axis::Z, Wstag >::ContinuityVelocity;

    intType ni = fvCoeffs.nCells[0];

    // Staggered indices, U momenutm is not staggered wrt to the line
    intType jV{ j + sCV::iCoupled }, kV{ k                 }; // V momentum
    intType jW{ j                 }, kW{ k + sCW::iCoupled }; // W momentum

    // Ghost cells
    intType jgV{ G(jV) }, kgV{ G(kV) };
    intType jgW{ G(jW) }, kgW{ G(kW) };
    intType   jg{ G(j) },   kg{ G(k) };

    FieldData<Tensor1D> lineConstants( Tensor1D( fvCoeffs.nCells(Axis::X) + 2*nGhost ).setZero() );

    CFD_PRAGMA_VECTORIZE
    for ( intType i = 0; i != ni; i++ ) {
        intType ig{ G(i) };

        // U momentum
        floatType newtonStencilX = 0.0f;
        if constexpr ( LI == Linearisation::Newton ) {
            newtonStencilX = - fvCoeffs.Mom[X].AU[Y][sCV::cLeft ]( ig, jg, kg ) * fields.U[Y]( ig, jg+sCV::iLeft , kg )
                             - fvCoeffs.Mom[X].AU[Y][sCV::cRight]( ig, jg, kg ) * fields.U[Y]( ig, jg+sCV::iRight, kg );
        }
        lineConstants.U[X](ig) = planeConstants.U[X](ig, jg)
                                + ( 
                                   - fvCoeffs.Mom[X].AU[X][n](ig, jg, kg) * fields.U[X]( ig  , jg+1, kg  )
                                   - fvCoeffs.Mom[X].AU[X][s](ig, jg, kg) * fields.U[X]( ig  , jg-1, kg  )

                                   + newtonStencilX
                                  );

        // V momentum
        floatType newtonStencilY = 0.0f;
        if constexpr ( LI == Linearisation::Newton ) {
            newtonStencilY = - fvCoeffs.Mom[Y].AU[X][e]( ig, jgV, kgV ) * fields.U[X]( ig+1, jgV  , kgV  )
                             - fvCoeffs.Mom[Y].AU[X][p]( ig, jgV, kgV ) * fields.U[X]( ig  , jgV  , kgV  )
                             - fvCoeffs.Mom[Y].AU[X][w]( ig, jgV, kgV ) * fields.U[X]( ig-1, jgV  , kgV  )

                             - fvCoeffs.Mom[Y].AU[Z][sCW::cCoupled]( ig, jgV, kgV ) * fields.U[Z]( ig, jgV, kgV+sCW::iCoupled );
        }
        lineConstants.U[Y](ig) = planeConstants.U[Y](ig, jgV)
                               + (
                                  - fvCoeffs.Mom[Y].AU[Y][n](ig, jgV, kgV) * fields.U[Y]( ig  , jgV+1, kgV  )  
                                  - fvCoeffs.Mom[Y].AU[Y][s](ig, jgV, kgV) * fields.U[Y]( ig  , jgV-1, kgV  ) 

                                  - fvCoeffs.Mom[Y].AP[sVP::cLeft ](jgV) * fields.P( ig, jgV + sVP::iLeft , kgV)
                                  - fvCoeffs.Mom[Y].AP[sVP::cRight](jgV) * fields.P( ig, jgV + sVP::iRight, kgV)

                                  + newtonStencilY
                                 );

        // W momentum
        floatType newtonStencilZ = 0.0f;
        lineConstants.U[Z](ig) = planeConstants.U[Z](ig, jgW)
                                + ( 
                                   - fvCoeffs.Mom[Z].AU[Z][n](ig, jgW, kgW) * fields.U[Z]( ig  , jgW+1, kgW  ) 
                                   - fvCoeffs.Mom[Z].AU[Z][s](ig, jgW, kgW) * fields.U[Z]( ig  , jgW-1, kgW  ) 
                                    
                                   + newtonStencilZ
                                  );

        // Continuity equation
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - fvCoeffs.Cont.AP[nn](ig, jg, kg) * fields.P( ig  , jg+2, kg  )
                                  - fvCoeffs.Cont.AP[ss](ig, jg, kg) * fields.P( ig  , jg-2, kg  );
        }

        lineConstants.P(ig) = planeConstants.P(ig, jg)

                              - fvCoeffs.Cont.AU[Y][sCV::cLeft ](jg) * fields.U[Y]( ig, jg + sCV::iLeft , kg)
                              - fvCoeffs.Cont.AU[Y][sCV::cRight](jg) * fields.U[Y]( ig, jg + sCV::iRight, kg)

                              - fvCoeffs.Cont.AP[n](ig, jg, kg) * fields.P( ig  , jg+1, kg  ) 
                              - fvCoeffs.Cont.AP[s](ig, jg, kg) * fields.P( ig  , jg-1, kg  )
                                
                              + pressureWideStencil;

    }

    return lineConstants;
}


}   // end namespace CFD    


#endif // STENCIL_CONSTANTS