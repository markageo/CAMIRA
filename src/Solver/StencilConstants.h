#ifndef STENCIL_CONSTANTS
#define STENCIL_CONSTANTS

#include "StaggerIndexing.h"

#include "../Core/Types.h"
#include "../Core/Macros.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/FiniteVolume.h"

namespace CFD
{

// Precalculate parts of stencil that are constant along a plane
template <TransportCoefficients::ENUMDATA Wstag, 
          MomentumInterpolation MI >
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
            planeConstants.U[X](ig, jg) = fvCoeffs.Mom.coeffs[X].F(ig, jg, kg)
                                        - fvCoeffs.Mom.coeffs[X].B(ig, jg, kg)

                                        - fvCoeffs.Mom.coeffs[X].AU[t](ig, jg, kg) * fields.U[X]( ig  , jg  , kg+1) 
                                        - fvCoeffs.Mom.coeffs[X].AU[b](ig, jg, kg) * fields.U[X]( ig  , jg  , kg-1);


            // V momentum
            planeConstants.U[Y](ig, jg) = fvCoeffs.Mom.coeffs[Y].F(ig, jg, kg)
                                        - fvCoeffs.Mom.coeffs[Y].B(ig, jg, kg)

                                        - fvCoeffs.Mom.coeffs[Y].AU[t](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg+1) 
                                        - fvCoeffs.Mom.coeffs[Y].AU[b](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg-1);
                                        

            // W momentum 
            planeConstants.U[Z](ig, jg) = fvCoeffs.Mom.coeffs[Z].F(ig, jg, kgW)
                                        - fvCoeffs.Mom.coeffs[Z].B(ig, jg, kgW)
                            
                                        - fvCoeffs.Mom.coeffs[Z].AU[t](ig, jg, kgW) * fields.U[Z]( ig  , jg  , kgW+1) 
                                        - fvCoeffs.Mom.coeffs[Z].AU[b](ig, jg, kgW) * fields.U[Z]( ig  , jg  , kgW-1)

                                        - fvCoeffs.Mom.coeffs[Z].AP[sWP::cLeft ](kgW) * fields.P( ig, jg, kgW + sWP::iLeft ) 
                                        - fvCoeffs.Mom.coeffs[Z].AP[sWP::cRight](kgW) * fields.P( ig, jg, kgW + sWP::iRight);


            // Continuity equation
            floatType pressureWideStencil = 0.0f;
            if constexpr ( MI == MomentumInterpolation::Implicit ) {
                pressureWideStencil = - fvCoeffs.Cont.coeffs.AP[tt](ig, jg, kg) * fields.P( ig  , jg  , kg+2) 
                                      - fvCoeffs.Cont.coeffs.AP[bb](ig, jg, kg) * fields.P( ig  , jg  , kg-2);
            }
            planeConstants.P(ig, jg) = fvCoeffs.Cont.coeffs.F(ig, jg, kg)
                                     - fvCoeffs.Cont.coeffs.B(ig, jg, kg)

                                     - fvCoeffs.Cont.coeffs.AU[Z][sCW::cLeft ](kg) * fields.U[Z]( ig, jg, kg + sCW::iLeft )
                                     - fvCoeffs.Cont.coeffs.AU[Z][sCW::cRight](kg) * fields.U[Z]( ig, jg, kg + sCW::iRight)

                                     - fvCoeffs.Cont.coeffs.AP[t](ig, jg, kg) * fields.P( ig  , jg  , kg+1) 
                                     - fvCoeffs.Cont.coeffs.AP[b](ig, jg, kg) * fields.P( ig  , jg  , kg-1)
    
                                     + pressureWideStencil;
            
        }
    }

    return planeConstants;
}



// Precalculate parts of stencil that are constant along a line
template < TransportCoefficients::ENUMDATA Vstag,
           TransportCoefficients::ENUMDATA Wstag,
           MomentumInterpolation MI >
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
        lineConstants.U[X](ig) = planeConstants.U[X](ig, jg)
                                + ( 
                                   - fvCoeffs.Mom.coeffs[X].AU[n](ig, jg, kg) * fields.U[X]( ig  , jg+1, kg  )
                                   - fvCoeffs.Mom.coeffs[X].AU[s](ig, jg, kg) * fields.U[X]( ig  , jg-1, kg  )
                                  );

        // V momentum
        lineConstants.U[Y](ig) = planeConstants.U[Y](ig, jgV)
                               + (
                                  - fvCoeffs.Mom.coeffs[Y].AU[n](ig, jgV, kgV) * fields.U[Y]( ig  , jgV+1, kgV  )  
                                  - fvCoeffs.Mom.coeffs[Y].AU[s](ig, jgV, kgV) * fields.U[Y]( ig  , jgV-1, kgV  ) 

                                  - fvCoeffs.Mom.coeffs[Y].AP[sVP::cLeft ](jgV) * fields.P( ig, jgV + sVP::iLeft , kgV)
                                  - fvCoeffs.Mom.coeffs[Y].AP[sVP::cRight](jgV) * fields.P( ig, jgV + sVP::iRight, kgV)
                                  
                                 );

        // W momentum
        lineConstants.U[Z](ig) = planeConstants.U[Z](ig, jgW)
                                + ( 
                                   - fvCoeffs.Mom.coeffs[Z].AU[n](ig, jgW, kgW) * fields.U[Z]( ig  , jgW+1, kgW  ) 
                                   - fvCoeffs.Mom.coeffs[Z].AU[s](ig, jgW, kgW) * fields.U[Z]( ig  , jgW-1, kgW  ) 
                                  );

        // Continuity equation
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - fvCoeffs.Cont.coeffs.AP[nn](ig, jg, kg) * fields.P( ig  , jg+2, kg  )
                                  - fvCoeffs.Cont.coeffs.AP[ss](ig, jg, kg) * fields.P( ig  , jg-2, kg  );
        }

        lineConstants.P(ig) = planeConstants.P(ig, jg)

                              - fvCoeffs.Cont.coeffs.AU[Y][sCV::cLeft ](jg) * fields.U[Y]( ig, jg + sCV::iLeft , kg)
                              - fvCoeffs.Cont.coeffs.AU[Y][sCV::cRight](jg) * fields.U[Y]( ig, jg + sCV::iRight, kg)

                              - fvCoeffs.Cont.coeffs.AP[n](ig, jg, kg) * fields.P( ig  , jg+1, kg  ) 
                              - fvCoeffs.Cont.coeffs.AP[s](ig, jg, kg) * fields.P( ig  , jg-1, kg  )
                                
                              + pressureWideStencil;

    }

    return lineConstants;
}


}   // end namespace CFD    


#endif // STENCIL_CONSTANTS