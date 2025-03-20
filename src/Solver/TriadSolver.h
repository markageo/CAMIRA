#ifndef TRIAD_SOLVER
#define TRIAD_SOLVER

#include "StaggerIndexing.h"

#include "../Core/Types.h"
#include "../Core/Macros.h"
#include "../Core/FVTools.h"
#include "../Core/FVLookups.h"
#include "../IO/InputProcessing.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"
#include "../FiniteVolume/FiniteVolume.h"

#include "../IO/ArrayIO.h"


namespace CFD
{

template < TransportCoefficients::ENUMDATA Ustag,
           TransportCoefficients::ENUMDATA Vstag,
           TransportCoefficients::ENUMDATA Wstag, 
           MomentumInterpolation MI >
class TriadSolver
{
    using TC = TransportCoefficients::ENUMDATA;
    using A  = Axis::ENUMDATA;

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
                 const FVCoefficients &fvCoeffs,
                 const InputData::SmootherSettings &smootherSettings ) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_mask( mask ),

                    momX_AU( fvCoeffs.Mom.coeffs[A::X].AU ),
                    momY_AU( fvCoeffs.Mom.coeffs[A::Y].AU ),
                    momZ_AU( fvCoeffs.Mom.coeffs[A::Z].AU ),
                    cont_AP( fvCoeffs.Cont.coeffs.AP ),
                    momX_AP( fvCoeffs.Mom.coeffs[A::X].AP ),
                    momY_AP( fvCoeffs.Mom.coeffs[A::Y].AP ),
                    momZ_AP( fvCoeffs.Mom.coeffs[A::Z].AP ),
                    cont_AUX( fvCoeffs.Cont.coeffs.AU[A::X] ),
                    cont_AUY( fvCoeffs.Cont.coeffs.AU[A::Y] ),
                    cont_AUZ( fvCoeffs.Cont.coeffs.AU[A::Z] ),
                    momX_B( fvCoeffs.Mom.coeffs[A::X].B ),
                    momX_F( fvCoeffs.Mom.coeffs[A::X].F ),
                    momX_AUpInv( fvCoeffs.Mom.coeffs[A::X].diagCoeffInv ),
                    momY_B( fvCoeffs.Mom.coeffs[A::Y].B ),
                    momY_F( fvCoeffs.Mom.coeffs[A::Y].F ),
                    momY_AUpInv( fvCoeffs.Mom.coeffs[A::Y].diagCoeffInv ),
                    momZ_B( fvCoeffs.Mom.coeffs[A::Z].B ),
                    momZ_F( fvCoeffs.Mom.coeffs[A::Z].F ),
                    momZ_AUpInv( fvCoeffs.Mom.coeffs[A::Z].diagCoeffInv ),
                    cont_B( fvCoeffs.Cont.coeffs.B ),
                    cont_F( fvCoeffs.Cont.coeffs.F ),

                    m_ni( fvCoeffs.nCells(0) ),
                    m_nj( fvCoeffs.nCells(1) ),
                    m_nk (fvCoeffs.nCells(2) ),
                    m_K( Tensor3D(m_ni, m_nj, m_nk).setZero() ),
                    m_relaxation( smootherSettings.relaxation )
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
        const floatType bU = ( lineConstants.U[X](igU)  

                             - momX_AU[e](igU, jgU, kgU) * m_fields.U[X]( igU+1, jgU  , kgU  )
                             - momX_AU[w](igU, jgU, kgU) * m_fields.U[X]( igU-1, jgU  , kgU  )

                             - momX_AP[sUP::cLeft ](igU) * m_fields.P( igU + sUP::iLeft , jgU, kgU)
                             - momX_AP[sUP::cRight](igU) * m_fields.P( igU + sUP::iRight, jgU, kgU) 

                             ) * momX_AUpInv(igU, jgU, kgU);


        // V momentum
        const floatType bV = ( lineConstants.U[Y](igV)

                             - momY_AU[e](igV, jgV, kgV) * m_fields.U[Y]( igV+1, jgV  , kgV  ) 
                             - momY_AU[w](igV, jgV, kgV) * m_fields.U[Y]( igV-1, jgV  , kgV  ) 

                             ) * momY_AUpInv(igV, jgV, kgV);


        // W momentum
        const floatType bW = ( lineConstants.U[Z](igW)

                             - momZ_AU[e](igW, jgW, kgW) * m_fields.U[Z]( igW+1, jgW  , kgW  ) 
                             - momZ_AU[w](igW, jgW, kgW) * m_fields.U[Z]( igW-1, jgW  , kgW  ) 

                             ) * momZ_AUpInv(igW, jgW, kgW);


        // Continuity for pressure
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - cont_AP[ee](ig, jg, kg) * m_fields.P( ig+2, jg  , kg  ) 
                                  - cont_AP[ww](ig, jg, kg) * m_fields.P( ig-2, jg  , kg  );
        }
        const floatType bP = lineConstants.P(ig)

                             - cont_AUX[sCU::cLeft ](ig) * m_fields.U[X]( ig + sCU::iLeft , jg, kg)
                             - cont_AUX[sCU::cRight](ig) * m_fields.U[X]( ig + sCU::iRight, jg, kg)

                             - cont_AP[e](ig, jg, kg) * m_fields.P( ig+1, jg  , kg  ) 
                             - cont_AP[w](ig, jg, kg) * m_fields.P( ig-1, jg  , kg  ) 

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
        const floatType bU = ( momX_F(igU, jgU, kgU)
                             - momX_B(igU, jgU, kgU)  

                             - momX_AU[n](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU+1, kgU  )
                             - momX_AU[e](igU, jgU, kgU) * m_fields.U[X]( igU+1, jgU  , kgU  )
                             - momX_AU[s](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU-1, kgU  )
                             - momX_AU[w](igU, jgU, kgU) * m_fields.U[X]( igU-1, jgU  , kgU  )
                             - momX_AU[t](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU  , kgU+1) 
                             - momX_AU[b](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU  , kgU-1)

                             - momX_AP[sUP::cLeft ](igU) * m_fields.P( igU + sUP::iLeft , jgU, kgU)
                             - momX_AP[sUP::cRight](igU) * m_fields.P( igU + sUP::iRight, jgU, kgU) 

                             ) * momX_AUpInv(igU, jgU, kgU);
    


        // V momentum
        const floatType bV = ( momY_F(igV, jgV, kgV)
                             - momY_B(igV, jgV, kgV)

                             - momY_AU[n](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV+1, kgV  ) 
                             - momY_AU[e](igV, jgV, kgV) * m_fields.U[Y]( igV+1, jgV  , kgV  ) 
                             - momY_AU[s](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV-1, kgV  ) 
                             - momY_AU[w](igV, jgV, kgV) * m_fields.U[Y]( igV-1, jgV  , kgV  ) 
                             - momY_AU[t](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV  , kgV+1) 
                             - momY_AU[b](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV  , kgV-1)

                             - momY_AP[sVP::cLeft ](jgV) * m_fields.P( igV, jgV + sVP::iLeft , kgV)
                             - momY_AP[sVP::cRight](jgV) * m_fields.P( igV, jgV + sVP::iRight, kgV)

                             ) * momY_AUpInv(igV, jgV, kgV);
        


        // W momentum
        const floatType bW = ( momZ_F(igW, jgW, kgW)
                             - momZ_B(igW, jgW, kgW)
                            
                             - momZ_AU[n](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW+1, kgW  ) 
                             - momZ_AU[e](igW, jgW, kgW) * m_fields.U[Z]( igW+1, jgW  , kgW  ) 
                             - momZ_AU[s](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW-1, kgW  ) 
                             - momZ_AU[w](igW, jgW, kgW) * m_fields.U[Z]( igW-1, jgW  , kgW  ) 
                             - momZ_AU[t](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW  , kgW+1) 
                             - momZ_AU[b](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW  , kgW-1)

                             - momZ_AP[sWP::cLeft ](kgW) * m_fields.P( igW, jgW, kgW + sWP::iLeft ) 
                             - momZ_AP[sWP::cRight](kgW) * m_fields.P( igW, jgW, kgW + sWP::iRight)

                             ) * momZ_AUpInv(igW, jgW, kgW);
    


        // Continuity for pressure
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - cont_AP[nn](ig, jg, kg) * m_fields.P( ig  , jg+2, kg  )
                                  - cont_AP[ee](ig, jg, kg) * m_fields.P( ig+2, jg  , kg  ) 
                                  - cont_AP[ss](ig, jg, kg) * m_fields.P( ig  , jg-2, kg  ) 
                                  - cont_AP[ww](ig, jg, kg) * m_fields.P( ig-2, jg  , kg  ) 
                                  - cont_AP[tt](ig, jg, kg) * m_fields.P( ig  , jg  , kg+2) 
                                  - cont_AP[bb](ig, jg, kg) * m_fields.P( ig  , jg  , kg-2);
        }
        const floatType bP = cont_F(ig, jg, kg)
                             - cont_B(ig, jg, kg)

                             - cont_AUX[sCU::cLeft ](ig) * m_fields.U[X]( ig + sCU::iLeft , jg, kg)
                             - cont_AUX[sCU::cRight](ig) * m_fields.U[X]( ig + sCU::iRight, jg, kg)

                             - cont_AUY[sCV::cLeft ](jg) * m_fields.U[Y]( ig, jg + sCV::iLeft , kg)
                             - cont_AUY[sCV::cRight](jg) * m_fields.U[Y]( ig, jg + sCV::iRight, kg)

                             - cont_AUZ[sCW::cLeft ](kg) * m_fields.U[Z]( ig, jg, kg + sCW::iLeft )
                             - cont_AUZ[sCW::cRight](kg) * m_fields.U[Z]( ig, jg, kg + sCW::iRight)

                             - cont_AP[n](ig, jg, kg) * m_fields.P( ig  , jg+1, kg  ) 
                             - cont_AP[e](ig, jg, kg) * m_fields.P( ig+1, jg  , kg  ) 
                             - cont_AP[s](ig, jg, kg) * m_fields.P( ig  , jg-1, kg  ) 
                             - cont_AP[w](ig, jg, kg) * m_fields.P( ig-1, jg  , kg  ) 
                             - cont_AP[t](ig, jg, kg) * m_fields.P( ig  , jg  , kg+1) 
                             - cont_AP[b](ig, jg, kg) * m_fields.P( ig  , jg  , kg-1)

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
        const floatType newP = ( 1 - m_relaxation.P ) * m_fields.P( ig, jg, kg )
                                 + m_relaxation.P * 
                                   ( bP 
                                   - cont_AUX[sCU::cCoupled](ig) * ( bU * maskU  +  (1.0f - maskU) * m_fields.U[X]( igU, jgU, kgU ) )   // This masking accounts for IB and ghost cells
                                   - cont_AUY[sCV::cCoupled](jg) * ( bV * maskV  +  (1.0f - maskV) * m_fields.U[Y]( igV, jgV, kgV ) )
                                   - cont_AUZ[sCW::cCoupled](kg) * ( bW * maskW  +  (1.0f - maskW) * m_fields.U[Z]( igW, jgW, kgW ) )
                                   ) * m_K(i, j, k);

        // Pressure  update
        m_fields.P( ig, jg, kg )       = (1.0f - maskP ) * m_fields.P( ig, jg, kg )        +  maskP * newP;

        // Update U from momentum
        const floatType newU = ( 1.0f - m_relaxation.U[X]) * m_fieldsOld.U[X]( igU, jgU, kgU )
                                + m_relaxation.U[X] * ( bU - momX_AP[sUP::cCoupled](igU) * m_fields.P( ig, jg, kg ) * momX_AUpInv(igU, jgU, kgU) );
        
        // Update V from momentum
        const floatType newV = ( 1.0f -m_relaxation.U[Y] ) * m_fieldsOld.U[Y]( igV, jgV, kgV )
                                + m_relaxation.U[Y] * ( bV - momY_AP[sVP::cCoupled](jgV) * m_fields.P( ig, jg, kg ) * momY_AUpInv(igV, jgV, kgV) );
        
        // Update W from momentum
        const floatType newW = ( 1.0f - m_relaxation.U[Z]) * m_fieldsOld.U[Z]( igW, jgW, kgW ) 
                                + m_relaxation.U[Z] * ( bW - momZ_AP[sWP::cCoupled](kgW) * m_fields.P( ig, jg, kg ) * momZ_AUpInv(igW, jgW, kgW) );

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

                    m_K(i, j, k) = cont_AP[p](G(i, j, k))
                                 - cont_AUX[sCU::cCoupled](G(i)) * momX_AP[sUP::cCoupled](igU) * momX_AUpInv(igU, jgU, kgU) * maskU
                                 - cont_AUY[sCV::cCoupled](G(j)) * momY_AP[sVP::cCoupled](jgV) * momY_AUpInv(igV, jgV, kgV) * maskV
                                 - cont_AUZ[sCW::cCoupled](G(k)) * momZ_AP[sWP::cCoupled](kgW) * momZ_AUpInv(igW, jgW, kgW) * maskW;
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
    const FieldData<Tensor3D> &m_fieldsOld;
    const Tensor3D &m_mask;
    // const FVCoefficients &m_fvCoeffs;
    const EnumVector<TransportCoefficients, Tensor3D> &momX_AU,
                                                      &momY_AU,
                                                      &momZ_AU,
                                                      &cont_AP;
    const EnumVector<TransportCoefficients, Tensor1D> &momX_AP, 
                                                      &momY_AP,
                                                      &momZ_AP,
                                                      &cont_AUX, 
                                                      &cont_AUY, 
                                                      &cont_AUZ;
    const Tensor3D &momX_B, &momX_F, &momX_AUpInv,
                   &momY_B, &momY_F, &momY_AUpInv,
                   &momZ_B, &momZ_F, &momZ_AUpInv,
                   &cont_B, &cont_F;
    const intType m_ni, m_nj, m_nk;
    Tensor3D m_K;
    const FieldData<floatType> m_relaxation;

};

}   // end namespace CFD    


#endif // TRIAD_SOLVER