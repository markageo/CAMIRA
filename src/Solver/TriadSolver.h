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


namespace CAMIRA
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

                    m_momX_AU( fvCoeffs.Mom[A::X].AU ),
                    m_momY_AU( fvCoeffs.Mom[A::Y].AU ),
                    m_momZ_AU( fvCoeffs.Mom[A::Z].AU ),
                    m_cont_AP( fvCoeffs.Cont.AP ),
                    m_momX_AP( fvCoeffs.Mom[A::X].AP ),
                    m_momY_AP( fvCoeffs.Mom[A::Y].AP ),
                    m_momZ_AP( fvCoeffs.Mom[A::Z].AP ),
                    m_cont_AUX( fvCoeffs.Cont.AU[A::X] ),
                    m_cont_AUY( fvCoeffs.Cont.AU[A::Y] ),
                    m_cont_AUZ( fvCoeffs.Cont.AU[A::Z] ),
                    m_momX_B( fvCoeffs.Mom[A::X].B ),
                    m_momX_F( fvCoeffs.Mom[A::X].F ),
                    m_momY_B( fvCoeffs.Mom[A::Y].B ),
                    m_momY_F( fvCoeffs.Mom[A::Y].F ),
                    m_momZ_B( fvCoeffs.Mom[A::Z].B ),
                    m_momZ_F( fvCoeffs.Mom[A::Z].F ),
                    m_cont_B( fvCoeffs.Cont.B ),
                    m_cont_F( fvCoeffs.Cont.F ),

                    m_ni( fvCoeffs.nCells(0) ),
                    m_nj( fvCoeffs.nCells(1) ),
                    m_nk (fvCoeffs.nCells(2) ),
                    m_relaxation( smootherSettings.relaxation )
    { UpdateGlobalConstants(); };




    // ------------------------------------ Standard triad solver functions ------------------------------------ //

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

                             - m_momX_AU[e](igU, jgU, kgU) * m_fields.U[X]( igU+1, jgU  , kgU  )
                             - m_momX_AU[w](igU, jgU, kgU) * m_fields.U[X]( igU-1, jgU  , kgU  )

                             - m_momX_AP[sUP::cLeft ](igU) * m_fields.P( igU + sUP::iLeft , jgU, kgU)
                             - m_momX_AP[sUP::cRight](igU) * m_fields.P( igU + sUP::iRight, jgU, kgU) 

                             ) / m_momX_AU[p](igU, jgU, kgU);


        // V momentum
        const floatType bV = ( lineConstants.U[Y](igV)

                             - m_momY_AU[e](igV, jgV, kgV) * m_fields.U[Y]( igV+1, jgV  , kgV  ) 
                             - m_momY_AU[w](igV, jgV, kgV) * m_fields.U[Y]( igV-1, jgV  , kgV  ) 

                             ) / m_momY_AU[p](igV, jgV, kgV);


        // W momentum
        const floatType bW = ( lineConstants.U[Z](igW)

                             - m_momZ_AU[e](igW, jgW, kgW) * m_fields.U[Z]( igW+1, jgW  , kgW  ) 
                             - m_momZ_AU[w](igW, jgW, kgW) * m_fields.U[Z]( igW-1, jgW  , kgW  ) 

                             ) / m_momZ_AU[p](igW, jgW, kgW);


        // Continuity for pressure
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - m_cont_AP[ee](ig, jg, kg) * m_fields.P( ig+2, jg  , kg  ) 
                                  - m_cont_AP[ww](ig, jg, kg) * m_fields.P( ig-2, jg  , kg  );
        }
        const floatType bP = lineConstants.P(ig)

                             - m_cont_AUX[sCU::cLeft ](ig) * m_fields.U[X]( ig + sCU::iLeft , jg, kg)
                             - m_cont_AUX[sCU::cRight](ig) * m_fields.U[X]( ig + sCU::iRight, jg, kg)

                             - m_cont_AP[e](ig, jg, kg) * m_fields.P( ig+1, jg  , kg  ) 
                             - m_cont_AP[w](ig, jg, kg) * m_fields.P( ig-1, jg  , kg  ) 

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
        const floatType bU = ( m_momX_F(igU, jgU, kgU)
                             - m_momX_B(igU, jgU, kgU)  

                             - m_momX_AU[n](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU+1, kgU  )
                             - m_momX_AU[e](igU, jgU, kgU) * m_fields.U[X]( igU+1, jgU  , kgU  )
                             - m_momX_AU[s](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU-1, kgU  )
                             - m_momX_AU[w](igU, jgU, kgU) * m_fields.U[X]( igU-1, jgU  , kgU  )
                             - m_momX_AU[t](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU  , kgU+1) 
                             - m_momX_AU[b](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU  , kgU-1)

                             - m_momX_AP[sUP::cLeft ](igU) * m_fields.P( igU + sUP::iLeft , jgU, kgU)
                             - m_momX_AP[sUP::cRight](igU) * m_fields.P( igU + sUP::iRight, jgU, kgU) 

                             ) / m_momX_AU[p](igU, jgU, kgU); 
    


        // V momentum
        const floatType bV = ( m_momY_F(igV, jgV, kgV)
                             - m_momY_B(igV, jgV, kgV)

                             - m_momY_AU[n](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV+1, kgV  ) 
                             - m_momY_AU[e](igV, jgV, kgV) * m_fields.U[Y]( igV+1, jgV  , kgV  ) 
                             - m_momY_AU[s](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV-1, kgV  ) 
                             - m_momY_AU[w](igV, jgV, kgV) * m_fields.U[Y]( igV-1, jgV  , kgV  ) 
                             - m_momY_AU[t](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV  , kgV+1) 
                             - m_momY_AU[b](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV  , kgV-1)

                             - m_momY_AP[sVP::cLeft ](jgV) * m_fields.P( igV, jgV + sVP::iLeft , kgV)
                             - m_momY_AP[sVP::cRight](jgV) * m_fields.P( igV, jgV + sVP::iRight, kgV)

                            ) / m_momY_AU[p](igV, jgV, kgV);
        


        // W momentum
        const floatType bW = ( m_momZ_F(igW, jgW, kgW)
                             - m_momZ_B(igW, jgW, kgW)
                            
                             - m_momZ_AU[n](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW+1, kgW  ) 
                             - m_momZ_AU[e](igW, jgW, kgW) * m_fields.U[Z]( igW+1, jgW  , kgW  ) 
                             - m_momZ_AU[s](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW-1, kgW  ) 
                             - m_momZ_AU[w](igW, jgW, kgW) * m_fields.U[Z]( igW-1, jgW  , kgW  ) 
                             - m_momZ_AU[t](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW  , kgW+1) 
                             - m_momZ_AU[b](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW  , kgW-1)

                             - m_momZ_AP[sWP::cLeft ](kgW) * m_fields.P( igW, jgW, kgW + sWP::iLeft ) 
                             - m_momZ_AP[sWP::cRight](kgW) * m_fields.P( igW, jgW, kgW + sWP::iRight)

                             ) / m_momZ_AU[p](igW, jgW, kgW);
    


        // Continuity for pressure
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - m_cont_AP[nn](ig, jg, kg) * m_fields.P( ig  , jg+2, kg  )
                                  - m_cont_AP[ee](ig, jg, kg) * m_fields.P( ig+2, jg  , kg  ) 
                                  - m_cont_AP[ss](ig, jg, kg) * m_fields.P( ig  , jg-2, kg  ) 
                                  - m_cont_AP[ww](ig, jg, kg) * m_fields.P( ig-2, jg  , kg  ) 
                                  - m_cont_AP[tt](ig, jg, kg) * m_fields.P( ig  , jg  , kg+2) 
                                  - m_cont_AP[bb](ig, jg, kg) * m_fields.P( ig  , jg  , kg-2);
        }
        const floatType bP = m_cont_F(ig, jg, kg)
                           - m_cont_B(ig, jg, kg)

                           - m_cont_AUX[sCU::cLeft ](ig) * m_fields.U[X]( ig + sCU::iLeft , jg, kg)
                           - m_cont_AUX[sCU::cRight](ig) * m_fields.U[X]( ig + sCU::iRight, jg, kg)

                           - m_cont_AUY[sCV::cLeft ](jg) * m_fields.U[Y]( ig, jg + sCV::iLeft , kg)
                           - m_cont_AUY[sCV::cRight](jg) * m_fields.U[Y]( ig, jg + sCV::iRight, kg)

                           - m_cont_AUZ[sCW::cLeft ](kg) * m_fields.U[Z]( ig, jg, kg + sCW::iLeft )
                           - m_cont_AUZ[sCW::cRight](kg) * m_fields.U[Z]( ig, jg, kg + sCW::iRight)

                           - m_cont_AP[n](ig, jg, kg) * m_fields.P( ig  , jg+1, kg  ) 
                           - m_cont_AP[e](ig, jg, kg) * m_fields.P( ig+1, jg  , kg  ) 
                           - m_cont_AP[s](ig, jg, kg) * m_fields.P( ig  , jg-1, kg  ) 
                           - m_cont_AP[w](ig, jg, kg) * m_fields.P( ig-1, jg  , kg  ) 
                           - m_cont_AP[t](ig, jg, kg) * m_fields.P( ig  , jg  , kg+1) 
                           - m_cont_AP[b](ig, jg, kg) * m_fields.P( ig  , jg  , kg-1)

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
        floatType K = m_cont_AP[p](ig, jg, kg)
                    - m_cont_AUX[sCU::cCoupled](ig) * m_momX_AP[sUP::cCoupled](igU) / m_momX_AU[p](igU, jgU, kgU) * maskU
                    - m_cont_AUY[sCV::cCoupled](jg) * m_momY_AP[sVP::cCoupled](jgV) / m_momY_AU[p](igV, jgV, kgV) * maskV
                    - m_cont_AUZ[sCW::cCoupled](kg) * m_momZ_AP[sWP::cCoupled](kgW) / m_momZ_AU[p](igW, jgW, kgW) * maskW;
        if ( K != 0.0f )
            K = 1.0f / K;

        const floatType newP = ( 1 - m_relaxation.P ) * m_fields.P( ig, jg, kg )
                                 + m_relaxation.P * 
                                   ( bP 
                                   - m_cont_AUX[sCU::cCoupled](ig) * ( bU * maskU  +  (1.0f - maskU) * m_fields.U[X]( igU, jgU, kgU ) )   // This masking accounts for IB and ghost cells
                                   - m_cont_AUY[sCV::cCoupled](jg) * ( bV * maskV  +  (1.0f - maskV) * m_fields.U[Y]( igV, jgV, kgV ) )
                                   - m_cont_AUZ[sCW::cCoupled](kg) * ( bW * maskW  +  (1.0f - maskW) * m_fields.U[Z]( igW, jgW, kgW ) )
                                   ) * K;

        // Pressure  update
        m_fields.P( ig, jg, kg )       = (1.0f - maskP ) * m_fields.P( ig, jg, kg )        +  maskP * newP;

        // Update U from momentum
        const floatType newU = ( 1.0f - m_relaxation.U[X]) * m_fieldsOld.U[X]( igU, jgU, kgU )
                                + m_relaxation.U[X] * ( bU - m_momX_AP[sUP::cCoupled](igU) * m_fields.P( ig, jg, kg ) / m_momX_AU[p](igU, jgU, kgU) );
        
        // Update V from momentum
        const floatType newV = ( 1.0f - m_relaxation.U[Y] ) * m_fieldsOld.U[Y]( igV, jgV, kgV )
                                + m_relaxation.U[Y] * ( bV - m_momY_AP[sVP::cCoupled](jgV) * m_fields.P( ig, jg, kg ) / m_momY_AU[p](igV, jgV, kgV) );
        
        // Update W from momentum
        const floatType newW = ( 1.0f - m_relaxation.U[Z]) * m_fieldsOld.U[Z]( igW, jgW, kgW ) 
                                + m_relaxation.U[Z] * ( bW - m_momZ_AP[sWP::cCoupled](kgW) * m_fields.P( ig, jg, kg ) / m_momZ_AU[p](igW, jgW, kgW));

        // Momentum update
        m_fields.U[X]( igU, jgU, kgU ) = (1.0f - maskU ) * m_fields.U[X]( igU, jgU, kgU )  +  maskU * newU;
        m_fields.U[Y]( igV, jgV, kgV ) = (1.0f - maskV ) * m_fields.U[Y]( igV, jgV, kgV )  +  maskV * newV;
        m_fields.U[Z]( igW, jgW, kgW ) = (1.0f - maskW ) * m_fields.U[Z]( igW, jgW, kgW )  +  maskW * newW;

    }



    // -------------------------------------- Full triad solver functions -------------------------------------- //

    // Core function which updates the local coupled system. Templated by staggering direction.
    // Velocities at central cell are also coupled into the molecule
    // This function evaluates the full stencil
    __attribute__((always_inline)) 
    inline void UpdateTriadFull( const intType i, 
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
        // U momentum - staggered
        const floatType bUS = ( m_momX_F(igU, jgU, kgU)
                              - m_momX_B(igU, jgU, kgU)  

                              - m_momX_AU[n](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU+1, kgU  )
                              - m_momX_AU[e](igU, jgU, kgU) * m_fields.U[X]( igU+1, jgU  , kgU  )
                              - m_momX_AU[s](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU-1, kgU  )
                              - m_momX_AU[w](igU, jgU, kgU) * m_fields.U[X]( igU-1, jgU  , kgU  )
                              - m_momX_AU[t](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU  , kgU+1) 
                              - m_momX_AU[b](igU, jgU, kgU) * m_fields.U[X]( igU  , jgU  , kgU-1)

                              - m_momX_AP[sUP::cLeft ](igU) * m_fields.P( igU + sUP::iLeft , jgU, kgU)
                              - m_momX_AP[sUP::cRight](igU) * m_fields.P( igU + sUP::iRight, jgU, kgU) 

                              ) / m_momX_AU[p](igU, jgU, kgU); 
    

        // U momentum - central
        const floatType bUC = ( m_momX_F(ig, jg, kg)
                              - m_momX_B(ig, jg, kg)  

                              - m_momX_AU[n](ig, jg, kg) * m_fields.U[X]( ig  , jg+1, kg  )
                              - m_momX_AU[e](ig, jg, kg) * m_fields.U[X]( ig+1, jg  , kg  )
                              - m_momX_AU[s](ig, jg, kg) * m_fields.U[X]( ig  , jg-1, kg  )
                              - m_momX_AU[w](ig, jg, kg) * m_fields.U[X]( ig-1, jg  , kg  )
                              - m_momX_AU[t](ig, jg, kg) * m_fields.U[X]( ig  , jg  , kg+1) 
                              - m_momX_AU[b](ig, jg, kg) * m_fields.U[X]( ig  , jg  , kg-1)

                              - m_momX_AP[e](ig) * m_fields.P( ig+1, jg, kg)
                              - m_momX_AP[w](ig) * m_fields.P( ig-1, jg, kg) 

                              ) / m_momX_AU[p](ig, jg, kg); 



        // V momentum - staggered
        const floatType bVS = ( m_momY_F(igV, jgV, kgV)
                              - m_momY_B(igV, jgV, kgV)

                              - m_momY_AU[n](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV+1, kgV  ) 
                              - m_momY_AU[e](igV, jgV, kgV) * m_fields.U[Y]( igV+1, jgV  , kgV  ) 
                              - m_momY_AU[s](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV-1, kgV  ) 
                              - m_momY_AU[w](igV, jgV, kgV) * m_fields.U[Y]( igV-1, jgV  , kgV  ) 
                              - m_momY_AU[t](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV  , kgV+1) 
                              - m_momY_AU[b](igV, jgV, kgV) * m_fields.U[Y]( igV  , jgV  , kgV-1)

                              - m_momY_AP[sVP::cLeft ](jgV) * m_fields.P( igV, jgV + sVP::iLeft , kgV)
                              - m_momY_AP[sVP::cRight](jgV) * m_fields.P( igV, jgV + sVP::iRight, kgV)

                             ) / m_momY_AU[p](igV, jgV, kgV);


        // V momentum - central
        const floatType bVC = ( m_momY_F(ig, jg, kg)
                              - m_momY_B(ig, jg, kg)

                              - m_momY_AU[n](ig, jg, kg) * m_fields.U[Y]( ig  , jg+1, kg  ) 
                              - m_momY_AU[e](ig, jg, kg) * m_fields.U[Y]( ig+1, jg  , kg  ) 
                              - m_momY_AU[s](ig, jg, kg) * m_fields.U[Y]( ig  , jg-1, kg  ) 
                              - m_momY_AU[w](ig, jg, kg) * m_fields.U[Y]( ig-1, jg  , kg  ) 
                              - m_momY_AU[t](ig, jg, kg) * m_fields.U[Y]( ig  , jg  , kg+1) 
                              - m_momY_AU[b](ig, jg, kg) * m_fields.U[Y]( ig  , jg  , kg-1)

                              - m_momY_AP[n](jg) * m_fields.P( ig, jg+1, kg)
                              - m_momY_AP[s](jg) * m_fields.P( ig, jg-1, kg)

                             ) / m_momY_AU[p](ig, jg, kg);
        


        // W momentum - staggered
        const floatType bWS = ( m_momZ_F(igW, jgW, kgW)
                              - m_momZ_B(igW, jgW, kgW)
                            
                              - m_momZ_AU[n](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW+1, kgW  ) 
                              - m_momZ_AU[e](igW, jgW, kgW) * m_fields.U[Z]( igW+1, jgW  , kgW  ) 
                              - m_momZ_AU[s](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW-1, kgW  ) 
                              - m_momZ_AU[w](igW, jgW, kgW) * m_fields.U[Z]( igW-1, jgW  , kgW  ) 
                              - m_momZ_AU[t](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW  , kgW+1) 
                              - m_momZ_AU[b](igW, jgW, kgW) * m_fields.U[Z]( igW  , jgW  , kgW-1)

                              - m_momZ_AP[sWP::cLeft ](kgW) * m_fields.P( igW, jgW, kgW + sWP::iLeft ) 
                              - m_momZ_AP[sWP::cRight](kgW) * m_fields.P( igW, jgW, kgW + sWP::iRight)

                              ) / m_momZ_AU[p](igW, jgW, kgW);


        // W momentum - central
        const floatType bWC = ( m_momZ_F(ig, jg, kg)
                              - m_momZ_B(ig, jg, kg)
                            
                              - m_momZ_AU[n](ig, jg, kg) * m_fields.U[Z]( ig  , jg+1, kg  ) 
                              - m_momZ_AU[e](ig, jg, kg) * m_fields.U[Z]( ig+1, jg  , kg  ) 
                              - m_momZ_AU[s](ig, jg, kg) * m_fields.U[Z]( ig  , jg-1, kg  ) 
                              - m_momZ_AU[w](ig, jg, kg) * m_fields.U[Z]( ig-1, jg  , kg  ) 
                              - m_momZ_AU[t](ig, jg, kg) * m_fields.U[Z]( ig  , jg  , kg+1) 
                              - m_momZ_AU[b](ig, jg, kg) * m_fields.U[Z]( ig  , jg  , kg-1)

                              - m_momZ_AP[t](kg) * m_fields.P( ig, jg, kg+1 ) 
                              - m_momZ_AP[b](kg) * m_fields.P( ig, jg, kg-1 )

                              ) / m_momZ_AU[p](ig, jg, kg);
    


        // Continuity for pressure
        floatType pressureWideStencil = 0.0f;
        if constexpr ( MI == MomentumInterpolation::Implicit ) {
            pressureWideStencil = - m_cont_AP[nn](ig, jg, kg) * m_fields.P( ig  , jg+2, kg  )
                                  - m_cont_AP[ee](ig, jg, kg) * m_fields.P( ig+2, jg  , kg  ) 
                                  - m_cont_AP[ss](ig, jg, kg) * m_fields.P( ig  , jg-2, kg  ) 
                                  - m_cont_AP[ww](ig, jg, kg) * m_fields.P( ig-2, jg  , kg  ) 
                                  - m_cont_AP[tt](ig, jg, kg) * m_fields.P( ig  , jg  , kg+2) 
                                  - m_cont_AP[bb](ig, jg, kg) * m_fields.P( ig  , jg  , kg-2);
        }
        const floatType bP = m_cont_F(ig, jg, kg)
                           - m_cont_B(ig, jg, kg)

                           - m_cont_AUX[sCU::cOpposite ](ig) * m_fields.U[X]( ig + sCU::iOpposite, jg, kg ) 

                           - m_cont_AUY[sCV::cOpposite ](jg) * m_fields.U[Y]( ig, jg + sCV::iOpposite, kg )

                           - m_cont_AUZ[sCW::cOpposite ](kg) * m_fields.U[Z]( ig, jg, kg + sCW::iOpposite )

                           - m_cont_AP[n](ig, jg, kg) * m_fields.P( ig  , jg+1, kg  ) 
                           - m_cont_AP[e](ig, jg, kg) * m_fields.P( ig+1, jg  , kg  ) 
                           - m_cont_AP[s](ig, jg, kg) * m_fields.P( ig  , jg-1, kg  ) 
                           - m_cont_AP[w](ig, jg, kg) * m_fields.P( ig-1, jg  , kg  ) 
                           - m_cont_AP[t](ig, jg, kg) * m_fields.P( ig  , jg  , kg+1) 
                           - m_cont_AP[b](ig, jg, kg) * m_fields.P( ig  , jg  , kg-1)

                           + pressureWideStencil;

        UpdateMoleculeFull( i, j, k, bUS, bUC, bVS, bVC, bWS, bWC, bP );

    }


    // Perform molecule update step with masking
    __attribute__((always_inline)) 
    inline void UpdateMoleculeFull( const intType i, 
                                    const intType j, 
                                    const intType k,
                                    const floatType bUS,
                                    const floatType bUC, 
                                    const floatType bVS,
                                    const floatType bVC, 
                                    const floatType bWS,
                                    const floatType bWC, 
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
        floatType K = m_cont_AP[p](ig, jg, kg)
                    - m_cont_AUX[sCU::cCoupled](ig) * m_momX_AP[sUP::cCoupled](igU) / m_momX_AU[p](igU, jgU, kgU) * maskU
                    - m_cont_AUX[p            ](ig) * m_momX_AP[p            ](ig ) / m_momX_AU[p](ig , jg,  kg ) * maskP

                    - m_cont_AUY[sCV::cCoupled](jg) * m_momY_AP[sVP::cCoupled](jgV) / m_momY_AU[p](igV, jgV, kgV) * maskV
                    - m_cont_AUY[p            ](jg) * m_momY_AP[p            ](jg ) / m_momY_AU[p](ig , jg , kg ) * maskP

                    - m_cont_AUZ[sCW::cCoupled](kg) * m_momZ_AP[sWP::cCoupled](kgW) / m_momZ_AU[p](igW, jgW, kgW) * maskW
                    - m_cont_AUZ[p            ](kg) * m_momZ_AP[p            ](kg ) / m_momZ_AU[p](ig , jg , kg ) * maskP;
        if ( K != 0.0f )
            K = 1.0f / K;

        const floatType newP = ( 1 - m_relaxation.P ) * m_fields.P( ig, jg, kg )
                                 + m_relaxation.P * 
                                   ( bP 
                                   - m_cont_AUX[sCU::cCoupled](ig) * ( bUS * maskU  +  (1.0f - maskU) * m_fields.U[X]( igU, jgU, kgU ) )   // This masking accounts for IB and ghost cells
                                   - m_cont_AUX[p            ](ig) * ( bUC * maskP  +  (1.0f - maskP) * m_fields.U[X]( ig , jg , kg  ) ) 

                                   - m_cont_AUY[sCV::cCoupled](jg) * ( bVS * maskV  +  (1.0f - maskV) * m_fields.U[Y]( igV, jgV, kgV ) )
                                   - m_cont_AUY[p            ](jg) * ( bVC * maskP  +  (1.0f - maskP) * m_fields.U[Y]( ig , jg , kg  ) )

                                   - m_cont_AUZ[sCW::cCoupled](kg) * ( bWS * maskW  +  (1.0f - maskW) * m_fields.U[Z]( igW, jgW, kgW ) )
                                   - m_cont_AUZ[p            ](kg) * ( bWC * maskP  +  (1.0f - maskP) * m_fields.U[Z]( ig , jg , kg  ) )
                                   ) * K;


        // Pressure  update
        m_fields.P( ig, jg, kg )       = (1.0f - maskP ) * m_fields.P( ig, jg, kg )        +  maskP * newP;

        // Update U from momentum
        const floatType newUS = ( 1.0f - m_relaxation.U[X]) * m_fieldsOld.U[X]( igU, jgU, kgU )
                                 + m_relaxation.U[X] * ( bUS - m_momX_AP[sUP::cCoupled](igU) * m_fields.P( ig, jg, kg ) / m_momX_AU[p](igU, jgU, kgU) );

        const floatType newUC = ( 1.0f - m_relaxation.U[X]) * m_fieldsOld.U[X]( ig , jg , kg  )
                                 + m_relaxation.U[X] * ( bUC - m_momX_AP[p            ](ig ) * m_fields.P( ig, jg, kg ) / m_momX_AU[p](ig , jg , kg ) );
        

        // Update V from momentum
        const floatType newVS = ( 1.0f - m_relaxation.U[Y] ) * m_fieldsOld.U[Y]( igV, jgV, kgV )
                                 + m_relaxation.U[Y] * ( bVS - m_momY_AP[sVP::cCoupled](jgV) * m_fields.P( ig, jg, kg ) / m_momY_AU[p](igV, jgV, kgV) );

        const floatType newVC = ( 1.0f - m_relaxation.U[Y] ) * m_fieldsOld.U[Y]( ig , jg , kg  )
                                 + m_relaxation.U[Y] * ( bVC - m_momY_AP[p            ](jg ) * m_fields.P( ig, jg, kg ) / m_momY_AU[p](ig , jg , kg ) );
        

        // Update W from momentum
        const floatType newWS = ( 1.0f - m_relaxation.U[Z]) * m_fieldsOld.U[Z]( igW, jgW, kgW ) 
                                 + m_relaxation.U[Z] * ( bWS - m_momZ_AP[sWP::cCoupled](kgW) * m_fields.P( ig, jg, kg ) / m_momZ_AU[p](igW, jgW, kgW));

        const floatType newWC = ( 1.0f - m_relaxation.U[Z]) * m_fieldsOld.U[Z]( ig , jg , kg  ) 
                                 + m_relaxation.U[Z] * ( bWC - m_momZ_AP[p            ](kg ) * m_fields.P( ig, jg, kg ) / m_momZ_AU[p](ig , jg , kg ));


        // Momentum update
        m_fields.U[X]( igU, jgU, kgU ) = (1.0f - maskU ) * m_fields.U[X]( igU, jgU, kgU )  +  maskU * newUS;
        m_fields.U[X]( ig , jg , kg  ) = (1.0f - maskP ) * m_fields.U[X]( ig , jg , kg  )  +  maskP * newUC;

        m_fields.U[Y]( igV, jgV, kgV ) = (1.0f - maskV ) * m_fields.U[Y]( igV, jgV, kgV )  +  maskV * newVS;
        m_fields.U[Y]( ig , jg , kg  ) = (1.0f - maskP ) * m_fields.U[Y]( ig , jg , kg  )  +  maskP * newVC;

        m_fields.U[Z]( igW, jgW, kgW ) = (1.0f - maskW ) * m_fields.U[Z]( igW, jgW, kgW )  +  maskW * newWS;
        m_fields.U[Z]( ig , jg , kg  ) = (1.0f - maskP ) * m_fields.U[Z]( ig , jg , kg  )  +  maskP * newWC;

    }



    // ----------------------------------------- Other solver functions ---------------------------------------- //

    // Constants which are global to the linear solver
    void UpdateGlobalConstants()
    {
        /* NULL */
    }


private:
    FieldData<Tensor3D> &m_fields;
    const FieldData<Tensor3D> &m_fieldsOld;
    const Tensor3D &m_mask;
    // const FVCoefficients &m_fvCoeffs;
    const EnumVector<TransportCoefficients, Tensor3D> &m_momX_AU,
                                                      &m_momY_AU,
                                                      &m_momZ_AU,
                                                      &m_cont_AP;
    const EnumVector<TransportCoefficients, Tensor1D> &m_momX_AP, 
                                                      &m_momY_AP,
                                                      &m_momZ_AP,
                                                      &m_cont_AUX, 
                                                      &m_cont_AUY, 
                                                      &m_cont_AUZ;
    const Tensor3D &m_momX_B, &m_momX_F, 
                   &m_momY_B, &m_momY_F, 
                   &m_momZ_B, &m_momZ_F, 
                   &m_cont_B, &m_cont_F;
    const intType m_ni, m_nj, m_nk;
    const FieldData<floatType> m_relaxation;

};

}   // end namespace CAMIRA    


#endif // TRIAD_SOLVER