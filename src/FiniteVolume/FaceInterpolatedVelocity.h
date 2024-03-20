#ifndef CFD_FACE_INTERPOLATED_VELOCITY
#define CFD_FACE_INTERPOLATED_VELOCITY

#include "FiniteVolume.h"
#include "../Tools/FVTools.h"

namespace CFD {


template< AdvectionSchemes advectionScheme,
          intType          advectionDirection >
floatType FaceInterpolatedVelocity( const Tensor3D &U,
                                    const MomentumEquation &momentumEquation,
                                    const Mesh &mesh,
                                    const Axis::ENUMDATA axis,
                                    const TensorIndex3D &hiIndex,
                                    const TensorIndex3D &loIndex )
{
    static_assert( (advectionDirection == +1) || (advectionDirection == -1) );

    using FVT::G;

    floatType advectedVelocity = 0.0f;
    intType fidx = hiIndex[axis]; 
    
    // First order upwind
    if        constexpr ( advectionScheme == AdvectionSchemes::Upwind ) {

        if constexpr ( advectionDirection == +1 ) {
            advectedVelocity = U( G(loIndex) );
        } else {
            advectedVelocity = U( G(hiIndex) );
        }

    // Central
    } else if constexpr ( advectionScheme == AdvectionSchemes::Central ) {

        advectedVelocity = ( 1.0f - mesh.interpFactors[axis](fidx) ) * U( G(loIndex) )
                         + mesh.interpFactors[axis](fidx)            * U( G(hiIndex) );

    // Second Order Upwind
    } else if constexpr ( advectionScheme == AdvectionSchemes::SOU ) {

        if constexpr ( advectionDirection == +1 ) {

            TensorIndex3D loloIndex = loIndex;
            loloIndex[axis] -= 1;

            advectedVelocity = momentumEquation.positiveFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * U( G(loIndex) ) 
                             + momentumEquation.positiveFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * U( G(loloIndex) );

        } else {

            TensorIndex3D hihiIndex = hiIndex;
            hihiIndex[axis] += 1;

            advectedVelocity = momentumEquation.negativeFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * U( G(hiIndex) ) 
                             + momentumEquation.negativeFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * U( G(hihiIndex) );

        }

        
    // QUICK
    } else if constexpr ( advectionScheme == AdvectionSchemes::QUICK ) {


        if constexpr ( advectionDirection == +1 ) {

            TensorIndex3D loloIndex = loIndex;
            loloIndex[axis] -= 1;

            advectedVelocity = U( G(loIndex) )
                             + momentumEquation.positiveFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * ( U( G(hiIndex) ) - U( G(loIndex) ) )
                             + momentumEquation.positiveFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * ( U( G(loIndex) ) - U( G(loloIndex) ) );

        } else {

            TensorIndex3D hihiIndex = hiIndex;
            hihiIndex[axis] += 1;

            advectedVelocity = U( G(hiIndex) )
                             + momentumEquation.negativeFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * ( U( G(loIndex) ) - U( G(hiIndex) ) )
                             + momentumEquation.negativeFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * ( U( G(hiIndex) ) - U( G(hihiIndex) ) );
        }
        

    }

    return advectedVelocity;
}


}


#endif  // CFD_FACE_INTERPOLATED_VELOCITY