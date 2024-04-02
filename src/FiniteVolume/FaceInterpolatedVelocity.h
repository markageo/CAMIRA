#ifndef CFD_FACE_INTERPOLATED_VELOCITY
#define CFD_FACE_INTERPOLATED_VELOCITY

#include "FiniteVolume.h"
#include "../Tools/FVTools.h"

namespace CFD {


template< AdvectionSchemes advectionScheme,
          intType          advectionDirection,  // +1 for flux in positive direciton, -1 for flux in negative direciton
          intType          ghostDirection = 0 > // +1 for ghost cell on positive side, -1 for ghost cell on negative side, 0 for no ghost cell
floatType FaceInterpolatedVelocity( const Tensor3D &U,
                                    const MomentumEquation &momentumEquation,
                                    const Mesh &mesh,
                                    const Axis::ENUMDATA axis,
                                    const TensorIndex3D &hiIndex,
                                    const TensorIndex3D &loIndex,
                                    const floatType ghostCellValue = 0.0f )
{
    using enum AdvectionSchemes;

    static_assert( (ghostDirection == +1    ) || (ghostDirection == -1    ) || (ghostDirection == 0));
    static_assert( (advectionDirection == +1) || (advectionDirection == -1) );
    if constexpr ( ghostDirection != 0 ) {  // Only schemes with wide stencils will need ghost cells
        static_assert( (advectionScheme == SOU  ) || (advectionScheme == QUICK) );
    } 
       
    using FVT::G;

    floatType advectedVelocity = 0.0f;
    intType fidx = hiIndex[axis]; 
    
    // First order upwind
    if        constexpr ( advectionScheme == Upwind ) {

        if constexpr ( advectionDirection == +1 ) {
            advectedVelocity = U( G(loIndex) );
        } else {
            advectedVelocity = U( G(hiIndex) );
        }

    // Central
    } else if constexpr ( advectionScheme == Central ) {

        advectedVelocity = ( 1.0f - mesh.interpFactors[axis](fidx) ) * U( G(loIndex) )
                         + mesh.interpFactors[axis](fidx)            * U( G(hiIndex) );

    // Second Order Upwind
    } else if constexpr ( advectionScheme == SOU ) {

        if constexpr ( advectionDirection == +1 ) {

            floatType farVelocityValue;
            if constexpr ( ghostDirection == -1 ) {
                farVelocityValue = ghostCellValue;
            } else {
                TensorIndex3D loloIndex = loIndex;
                loloIndex[axis] -= 1;
                farVelocityValue = U( G(loloIndex) );
            }

            advectedVelocity = momentumEquation.positiveFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * U( G(loIndex) ) 
                             + momentumEquation.positiveFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * farVelocityValue;

        } else {

            floatType farVelocityValue;
            if constexpr ( ghostDirection == +1 ) {
                farVelocityValue = ghostCellValue;
            } else {
                TensorIndex3D hihiIndex = hiIndex;
                hihiIndex[axis] += 1;
                farVelocityValue = U( G(hihiIndex) );
            }
        
            advectedVelocity = momentumEquation.negativeFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * U( G(hiIndex) ) 
                             + momentumEquation.negativeFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * farVelocityValue;

        }

        
    // QUICK
    } else if constexpr ( advectionScheme == QUICK ) {

        if constexpr ( advectionDirection == +1 ) {

            floatType farVelocityValue;
            if constexpr ( ghostDirection == -1 ) {
                farVelocityValue = ghostCellValue;
            } else {
                TensorIndex3D loloIndex = loIndex;
                loloIndex[axis] -= 1;
                farVelocityValue = U( G(loloIndex) );
            }

            advectedVelocity = U( G(loIndex) )
                             + momentumEquation.positiveFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * ( U( G(hiIndex) ) - U( G(loIndex) ) )
                             + momentumEquation.positiveFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * ( U( G(loIndex) ) - farVelocityValue );

        } else {

            floatType farVelocityValue;
            if constexpr ( ghostDirection == +1 ) {
                farVelocityValue = ghostCellValue;
            } else {
                TensorIndex3D hihiIndex = hiIndex;
                hihiIndex[axis] += 1;
                farVelocityValue = U( G(hihiIndex) );
            }

            advectedVelocity = U( G(hiIndex) )
                             + momentumEquation.negativeFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * ( U( G(loIndex) ) - U( G(hiIndex) ) )
                             + momentumEquation.negativeFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * ( U( G(hiIndex) ) - farVelocityValue);
        }
        
    }

    return advectedVelocity;
}


}


#endif  // CFD_FACE_INTERPOLATED_VELOCITY