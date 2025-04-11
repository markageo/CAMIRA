#ifndef CFD_FACE_INTERPOLATED_VELOCITY
#define CFD_FACE_INTERPOLATED_VELOCITY

#include "FiniteVolume.h"
#include "../Core/FVTools.h"
#include "FiniteVolumeStructures.h"
#include <sys/cdefs.h>

namespace CFD {


template< AdvectionSchemes advectionScheme,
          intType          advectionDirection,  // +1 for flux in positive direciton, -1 for flux in negative direciton
          intType          ghostDirection = 0 > // +1 for ghost cell on positive side, -1 for ghost cell on negative side, 0 for no ghost cell
__attribute__((always_inline, flatten))
inline floatType FaceInterpolatedVelocity( const Tensor3D &U,
                                           const FVCoefficients &fvCoeffs,
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
    
    TensorIndex3D hiIndexG = G(hiIndex),
                  loIndexG = G(loIndex);

    // First order upwind
    if        constexpr ( advectionScheme == Upwind ) {

        if constexpr ( advectionDirection == +1 ) {
            advectedVelocity = U( loIndexG );
        } else {
            advectedVelocity = U( hiIndexG );
        }

    // Central
    } else if constexpr ( advectionScheme == Central ) {

        advectedVelocity = ( 1.0f - mesh.interpFactors[axis](fidx) ) * U( loIndexG )
                         + mesh.interpFactors[axis](fidx)            * U( hiIndexG );

    // Second Order Upwind
    } else if constexpr ( advectionScheme == SOU ) {

        if constexpr ( advectionDirection == +1 ) {

            floatType farVelocityValue;
            if constexpr ( ghostDirection == -1 ) {
                farVelocityValue = ghostCellValue;
            } else {
                TensorIndex3D loloIndexG = loIndexG;
                loloIndexG[axis] -= 1;
                farVelocityValue = U( loloIndexG );
            }

            advectedVelocity = fvCoeffs.positiveFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * U( loIndexG ) 
                             + fvCoeffs.positiveFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * farVelocityValue;

        } else {

            floatType farVelocityValue;
            if constexpr ( ghostDirection == +1 ) {
                farVelocityValue = ghostCellValue;
            } else {
                TensorIndex3D hihiIndexG = hiIndexG;
                hihiIndexG[axis] += 1;
                farVelocityValue = U( hihiIndexG );
            }
        
            advectedVelocity = fvCoeffs.negativeFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * U( hiIndexG ) 
                             + fvCoeffs.negativeFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * farVelocityValue;

        }

        
    // QUICK
    } else if constexpr ( advectionScheme == QUICK ) {

        if constexpr ( advectionDirection == +1 ) {

            floatType farVelocityValue;
            if constexpr ( ghostDirection == -1 ) {
                farVelocityValue = ghostCellValue;
            } else {
                TensorIndex3D loloIndexG = loIndexG;
                loloIndexG[axis] -= 1;
                farVelocityValue = U( loloIndexG );
            }

            advectedVelocity = U( loIndexG )
                             + fvCoeffs.positiveFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * ( U( hiIndexG ) - U( loIndexG ) )
                             + fvCoeffs.positiveFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * ( U( loIndexG ) - farVelocityValue );

        } else {

            floatType farVelocityValue;
            if constexpr ( ghostDirection == +1 ) {
                farVelocityValue = ghostCellValue;
            } else {
                TensorIndex3D hihiIndexG = hiIndexG;
                hihiIndexG[axis] += 1;
                farVelocityValue = U( hihiIndexG );
            }

            advectedVelocity = U( hiIndexG )
                             + fvCoeffs.negativeFluxHiOrderAdvectionCoeffs.g1[axis]( fidx ) * ( U( loIndexG ) - U( hiIndexG ) )
                             + fvCoeffs.negativeFluxHiOrderAdvectionCoeffs.g2[axis]( fidx ) * ( U( hiIndexG ) - farVelocityValue);
        }
        
    }

    return advectedVelocity;
}


}


#endif  // CFD_FACE_INTERPOLATED_VELOCITY