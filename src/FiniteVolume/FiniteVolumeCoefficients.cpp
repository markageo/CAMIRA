#include "FiniteVolume.h"
#include "../Core/Macros.h"
#include "../Core/FVTools.h"
#include "../Core/FVLookups.h"
#include "FaceInterpolatedVelocity.h"
#include "../Parallel/Parallel.h"

#include <RAJA/index/RangeSegment.hpp>
#include <RAJA/pattern/kernel.hpp>
#include <RAJA/policy/openmp/kernel/Collapse.hpp>
#include <RAJA/policy/openmp/policy.hpp>
#include <RAJA/policy/sequential/policy.hpp>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>

#include <string>
#include "../IO/ArrayIO.h"
#include "FiniteVolumeFunctions.h"
#include "FiniteVolumeStructures.h"

namespace CFD
{

using namespace FVT;
 
namespace
{


/*---------------------------------------------------------------------------------------------------------------*\
                                                    Diffusion
\*---------------------------------------------------------------------------------------------------------------*/


// Set diffusion coefficients for a given momentum equation
void SetDiffusionCoeffients( FVCoefficients &fvCoeffs, 
                             const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    
    auto &diff = fvCoeffs.diff;

    // Diffusion in each axis is calculated in the same way
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis],
                                        west = LUT::LoCoeff[axis];     

        // Internal faces
        for (intType i = 1; i != mesh.nCells(axis); i++) {
            
            // Cell on west side of face
            diff[axis][p   ](i-1) +=   mesh.cellCenterDiffInv[axis](i);
            diff[axis][east](i-1) += - mesh.cellCenterDiffInv[axis](i);

            // Cell on east side of face
            diff[axis][p   ](i)   +=   mesh.cellCenterDiffInv[axis](i);
            diff[axis][west](i)   += - mesh.cellCenterDiffInv[axis](i);

        }

        // Lo boundary
        diff[axis][p   ](0)   +=   mesh.cellCenterDiffInv[axis](0);
        diff[axis][west](0)   += - mesh.cellCenterDiffInv[axis](0);

        // Hi boundary
        diff[axis][p   ]( mesh.nCells(axis)-1 ) +=   mesh.cellCenterDiffInv[axis]( mesh.nFacesNormal[axis](axis)-1 );
        diff[axis][east]( mesh.nCells(axis)-1 ) += - mesh.cellCenterDiffInv[axis]( mesh.nFacesNormal[axis](axis)-1 );


        // Multiply by inverse cell length
        for (intType i = 0; i != mesh.nCells(axis); i++) {
            diff[axis][p   ](i) *= mesh.cellLengthsInv[axis](i);
            diff[axis][east](i) *= mesh.cellLengthsInv[axis](i);
            diff[axis][west](i) *= mesh.cellLengthsInv[axis](i);
        }

        // Multiply by viscosity
        for (intType i = 0; i != mesh.nCells(axis); i++) {
            diff[axis][p   ](i) *= fvCoeffs.nu;
            diff[axis][east](i) *= fvCoeffs.nu;
            diff[axis][west](i) *= fvCoeffs.nu;
        }

    } );
}



/*---------------------------------------------------------------------------------------------------------------*\
                                           Momentum Advection Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


void SetHighOrderAdvectionCoefficients( FVCoefficients &fvCoeffs,
                                        const Mesh &mesh )
{
    auto &negativeFluxHiOrderAdvectionCoeffs = fvCoeffs.negativeFluxHiOrderAdvectionCoeffs;
    auto &positiveFluxHiOrderAdvectionCoeffs = fvCoeffs.positiveFluxHiOrderAdvectionCoeffs;

    switch ( fvCoeffs.advectionScheme ) {

        case AdvectionSchemes::SOU:
            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                negativeFluxHiOrderAdvectionCoeffs.g1[axis] = Tensor1D( mesh.nFacesNormal[axis][axis] );
                negativeFluxHiOrderAdvectionCoeffs.g2[axis] = Tensor1D( mesh.nFacesNormal[axis][axis] );

                positiveFluxHiOrderAdvectionCoeffs.g1[axis] = Tensor1D( mesh.nFacesNormal[axis][axis] );
                positiveFluxHiOrderAdvectionCoeffs.g2[axis] = Tensor1D( mesh.nFacesNormal[axis][axis] );

                for ( intType i = 1; i != mesh.nFacesNormal[axis][axis] - 1; i++ ) {

                    const floatType xf = mesh.cellFaces[axis]( i );
                    floatType xU = 0.0f, xUU = 0.0f;

                    // Flux in positive direction across the face
                    xU = mesh.cellCenters[axis]( i-1 );
                    if ( i == 1 ) {
                        xUU = xU - mesh.cellLengths[axis]( i-1 );
                    } else {
                        xUU = mesh.cellCenters[axis]( i-2 );
                    }
                    positiveFluxHiOrderAdvectionCoeffs.g1[axis](i) = ( xf - xUU ) / ( xU - xUU );
                    positiveFluxHiOrderAdvectionCoeffs.g2[axis](i) = ( xf - xU  ) / ( xUU - xU );

                    // Flux in negative direction across the face
                    xU  = mesh.cellCenters[axis]( i );
                    if ( i == mesh.nFacesNormal[axis][axis] - 2 ) {
                        xUU = xU + mesh.cellLengths[axis]( i );
                    } else {
                        xUU = mesh.cellCenters[axis]( i+1 );
                    }
                    negativeFluxHiOrderAdvectionCoeffs.g1[axis](i) = ( xf - xUU ) / ( xU - xUU );
                    negativeFluxHiOrderAdvectionCoeffs.g2[axis](i) = ( xf - xU  ) / ( xUU - xU );
                    
                }
            } );
            break;


        case AdvectionSchemes::QUICK:
            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                negativeFluxHiOrderAdvectionCoeffs.g1[axis] = Tensor1D( mesh.nFacesNormal[axis][axis] );
                negativeFluxHiOrderAdvectionCoeffs.g2[axis] = Tensor1D( mesh.nFacesNormal[axis][axis] );

                positiveFluxHiOrderAdvectionCoeffs.g1[axis] = Tensor1D( mesh.nFacesNormal[axis][axis] );
                positiveFluxHiOrderAdvectionCoeffs.g2[axis] = Tensor1D( mesh.nFacesNormal[axis][axis] );

                for ( intType i = 1; i != mesh.nFacesNormal[axis][axis] - 1; i++ ) {

                    const floatType xf = mesh.cellFaces[axis]( i );
                    floatType xU = 0.0f, xD = 0.0f, xUU = 0.0f;

                    // Flux in positive direction across the face
                    xU = mesh.cellCenters[axis]( i-1 ),
                    xD = mesh.cellCenters[axis]( i );
                    if ( i == 1 ) {
                        xUU = xU - mesh.cellLengths[axis]( i-1 );
                    } else {
                        xUU = mesh.cellCenters[axis]( i-2 );
                    }
                    positiveFluxHiOrderAdvectionCoeffs.g1[axis](i) = ( xf - xU ) * ( xf - xUU ) / ( xD - xU  ) / ( xD - xUU );
                    positiveFluxHiOrderAdvectionCoeffs.g2[axis](i) = ( xf - xU ) * ( xD - xf  ) / ( xU - xUU ) / ( xD - xUU );

                    // Flux in negative direction across the face
                    xU = mesh.cellCenters[axis]( i ),
                    xD = mesh.cellCenters[axis]( i-1 );
                    if ( i == mesh.nFacesNormal[axis][axis] - 2 ) {
                        xUU = xU + mesh.cellLengths[axis]( i );
                    } else {
                        xUU = mesh.cellCenters[axis]( i+1 );
                    }
                    negativeFluxHiOrderAdvectionCoeffs.g1[axis](i) = ( xf - xU ) * ( xf - xUU ) / ( xD - xU  ) / ( xD - xUU );
                    negativeFluxHiOrderAdvectionCoeffs.g2[axis](i) = ( xf - xU ) * ( xD - xf  ) / ( xU - xUU ) / ( xD - xUU );
                    
                }
            } );
            break;

        default:
            /* NULL */
            break;
    }
}



// Upwind implicit coefficients with deffered correction for higher order schemes
// Assumes that the 'p' coefficient has been set to zero
template< AdvectionSchemes advectionScheme > 
[[maybe_unused]]
__attribute__((flatten))
void InteriorAdvectionCoefficients( FVCoefficients &fvCoeffs, 
                                    const FieldData<Tensor3D> &fields,
                                    const EnumVector<Axis, Tensor3D> &faceFluxes, 
                                    const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    constexpr bool hasDeferredCorrection = ( advectionScheme != AdvectionSchemes::Upwind );

    auto &AU     = fvCoeffs.Mom[X].AU;  // All equations share the same underlying data
    
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

        const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis], 
                                              west = LUT::LoCoeff[axis];

        for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
            for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
                for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                    
                    TensorIndex3D hiIndex = { i, j, k },
                                  loIndex = { i, j, k };
                    loIndex[axis] -= 1;

                    const floatType faceFlux = faceFluxes[ axis ](i, j, k);
                    const bool fluxIsPositive = ( faceFlux >= 0.0f );

                    if ( fluxIsPositive ) {
                        AU[p   ](G(loIndex)) +=   faceFlux * mesh.cellLengthsInv[axis]( loIndex[axis] );
                        AU[east](G(loIndex))  =   0.0f;
                        AU[p   ](G(hiIndex)) +=   0.0f;
                        AU[west](G(hiIndex))  = - faceFlux * mesh.cellLengthsInv[axis]( hiIndex[axis] );
                    } else {
                        AU[p   ](G(loIndex)) +=   0.0f;
                        AU[east](G(loIndex))  =   faceFlux * mesh.cellLengthsInv[axis]( loIndex[axis] );
                        AU[p   ](G(hiIndex)) += - faceFlux * mesh.cellLengthsInv[axis]( hiIndex[axis] );
                        AU[west](G(hiIndex))  =   0.0f;
                    }


                    if constexpr ( !hasDeferredCorrection )
                        continue;

                    
                    // Source term is different for each momentum equation (due to advected velocity)
                    EnumFor<Axis>( [&] (Axis::ENUMDATA comp) {

                        const auto &U   = fields.U[comp];
                        auto &sourceTerm = fvCoeffs.Mom[comp].B;

                        floatType highOrderAdvectedVelocity{0.0f},
                                  upwindAdvectedVelocity{0.0f};

                        if ( fluxIsPositive ) {
                            highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme         , +1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                            upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, +1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                        } else {
                            highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme         , -1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                            upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, -1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                        }

                        // Deferred correction term
                        sourceTerm(G(loIndex)) +=   fvCoeffs.advectionBlendingFactor 
                                                *   faceFlux
                                                *   ( highOrderAdvectedVelocity - upwindAdvectedVelocity ) 
                                                *   mesh.cellLengthsInv[axis]( loIndex[axis] );
                            
                        sourceTerm(G(hiIndex)) += - fvCoeffs.advectionBlendingFactor 
                                                *   faceFlux
                                                *   ( highOrderAdvectedVelocity - upwindAdvectedVelocity ) 
                                                *   mesh.cellLengthsInv[axis]( hiIndex[axis] );


                    } );

                    
                }
            }
        }

    } );

}



// Direction is a template parameter to allow for compiler optimisations
// Assumes that the 'p' coefficient has been set to zero
template< AdvectionSchemes advectionScheme,
          Axis::ENUMDATA axis >
__attribute__((flatten))
void InteriorAdvectionCoefficientsDirection( FVCoefficients &fvCoeffs, 
                                             const FieldData<Tensor3D> &fields,
                                             const EnumVector<Axis, Tensor3D> &faceFluxes, 
                                             const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    constexpr bool hasDeferredCorrection = ( advectionScheme != AdvectionSchemes::Upwind );

    auto &AU     = fvCoeffs.Mom[X].AU;  // All equations share the same underlying data

    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    constexpr TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis], 
                                              west = LUT::LoCoeff[axis];

    floatType loCellLengthInv, hiCellLengthInv;
    TensorIndex3D loIndex, hiIndex;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {

        hiIndex[Z] = k;
        loIndex[Z] = k;

        if constexpr ( axis == Z ) {
            loCellLengthInv = mesh.cellLengthsInv[axis]( k-1 );
            hiCellLengthInv = mesh.cellLengthsInv[axis]( k   );
            loIndex[axis] -= 1;
        }

        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {

            hiIndex[Y] = j;
            loIndex[Y] = j;

            if constexpr ( axis == Y ) {
                loCellLengthInv = mesh.cellLengthsInv[axis]( j-1 );
                hiCellLengthInv = mesh.cellLengthsInv[axis]( j   );
                loIndex[axis] -= 1;
            }

            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                hiIndex[X] = i;
                loIndex[X] = i;

                if constexpr ( axis == X ) {
                    loCellLengthInv = mesh.cellLengthsInv[axis]( i-1 );
                    hiCellLengthInv = mesh.cellLengthsInv[axis]( i   );
                    loIndex[axis] -= 1;
                }

                loCellLengthInv = mesh.cellLengthsInv[axis]( loIndex[axis] );
                hiCellLengthInv = mesh.cellLengthsInv[axis]( hiIndex[axis] );

                const floatType faceFlux = faceFluxes[ axis ](i, j, k);
                const bool fluxIsPositive = ( faceFlux >= 0.0f );

                if ( fluxIsPositive ) {
                    AU[p   ](G(loIndex)) +=   faceFlux * loCellLengthInv;
                    AU[east](G(loIndex))  =   0.0f;
                    AU[p   ](G(hiIndex)) +=   0.0f;
                    AU[west](G(hiIndex))  = - faceFlux * hiCellLengthInv;
                } else {
                    AU[p   ](G(loIndex)) +=   0.0f;
                    AU[east](G(loIndex))  =   faceFlux * loCellLengthInv;
                    AU[p   ](G(hiIndex)) += - faceFlux * hiCellLengthInv;
                    AU[west](G(hiIndex))  =   0.0f;
                }


                if constexpr ( !hasDeferredCorrection )
                    continue;

                // Source term is different for each momentum equation (due to advected velocity)
                EnumFor<Axis>( [&] (Axis::ENUMDATA comp) {

                    const auto &U    = fields.U[comp];
                    auto &sourceTerm = fvCoeffs.Mom[comp].B;

                    floatType highOrderAdvectedVelocity{0.0f},
                              upwindAdvectedVelocity{0.0f};

                    if ( fluxIsPositive ) {
                        highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme         , +1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                        upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, +1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                    } else {
                        highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme         , -1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                        upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, -1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                    }

                    // Deferred correction term
                    sourceTerm(G(loIndex)) +=   fvCoeffs.advectionBlendingFactor 
                                            *   faceFlux
                                            *   ( highOrderAdvectedVelocity - upwindAdvectedVelocity ) 
                                            *   loCellLengthInv;
                        
                    sourceTerm(G(hiIndex)) += - fvCoeffs.advectionBlendingFactor 
                                            *   faceFlux
                                            *   ( highOrderAdvectedVelocity - upwindAdvectedVelocity ) 
                                            *   hiCellLengthInv;


                } );

            }
        }
    }

}



template< Axis::ENUMDATA axis>
struct AdvectionExecutionPolicy;
template<>
struct AdvectionExecutionPolicy< Axis::X > {
    using type = RAJA::KernelPolicy< RAJA::statement::Collapse<RAJA::omp_parallel_collapse_exec,
                                         RAJA::ArgList<2, 1>,
                                         RAJA::statement::For< 0, RAJA::seq_exec, 
                                             RAJA::statement::Lambda<0>
                                         >  
                                         
                                     > 
                                   >;  
};

template<>
struct AdvectionExecutionPolicy< Axis::Y > {
    using type = RAJA::KernelPolicy< RAJA::statement::For< 2, RAJA::omp_parallel_for_exec,
                                         RAJA::statement::For< 1, RAJA::seq_exec, 
                                             RAJA::statement::For< 0, RAJA::seq_exec, 
                                                RAJA::statement::Lambda<0>
                                             > 
                                         > 
                                     > 
                                   >;  
};

template<>
struct AdvectionExecutionPolicy< Axis::Z > {
    using type = RAJA::KernelPolicy< RAJA::statement::For< 2, RAJA::seq_exec,
                                         RAJA::statement::Collapse<RAJA::omp_parallel_collapse_exec, 
                                             RAJA::ArgList<1, 0>,
                                                RAJA::statement::Lambda<0>
                                         > 
                                     > 
                                   >;  
};


// Upwind implicit coefficients with deffered correction for higher order schemes
// Assumes that the 'p' coefficient has been set to zero
template< AdvectionSchemes advectionScheme,
          Axis::ENUMDATA axis >
__attribute__((flatten))
void InteriorAdvectionCoefficientsParallelDirection( FVCoefficients &fvCoeffs, 
                                                     const FieldData<Tensor3D> &fields,
                                                     const EnumVector<Axis, Tensor3D> &faceFluxes, 
                                                     const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    using execPolicy = AdvectionExecutionPolicy<axis>::type;  

    constexpr bool hasDeferredCorrection = ( advectionScheme != AdvectionSchemes::Upwind );

    auto &AU     = fvCoeffs.Mom[X].AU;  // All equations share the same underlying data
    
    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    constexpr TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis], 
                                              west = LUT::LoCoeff[axis];

    RAJA::kernel<execPolicy>( 
        RAJA::make_tuple( RAJA::TypedRangeSegment<intType>( startIndex[X], nFaces[X] ),
                          RAJA::TypedRangeSegment<intType>( startIndex[Y], nFaces[Y] ),
                          RAJA::TypedRangeSegment<intType>( startIndex[Z], nFaces[Z] ) 
                        ),
        [&] ( intType i, intType j, intType k ) {

            TensorIndex3D hiIndex = { i, j, k },
                          loIndex = { i, j, k };
            loIndex[axis] -= 1;


            const TensorIndex3D hiIndexG = G(hiIndex),
                                loIndexG = G(loIndex);

            const floatType faceFlux = faceFluxes[ axis ](i, j, k);
            const bool fluxIsPositive = ( faceFlux >= 0.0f );

            if ( fluxIsPositive ) {
                AU[p   ](loIndexG) +=   faceFlux * mesh.cellLengthsInv[axis]( loIndex[axis] );
                AU[east](loIndexG)  =   0.0f;
                AU[p   ](hiIndexG) +=   0.0f;
                AU[west](hiIndexG)  = - faceFlux * mesh.cellLengthsInv[axis]( hiIndex[axis] );
            } else {
                AU[p   ](loIndexG) +=   0.0f;
                AU[east](loIndexG)  =   faceFlux * mesh.cellLengthsInv[axis]( loIndex[axis] );
                AU[p   ](hiIndexG) += - faceFlux * mesh.cellLengthsInv[axis]( hiIndex[axis] );
                AU[west](hiIndexG)  =   0.0f;
            }


            if constexpr ( hasDeferredCorrection ) {
            
                // Source term is different for each momentum equation (due to advected velocity)
                EnumFor<Axis>( [&] (Axis::ENUMDATA comp) {

                    const auto &U   = fields.U[comp];
                    auto &sourceTerm = fvCoeffs.Mom[comp].B;

                    floatType highOrderAdvectedVelocity{0.0f},
                                upwindAdvectedVelocity{0.0f};

                    if ( fluxIsPositive ) {
                        highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme         , +1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                        upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, +1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                    } else {
                        highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme         , -1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                        upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, -1>(U, fvCoeffs, mesh, axis, hiIndex, loIndex);
                    }

                    // Deferred correction term
                    sourceTerm(loIndexG) +=   fvCoeffs.advectionBlendingFactor 
                                          *   faceFlux
                                          *   ( highOrderAdvectedVelocity - upwindAdvectedVelocity ) 
                                          *   mesh.cellLengthsInv[axis]( loIndex[axis] );
                        
                    sourceTerm(hiIndexG) += - fvCoeffs.advectionBlendingFactor 
                                          *   faceFlux
                                          *   ( highOrderAdvectedVelocity - upwindAdvectedVelocity ) 
                                          *   mesh.cellLengthsInv[axis]( hiIndex[axis] );

                } );

            }
            
        }  
    ); 

}



void BoundaryAdvectionCoefficients( FVCoefficients &fvCoeffs, 
                                    const EnumVector<Axis, Tensor3D> &faceFluxes, 
                                    const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    auto &AU     = fvCoeffs.Mom[X].AU;  // All equations share the same underlying data

    EnumFor<Axis>( [&] (Axis::ENUMDATA normal) {

        Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[normal],
                       axis2 = LUT::HiOrthogonalAxis[normal];

        TransportCoefficients::ENUMDATA east = LUT::HiCoeff[normal],
                                        west = LUT::LoCoeff[normal];

        // Iterate plane
        for ( intType idx1 = 0; idx1 != mesh.nCells(axis1); idx1++ ) {
            for ( intType idx2 = 0; idx2 != mesh.nCells(axis2); idx2++ ) {
                
                TensorIndex3D loFaceIndex;
                loFaceIndex[normal] = 0;
                loFaceIndex[axis1]  = idx1;
                loFaceIndex[axis2]  = idx2;

                TensorIndex3D hiFaceIndex = loFaceIndex;
                hiFaceIndex[normal] = mesh.nFacesNormal[normal](normal)-1;

                TensorIndex3D loCellIndex = loFaceIndex,
                              hiCellIndex = hiFaceIndex;
                hiCellIndex[normal] -= 1;

                // Hi side boundary
                AU[p   ](G(hiCellIndex)) += faceFluxes[normal](hiFaceIndex) 
                                          * ( 1.0f - mesh.interpFactors[normal]( hiFaceIndex[normal] ) ) 
                                          * mesh.cellLengthsInv[normal]( hiCellIndex[normal] );

                AU[east](G(hiCellIndex))  = faceFluxes[normal](hiFaceIndex) 
                                          * mesh.interpFactors[normal]( hiFaceIndex[normal] ) 
                                          * mesh.cellLengthsInv[normal]( hiCellIndex[normal] );


                // Lo side boundary
                AU[p   ](G(loCellIndex)) += - faceFluxes[normal](loFaceIndex) 
                                          *   mesh.interpFactors[normal]( loFaceIndex[normal] ) 
                                          *   mesh.cellLengthsInv[normal]( loCellIndex[normal] );

                AU[west](G(loCellIndex))  = - faceFluxes[normal](loFaceIndex) 
                                          *   ( 1.0f - mesh.interpFactors[normal]( loFaceIndex[normal] ) ) 
                                          *   mesh.cellLengthsInv[normal]( loCellIndex[normal] );

            }
        }       

    } );

    
}

[[maybe_unused]]
void SetAdvectionCoefficients( FVCoefficients &fvCoeffs,
                               const FieldData<Tensor3D> &fields,
                               const EnumVector<Axis, Tensor3D> &faceFluxes,
                               const Mesh &mesh )
{
    using enum TransportCoefficients::ENUMDATA;
    using enum Axis::ENUMDATA;
    
    // Boundary terms
    BoundaryAdvectionCoefficients(fvCoeffs, faceFluxes, mesh);

    // Interior terms
    switch ( fvCoeffs.advectionScheme ) {

        case AdvectionSchemes::Upwind:
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::Upwind, X>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::Upwind, Y>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::Upwind, Z>(fvCoeffs, fields, faceFluxes, mesh);
            break;

        case AdvectionSchemes::Central:
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::Central, X>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::Central, Y>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::Central, Z>(fvCoeffs, fields, faceFluxes, mesh);
            break;

        case AdvectionSchemes::SOU:
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::SOU, X>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::SOU, Y>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::SOU, Z>(fvCoeffs, fields, faceFluxes, mesh);
            break;

        case AdvectionSchemes::QUICK:
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::QUICK, X>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::QUICK, Y>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsDirection<AdvectionSchemes::QUICK, Z>(fvCoeffs, fields, faceFluxes, mesh);
            break;
    }      
}



void SetAdvectionCoefficientsParallel( FVCoefficients &fvCoeffs,
                                       const FieldData<Tensor3D> &fields,
                                       const EnumVector<Axis, Tensor3D> &faceFluxes,
                                       const Mesh &mesh )
{
    using enum TransportCoefficients::ENUMDATA;
    using enum Axis::ENUMDATA;
    
    // Boundary terms
    BoundaryAdvectionCoefficients(fvCoeffs, faceFluxes, mesh);

    // Interior terms
    switch ( fvCoeffs.advectionScheme ) {

        case AdvectionSchemes::Upwind:
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::Upwind, X>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::Upwind, Y>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::Upwind, Z>(fvCoeffs, fields, faceFluxes, mesh);
            break;

        case AdvectionSchemes::Central:
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::Central, X>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::Central, Y>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::Central, Z>(fvCoeffs, fields, faceFluxes, mesh);
            break;

        case AdvectionSchemes::SOU:
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::SOU, X>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::SOU, Y>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::SOU, Z>(fvCoeffs, fields, faceFluxes, mesh);
            break;

        case AdvectionSchemes::QUICK:
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::QUICK, X>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::QUICK, Y>(fvCoeffs, fields, faceFluxes, mesh);
            InteriorAdvectionCoefficientsParallelDirection<AdvectionSchemes::QUICK, Z>(fvCoeffs, fields, faceFluxes, mesh);
            break;
    }    
}


/*---------------------------------------------------------------------------------------------------------------*\
                                           Add Diffusion Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

[[maybe_unused]]
void AddDiffusion( FVCoefficients &fvCoeffs,
                   const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    auto &velCoeffs   = fvCoeffs.Mom[X].AU; // Shared data by all momentum equations
    const auto &diffCoeffs = fvCoeffs.diff;

    for (intType k = 0; k != mesh.nCells(Z); k++) {

        const intType kg = G(k);

        floatType zpk = diffCoeffs[Z][p](k),
                  ztk = diffCoeffs[Z][t](k),
                  zbk = diffCoeffs[Z][b](k);
        
        for (intType j = 0; j != mesh.nCells(Y); j++) {

            const intType jg = G(j);

            const floatType ypj = diffCoeffs[Y][p](j),
                            ynj = diffCoeffs[Y][n](j),
                            ysj = diffCoeffs[Y][s](j);

            CFD_PRAGMA_VECTORIZE
            for (intType i = 0; i != mesh.nCells(X); i++) {

                const intType ig = G(i);

                velCoeffs[p](ig, jg, kg) += diffCoeffs[X][p](i) + ypj + zpk;

                velCoeffs[e](ig, jg, kg) += diffCoeffs[X][e](i);
                velCoeffs[w](ig, jg, kg) += diffCoeffs[X][w](i);
                
                velCoeffs[n](ig, jg, kg) += ynj;
                velCoeffs[s](ig, jg, kg) += ysj;

                velCoeffs[t](ig, jg, kg) += ztk;
                velCoeffs[b](ig, jg, kg) += zbk;
            }
        }
    }

}




void AddDiffusionParallel( FVCoefficients &fvCoeffs,
                           const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    auto &velCoeffs   = fvCoeffs.Mom[X].AU; // Shared data by all momentum equations
    const auto &diffCoeffs = fvCoeffs.diff;

    using execPolicy = RAJA::KernelPolicy< RAJA::statement::Collapse<RAJA::omp_parallel_collapse_exec,
                                           RAJA::ArgList<2, 1, 0>,
                                               RAJA::statement::Lambda<0>
                                            > 
                                         >;  

    RAJA::kernel<execPolicy>( 
        RAJA::make_tuple( RAJA::TypedRangeSegment<intType>( 0, mesh.nCells(X) ),
                          RAJA::TypedRangeSegment<intType>( 0, mesh.nCells(Y) ),
                          RAJA::TypedRangeSegment<intType>( 0, mesh.nCells(Z) ) 
                        ),
                            [&] ( intType i, intType j, intType k ) {
                            const intType ig = G(i),
                                          jg = G(j),
                                          kg = G(k);

                            velCoeffs[p](ig, jg, kg) += diffCoeffs[X][p](i) + diffCoeffs[Y][p](j) + diffCoeffs[Z][p](k);

                            velCoeffs[e](ig, jg, kg) += diffCoeffs[X][e](i);
                            velCoeffs[w](ig, jg, kg) += diffCoeffs[X][w](i);
                            
                            velCoeffs[n](ig, jg, kg) += diffCoeffs[Y][n](j);
                            velCoeffs[s](ig, jg, kg) += diffCoeffs[Y][s](j);

                            velCoeffs[t](ig, jg, kg) += diffCoeffs[Z][t](k);
                            velCoeffs[b](ig, jg, kg) += diffCoeffs[Z][b](k);
                            }
                        );

}



/*---------------------------------------------------------------------------------------------------------------*\
                                        Linear Interpolated Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


void SetCellGradientCoefficients( FVCoefficients &fvCoeffs,
                                  const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        auto &coeffs = fvCoeffs.Cont.AU[axis];      // Shared by momentum pressure coeffs

        const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis],    // These are just names, they can be north, south etc.
                                              west = LUT::LoCoeff[axis];  

        // Internal faces
        for (intType i = 1; i != mesh.nCells(axis); i++) {
            
            // Cell on west side
            coeffs[p   ](G(i-1)) += 1 - mesh.interpFactors[axis](i);
            coeffs[east](G(i-1)) += mesh.interpFactors[axis](i);

            // Cell on east side
            coeffs[p   ](G(i)) += - mesh.interpFactors[axis](i);
            coeffs[west](G(i)) += - ( 1 - mesh.interpFactors[axis](i) ); 

        }

        // Lo boundary
        coeffs[p   ](G(0)) += - mesh.interpFactors[axis](0);
        coeffs[west](G(0)) += - ( 1 - mesh.interpFactors[axis](0) ); 

        // Hi boundary
        coeffs[p   ](G(mesh.nCells(axis)-1)) += 1 - mesh.interpFactors[axis]( mesh.nFacesNormal[axis](axis)-1 );
        coeffs[east](G(mesh.nCells(axis)-1)) += mesh.interpFactors[axis]( mesh.nFacesNormal[axis](axis)-1 );

        
        // Multiply by inverse of cell length
        for (intType i = 0; i != mesh.nCells(axis); i++) {
            coeffs[p   ](G(i)) *= mesh.cellLengthsInv[axis](i);
            coeffs[east](G(i)) *= mesh.cellLengthsInv[axis](i);
            coeffs[west](G(i)) *= mesh.cellLengthsInv[axis](i); 
        }

    } );

}



/*---------------------------------------------------------------------------------------------------------------*\
                        Momentum Weighted Interpolation (Rhie-Chow Interpolation) Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

// Unweighted constants that appear in the sparse gradient part of MWI
void SetMomentumInterpolationSparseConstants( FVCoefficients &fvCoeffs,
                                              const Mesh &mesh)
{
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        auto &mwiSparseCoeffs = fvCoeffs.mwiSparseCoeffs[axis];
        const auto &momentumPressureCoeffs = fvCoeffs.Mom[axis].AP;

        const TransportCoefficients::ENUMDATA west = LUT::LoCoeff[axis],
                                              east = LUT::HiCoeff[axis];


        // Internal faces
        for ( intType i = 1; i != mesh.cellFaces[axis].size()-1; i++ ) {

            mwiSparseCoeffs[0](i) = (1 - mesh.interpFactors[axis](i))   * momentumPressureCoeffs[west](G(i-1));

            mwiSparseCoeffs[1](i) = ( (1 - mesh.interpFactors[axis](i)) * momentumPressureCoeffs[p   ](G(i-1))
                                  +    mesh.interpFactors[axis](i)      * momentumPressureCoeffs[west](G(i  )) );

            mwiSparseCoeffs[2](i) = ( (1 - mesh.interpFactors[axis](i)) * momentumPressureCoeffs[east](G(i-1))
                                  +   mesh.interpFactors[axis](i)       * momentumPressureCoeffs[p   ](G(i  )) );

            mwiSparseCoeffs[3](i) = mesh.interpFactors[axis](i) * momentumPressureCoeffs[east](G(i));
        }

    } );
}



// Unweighted constants that appear in the compact gradient part of MWI
void SetMomentumInterpolationCompactConstants( FVCoefficients &fvCoeffs,
                                               const Mesh &mesh )
{
    using enum TransportCoefficients::ENUMDATA;

    const floatType rhoInv = 1 / fvCoeffs.rho;

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        auto &mwiCompactCoeffs = fvCoeffs.mwiCompactCoeffs[axis];

        // Internal faces
        for ( intType i = 1; i != mesh.cellFaces[axis].size()-1; i++ ) {

            mwiCompactCoeffs[0](i) =   mesh.cellCenterDiffInv[axis](i) * rhoInv;

            mwiCompactCoeffs[1](i) = - mesh.cellCenterDiffInv[axis](i) * rhoInv;
        }

    } );
}



// Cell weighting coefficient for MWI.
__attribute__((always_inline)) 
inline floatType MWIWeightingCoeff( const TensorIndex3D &loIndex,
                                    const TensorIndex3D &hiIndex,
                                    const Tensor3D &AUp, 
                                    const Mesh& mesh,
                                    const Axis::ENUMDATA axis)
{
    const intType idx = hiIndex[axis];    // Axis index of the face
    const floatType interpFactor = mesh.interpFactors[axis]( idx );
    const TensorIndex3D loIndexG = G(loIndex),  // This is faster than using 'G' inline for some reason
                        hiIndexG = G(hiIndex);
    return  ( 1.0f - interpFactor ) * ( 1.0f / AUp( loIndexG ) )  
         +  interpFactor            * ( 1.0f / AUp( hiIndexG ) );
}



// Fully implicit momentum interpolation coefficient for internal faces
[[maybe_unused]]
__attribute__((flatten))
void MWInterpolationInteriorImplicit( FVCoefficients &fvCoeffs,
                                      const Mesh &mesh,
                                      const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    // Unpack
    EnumVector<TransportCoefficients, Tensor3D> &continuityPressureCoeffs = fvCoeffs.Cont.AP;
    const Tensor3D &momentumDiagCoeff                                     = fvCoeffs.Mom[X].AU[p];
    const std::array< Tensor1D, 4 > &mwiSparseCoeffs                      = fvCoeffs.mwiSparseCoeffs[axis];
    const std::array< Tensor1D, 2 > &mwiCompactCoeffs                     = fvCoeffs.mwiCompactCoeffs[axis];

    // Cell indexing
    const TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                          eeast = LUT::HiHiCoeff[axis],
                                          west  = LUT::LoCoeff[axis],
                                          wwest = LUT::LoLoCoeff[axis];

    // Set the first most plane to zero only. We do this so that the high coefficients can be set in place, and less coefficients have to 
    // be zeroed upon re-linearisation.
    continuityPressureCoeffs[east ].chip(G(0), axis).setZero();
    continuityPressureCoeffs[west ].chip(G(0), axis).setZero();

    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D hiIndex = { i, j, k },
                              loIndex = { i, j, k };
                loIndex[axis] -= 1;

                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeff, mesh, axis );

                // Coefficients for westmost to eastmost cell
                const intType idx = hiIndex[axis];
                floatType coeff0 = d * mwiSparseCoeffs[0](idx),
                          coeff1 = d * ( mwiSparseCoeffs[1](idx) + mwiCompactCoeffs[0](idx) ),
                          coeff2 = d * ( mwiSparseCoeffs[2](idx) + mwiCompactCoeffs[1](idx) ),
                          coeff3 = d * mwiSparseCoeffs[3](idx);

                // Cell on west side 
                const floatType LoCellLengthInv = mesh.cellLengthsInv[axis]( loIndex[axis] );
                loIndex = G(loIndex);
                continuityPressureCoeffs[west ](loIndex) += coeff0 * LoCellLengthInv;
                continuityPressureCoeffs[p    ](loIndex) += coeff1 * LoCellLengthInv;
                continuityPressureCoeffs[east ](loIndex) += coeff2 * LoCellLengthInv;
                continuityPressureCoeffs[eeast](loIndex)  = coeff3 * LoCellLengthInv;

                // Cell on east side
                const floatType HiCellLengthInv = mesh.cellLengthsInv[axis]( hiIndex[axis] );
                hiIndex = G(hiIndex);
                continuityPressureCoeffs[wwest](hiIndex)  = - coeff0 * HiCellLengthInv;
                continuityPressureCoeffs[west ](hiIndex)  = - coeff1 * HiCellLengthInv;
                continuityPressureCoeffs[p    ](hiIndex) += - coeff2 * HiCellLengthInv;
                continuityPressureCoeffs[east ](hiIndex)  = - coeff3 * HiCellLengthInv;

            }
        }
    }

}



// Fully implicit momentum interpolation coefficient for internal faces
template<Axis::ENUMDATA axis>
[[maybe_unused]] 
__attribute__((flatten))
void MWInterpolationInteriorImplicit( FVCoefficients &fvCoeffs,
                                      const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    // Unpack
    EnumVector<TransportCoefficients, Tensor3D> &continuityPressureCoeffs = fvCoeffs.Cont.AP;
    const Tensor3D &momentumDiagCoeff                                     = fvCoeffs.Mom[X].AU[p];
    const std::array< Tensor1D, 4 > &mwiSparseCoeffs                      = fvCoeffs.mwiSparseCoeffs[axis];
    const std::array< Tensor1D, 2 > &mwiCompactCoeffs                     = fvCoeffs.mwiCompactCoeffs[axis];

    // Cell indexing
    constexpr TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                              eeast = LUT::HiHiCoeff[axis],
                                              west  = LUT::LoCoeff[axis],
                                              wwest = LUT::LoLoCoeff[axis];

    // Set the first most plane to zero only. We do this so that the high coefficients can be set in place, and less coefficients have to 
    // be zeroed upon re-linearisation.
    continuityPressureCoeffs[east ].chip(G(0), axis).setZero();
    continuityPressureCoeffs[west ].chip(G(0), axis).setZero();

    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    TensorIndex3D loIndex, hiIndex;
    intType idx;
    floatType loCellLengthInv, hiCellLengthInv;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {

        hiIndex[Z] = k;
        loIndex[Z] = k;

        if constexpr ( axis == Z ) {
            loIndex[axis] -= 1;
            idx = k;
            loCellLengthInv = mesh.cellLengthsInv[axis]( k-1 );
            hiCellLengthInv = mesh.cellLengthsInv[axis]( k   );
        }

        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {

            hiIndex[Y] = j;
            loIndex[Y] = j;

            if constexpr ( axis == Y ) {
                loIndex[axis] -= 1;
                idx = j;
                loCellLengthInv = mesh.cellLengthsInv[axis]( j-1 );
                hiCellLengthInv = mesh.cellLengthsInv[axis]( j   );
            }

            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                hiIndex[X] = i;
                loIndex[X] = i;
                
                if constexpr ( axis == X ) {
                    loIndex[axis] -= 1;
                    idx = i;
                    loCellLengthInv = mesh.cellLengthsInv[axis]( i-1 );
                    hiCellLengthInv = mesh.cellLengthsInv[axis]( i   );
                }

                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeff, mesh, axis );

                // Coefficients for westmost to eastmost cell
                floatType coeff0 = d *   mwiSparseCoeffs[0](idx),
                          coeff1 = d * ( mwiSparseCoeffs[1](idx) + mwiCompactCoeffs[0](idx) ),
                          coeff2 = d * ( mwiSparseCoeffs[2](idx) + mwiCompactCoeffs[1](idx) ),
                          coeff3 = d *   mwiSparseCoeffs[3](idx);

                // Cell on west side 
                const TensorIndex3D loIndexG = G(loIndex);
                continuityPressureCoeffs[west ](loIndexG) += coeff0 * loCellLengthInv;
                continuityPressureCoeffs[p    ](loIndexG) += coeff1 * loCellLengthInv;
                continuityPressureCoeffs[east ](loIndexG) += coeff2 * loCellLengthInv;
                continuityPressureCoeffs[eeast](loIndexG)  = coeff3 * loCellLengthInv;

                // Cell on east side
                const TensorIndex3D hiIndexG = G(hiIndex);
                continuityPressureCoeffs[wwest](hiIndexG)  = - coeff0 * hiCellLengthInv;
                continuityPressureCoeffs[west ](hiIndexG)  = - coeff1 * hiCellLengthInv;
                continuityPressureCoeffs[p    ](hiIndexG) += - coeff2 * hiCellLengthInv;
                continuityPressureCoeffs[east ](hiIndexG)  = - coeff3 * hiCellLengthInv;

            }
        }
    }

}



template< Axis::ENUMDATA axis>
struct MWIExecutionPolicy;
template<>
struct MWIExecutionPolicy< Axis::X > {
    using type = RAJA::KernelPolicy< RAJA::statement::Collapse<RAJA::omp_parallel_collapse_exec,
                                         RAJA::ArgList<2, 1>,
                                         RAJA::statement::For< 0, RAJA::seq_exec, 
                                             RAJA::statement::Lambda<0>
                                         >  
                                         
                                     > 
                                   >;  
};

template<>
struct MWIExecutionPolicy< Axis::Y > {
    using type = RAJA::KernelPolicy< RAJA::statement::For< 2, RAJA::omp_parallel_for_exec,
                                         RAJA::statement::For< 1, RAJA::seq_exec, 
                                             RAJA::statement::For< 0, RAJA::seq_exec, 
                                                RAJA::statement::Lambda<0>
                                             > 
                                         > 
                                     > 
                                   >;  
};

template<>
struct MWIExecutionPolicy< Axis::Z > {
    using type = RAJA::KernelPolicy< RAJA::statement::For< 2, RAJA::seq_exec,
                                         RAJA::statement::Collapse<RAJA::omp_parallel_collapse_exec, 
                                             RAJA::ArgList<1, 0>,
                                                RAJA::statement::Lambda<0>
                                         > 
                                     > 
                                   >;  
};


// Fully implicit momentum interpolation coefficient for internal faces
template< Axis::ENUMDATA axis >
[[maybe_unused]]
__attribute__((flatten))
void MWInterpolationInteriorImplicitParallel( FVCoefficients &fvCoeffs,
                                              const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    using execPolicy = MWIExecutionPolicy<axis>::type;  

    // Unpack
    EnumVector<TransportCoefficients, Tensor3D> &continuityPressureCoeffs = fvCoeffs.Cont.AP;
    const Tensor3D &momentumDiagCoeff                                     = fvCoeffs.Mom[X].AU[p];
    const std::array< Tensor1D, 4 > &mwiSparseCoeffs                      = fvCoeffs.mwiSparseCoeffs[axis];
    const std::array< Tensor1D, 2 > &mwiCompactCoeffs                     = fvCoeffs.mwiCompactCoeffs[axis];

    // Cell indexing
    constexpr TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                              eeast = LUT::HiHiCoeff[axis],
                                              west  = LUT::LoCoeff[axis],
                                              wwest = LUT::LoLoCoeff[axis];

    // Set the first most plane to zero only. We do this so that the high coefficients can be set in place, and less coefficients have to 
    // be zeroed upon re-linearisation.
    continuityPressureCoeffs[east ].chip(G(0), axis).setZero();
    continuityPressureCoeffs[west ].chip(G(0), axis).setZero();

    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    RAJA::kernel<execPolicy>( 
        RAJA::make_tuple( RAJA::TypedRangeSegment<intType>( startIndex[X], nFaces[X] ),
                          RAJA::TypedRangeSegment<intType>( startIndex[Y], nFaces[Y] ),
                          RAJA::TypedRangeSegment<intType>( startIndex[Z], nFaces[Z] ) 
                        ),
        [&] ( intType i, intType j, intType k ) {

                TensorIndex3D hiIndex = { i, j, k },
                              loIndex = { i, j, k };
                loIndex[axis] -= 1;

                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeff, mesh, axis );

                // Coefficients for westmost to eastmost cell
                const intType idx = hiIndex[axis];
                floatType coeff0 = d * mwiSparseCoeffs[0](idx),
                          coeff1 = d * ( mwiSparseCoeffs[1](idx) + mwiCompactCoeffs[0](idx) ),
                          coeff2 = d * ( mwiSparseCoeffs[2](idx) + mwiCompactCoeffs[1](idx) ),
                          coeff3 = d * mwiSparseCoeffs[3](idx);

                // Cell on west side 
                const floatType LoCellLengthInv = mesh.cellLengthsInv[axis]( loIndex[axis] );
                loIndex = G(loIndex);
                continuityPressureCoeffs[west ](loIndex) += coeff0 * LoCellLengthInv;
                continuityPressureCoeffs[p    ](loIndex) += coeff1 * LoCellLengthInv;
                continuityPressureCoeffs[east ](loIndex) += coeff2 * LoCellLengthInv;
                continuityPressureCoeffs[eeast](loIndex)  = coeff3 * LoCellLengthInv;

                // Cell on east side
                const floatType HiCellLengthInv = mesh.cellLengthsInv[axis]( hiIndex[axis] );
                hiIndex = G(hiIndex);
                continuityPressureCoeffs[wwest](hiIndex)  = - coeff0 * HiCellLengthInv;
                continuityPressureCoeffs[west ](hiIndex)  = - coeff1 * HiCellLengthInv;
                continuityPressureCoeffs[p    ](hiIndex) += - coeff2 * HiCellLengthInv;
                continuityPressureCoeffs[east ](hiIndex)  = - coeff3 * HiCellLengthInv;

        }
    );

}



// Semi explicit momentum interpolation coefficient for internal faces
[[maybe_unused]]
__attribute__((flatten))
void MWInterpolationInteriorSemiExplicit( FVCoefficients &fvCoeffs, 
                                          const Tensor3D &P,
                                          const Mesh &mesh,
                                          const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    // Unpack
    EnumVector<TransportCoefficients, Tensor3D> &continuityPressureCoeffs = fvCoeffs.Cont.AP;
    Tensor3D &continuitySourceTerm                                        = fvCoeffs.Cont.B;
    const Tensor3D &momentumDiagCoeff                                     = fvCoeffs.Mom[X].AU[p];
    const std::array< Tensor1D, 4 > &mwiSparseCoeffs                      = fvCoeffs.mwiSparseCoeffs[axis];
    const std::array< Tensor1D, 2 > &mwiCompactCoeffs                     = fvCoeffs.mwiCompactCoeffs[axis];

    // For getting the index of a neighbouring cell
    auto NeighbourIndex = [] ( TensorIndex3D index, intType shift, Axis::ENUMDATA shiftAxis ) { 
        index[shiftAxis] += shift; 
        return index; 
    };

    // Cell indexing
    const TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                          west  = LUT::LoCoeff[axis];

    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D hiIndex = { i, j, k },
                              loIndex = { i, j, k };
                loIndex[axis] -= 1;
                 
                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeff, mesh, axis );

                const floatType LoCellLengthInv = mesh.cellLengthsInv[axis]( loIndex[axis] ),
                          HiCellLengthInv = mesh.cellLengthsInv[axis]( hiIndex[axis] );


                // Implicit compact difference --------------------------------------------------------------------------
                                

                // Coefficients for westmost to eastmost cell
                const intType idx = hiIndex[axis];
                floatType coeffCompact0 = d * mwiCompactCoeffs[0](idx),
                          coeffCompact1 = d * mwiCompactCoeffs[1](idx);

                // Cell on west side 
                continuityPressureCoeffs[p    ](G(loIndex)) +=   coeffCompact0 * LoCellLengthInv;
                continuityPressureCoeffs[east ](G(loIndex))  =   coeffCompact1 * LoCellLengthInv;

                // Cell on east side
                continuityPressureCoeffs[west ](G(hiIndex))  = - coeffCompact0 * HiCellLengthInv;
                continuityPressureCoeffs[p    ](G(hiIndex)) += - coeffCompact1 * HiCellLengthInv;


                // Explicit sparse difference ---------------------------------------------------------------------------

                const TensorIndex3D loWest  = NeighbourIndex( loIndex, -1, axis ),
                                    loEast  = NeighbourIndex( loIndex,  1, axis ),
                                    loEEast = NeighbourIndex( loIndex,  2, axis ),

                                    hiWWest = NeighbourIndex( hiIndex, -2, axis ),
                                    hiWest  = NeighbourIndex( hiIndex, -1, axis ),
                                    hiEast  = NeighbourIndex( hiIndex,  1, axis );

                const floatType coeffSparse0 = d * mwiSparseCoeffs[0](idx),
                                coeffSparse1 = d * mwiSparseCoeffs[1](idx),
                                coeffSparse2 = d * mwiSparseCoeffs[2](idx),
                                coeffSparse3 = d * mwiSparseCoeffs[3](idx);

                // Cell on west side 
                continuitySourceTerm(G(loIndex)) += ( coeffSparse0 * P( G(loWest)  )
                                                    + coeffSparse1 * P( G(loIndex) )
                                                    + coeffSparse2 * P( G(loEast)  )
                                                    + coeffSparse3 * P( G(loEEast) )
                                                    ) * LoCellLengthInv;

                // Cell on east side
                continuitySourceTerm(G(hiIndex)) -= ( coeffSparse0 * P( G(hiWWest) )
                                                    + coeffSparse1 * P( G(hiWest)  )
                                                    + coeffSparse2 * P( G(hiIndex) )
                                                    + coeffSparse3 * P( G(hiEast)  )
                                                    ) * HiCellLengthInv;


                // ------------------------------------------------------------------------------------------------------

            }
        }
    }

}



// Semi explicit momentum interpolation coefficient for internal faces
template<Axis::ENUMDATA axis>
[[maybe_unused]] 
__attribute__((flatten))
void MWInterpolationInteriorSemiExplicit( FVCoefficients &fvCoeffs, 
                                          const Tensor3D &P,
                                          const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    // Unpack
    EnumVector<TransportCoefficients, Tensor3D> &continuityPressureCoeffs = fvCoeffs.Cont.AP;
    Tensor3D &continuitySourceTerm                                        = fvCoeffs.Cont.B;
    const Tensor3D &momentumDiagCoeff                                     = fvCoeffs.Mom[X].AU[p];
    const std::array< Tensor1D, 4 > &mwiSparseCoeffs                      = fvCoeffs.mwiSparseCoeffs[axis];
    const std::array< Tensor1D, 2 > &mwiCompactCoeffs                     = fvCoeffs.mwiCompactCoeffs[axis];

    // For getting the index of a neighbouring cell
    auto NeighbourIndex = [] ( TensorIndex3D index, intType shift, Axis::ENUMDATA shiftAxis ) { 
        index[shiftAxis] += shift; 
        return index; 
    };

    // Cell indexing
    constexpr TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                              west  = LUT::LoCoeff[axis];

    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    TensorIndex3D loIndex, hiIndex;
    intType idx;
    floatType loCellLengthInv, hiCellLengthInv;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {

        hiIndex[Z] = k;
        loIndex[Z] = k;

        if constexpr ( axis == Z ) {
            loIndex[axis] -= 1;
            idx = k;
            loCellLengthInv = mesh.cellLengthsInv[axis]( k-1 );
            hiCellLengthInv = mesh.cellLengthsInv[axis]( k   );
        }

        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {

            hiIndex[Y] = j;
            loIndex[Y] = j;

            if constexpr ( axis == Y ) {
                loIndex[axis] -= 1;
                idx = j;
                loCellLengthInv = mesh.cellLengthsInv[axis]( j-1 );
                hiCellLengthInv = mesh.cellLengthsInv[axis]( j   );
            }

            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                hiIndex[X] = i;
                loIndex[X] = i;
                
                if constexpr ( axis == X ) {
                    loIndex[axis] -= 1;
                    idx = i;
                    loCellLengthInv = mesh.cellLengthsInv[axis]( i-1 );
                    hiCellLengthInv = mesh.cellLengthsInv[axis]( i   );
                }

                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeff, mesh, axis );

                // Implicit compact difference --------------------------------------------------------------------------
                                
                // Coefficients for westmost to eastmost cell
                floatType coeffCompact0 = d * mwiCompactCoeffs[0](idx),
                          coeffCompact1 = d * mwiCompactCoeffs[1](idx);

                // Cell on west side 
                continuityPressureCoeffs[p    ](G(loIndex)) +=   coeffCompact0 * loCellLengthInv;
                continuityPressureCoeffs[east ](G(loIndex))  =   coeffCompact1 * loCellLengthInv;

                // Cell on east side
                continuityPressureCoeffs[west ](G(hiIndex))  = - coeffCompact0 * hiCellLengthInv;
                continuityPressureCoeffs[p    ](G(hiIndex)) += - coeffCompact1 * hiCellLengthInv;


                // Explicit sparse difference ---------------------------------------------------------------------------

                const TensorIndex3D loWest  = NeighbourIndex( loIndex, -1, axis ),
                                    loEast  = NeighbourIndex( loIndex,  1, axis ),
                                    loEEast = NeighbourIndex( loIndex,  2, axis ),

                                    hiWWest = NeighbourIndex( hiIndex, -2, axis ),
                                    hiWest  = NeighbourIndex( hiIndex, -1, axis ),
                                    hiEast  = NeighbourIndex( hiIndex,  1, axis );

                const floatType coeffSparse0 = d * mwiSparseCoeffs[0](idx),
                                coeffSparse1 = d * mwiSparseCoeffs[1](idx),
                                coeffSparse2 = d * mwiSparseCoeffs[2](idx),
                                coeffSparse3 = d * mwiSparseCoeffs[3](idx);

                // Cell on west side 
                continuitySourceTerm(G(loIndex)) += ( coeffSparse0 * P( G(loWest)  )
                                                    + coeffSparse1 * P( G(loIndex) )
                                                    + coeffSparse2 * P( G(loEast)  )
                                                    + coeffSparse3 * P( G(loEEast) )
                                                    ) * loCellLengthInv;

                // Cell on east side
                continuitySourceTerm(G(hiIndex)) -= ( coeffSparse0 * P( G(hiWWest) )
                                                    + coeffSparse1 * P( G(hiWest)  )
                                                    + coeffSparse2 * P( G(hiIndex) )
                                                    + coeffSparse3 * P( G(hiEast)  )
                                                    ) * hiCellLengthInv;


                // ------------------------------------------------------------------------------------------------------

            }
        }
    }

}



template< Axis::ENUMDATA axis >
[[maybe_unused]]
__attribute__((flatten))
void MWInterpolationInteriorSemiExplicitParallel( FVCoefficients &fvCoeffs, 
                                                  const Tensor3D &P,
                                                  const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    using execPolicy = AdvectionExecutionPolicy<axis>::type;  

    // Unpack
    EnumVector<TransportCoefficients, Tensor3D> &continuityPressureCoeffs = fvCoeffs.Cont.AP;
    Tensor3D &continuitySourceTerm                                        = fvCoeffs.Cont.B;
    const Tensor3D &momentumDiagCoeff                                     = fvCoeffs.Mom[X].AU[p];
    const std::array< Tensor1D, 4 > &mwiSparseCoeffs                      = fvCoeffs.mwiSparseCoeffs[axis];
    const std::array< Tensor1D, 2 > &mwiCompactCoeffs                     = fvCoeffs.mwiCompactCoeffs[axis];

    // For getting the index of a neighbouring cell
    auto NeighbourIndex = [] ( TensorIndex3D index, intType shift, Axis::ENUMDATA shiftAxis ) { 
        index[shiftAxis] += shift; 
        return index; 
    };

    // Cell indexing
    constexpr TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                              west  = LUT::LoCoeff[axis];

    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    RAJA::kernel<execPolicy>( 
        RAJA::make_tuple( RAJA::TypedRangeSegment<intType>( startIndex[X], nFaces[X] ),
                          RAJA::TypedRangeSegment<intType>( startIndex[Y], nFaces[Y] ),
                          RAJA::TypedRangeSegment<intType>( startIndex[Z], nFaces[Z] ) 
                        ),
        [&] ( intType i, intType j, intType k ) {

                TensorIndex3D hiIndex = { i, j, k },
                              loIndex = { i, j, k };
                loIndex[axis] -= 1;
                 
                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeff, mesh, axis );

                const floatType LoCellLengthInv = mesh.cellLengthsInv[axis]( loIndex[axis] ),
                                HiCellLengthInv = mesh.cellLengthsInv[axis]( hiIndex[axis] );


                // Implicit compact difference --------------------------------------------------------------------------
                                

                // Coefficients for westmost to eastmost cell
                const intType idx = hiIndex[axis];
                floatType coeffCompact0 = d * mwiCompactCoeffs[0](idx),
                          coeffCompact1 = d * mwiCompactCoeffs[1](idx);

                // Cell on west side 
                continuityPressureCoeffs[p    ](G(loIndex)) +=   coeffCompact0 * LoCellLengthInv;
                continuityPressureCoeffs[east ](G(loIndex))  =   coeffCompact1 * LoCellLengthInv;

                // Cell on east side
                continuityPressureCoeffs[west ](G(hiIndex))  = - coeffCompact0 * HiCellLengthInv;
                continuityPressureCoeffs[p    ](G(hiIndex)) += - coeffCompact1 * HiCellLengthInv;


                // Explicit sparse difference ---------------------------------------------------------------------------

                const TensorIndex3D loWest  = NeighbourIndex( loIndex, -1, axis ),
                                    loEast  = NeighbourIndex( loIndex,  1, axis ),
                                    loEEast = NeighbourIndex( loIndex,  2, axis ),

                                    hiWWest = NeighbourIndex( hiIndex, -2, axis ),
                                    hiWest  = NeighbourIndex( hiIndex, -1, axis ),
                                    hiEast  = NeighbourIndex( hiIndex,  1, axis );

                const floatType coeffSparse0 = d * mwiSparseCoeffs[0](idx),
                                coeffSparse1 = d * mwiSparseCoeffs[1](idx),
                                coeffSparse2 = d * mwiSparseCoeffs[2](idx),
                                coeffSparse3 = d * mwiSparseCoeffs[3](idx);

                // Cell on west side 
                continuitySourceTerm(G(loIndex)) += ( coeffSparse0 * P( G(loWest)  )
                                                    + coeffSparse1 * P( G(loIndex) )
                                                    + coeffSparse2 * P( G(loEast)  )
                                                    + coeffSparse3 * P( G(loEEast) )
                                                    ) * LoCellLengthInv;

                // Cell on east side
                continuitySourceTerm(G(hiIndex)) -= ( coeffSparse0 * P( G(hiWWest) )
                                                    + coeffSparse1 * P( G(hiWest)  )
                                                    + coeffSparse2 * P( G(hiIndex) )
                                                    + coeffSparse3 * P( G(hiEast)  )
                                                    ) * HiCellLengthInv;


                // ------------------------------------------------------------------------------------------------------

        }
    );

}



// Set momentum interpolation coefficients
[[maybe_unused]]
void SetMomentumInterpolationCoefficientsParallel( FVCoefficients &fvCoeffs,
                                                   const Mesh &mesh,
                                                   const Tensor3D &P )
{
    switch ( fvCoeffs.momentumInterpolation ) {
        case MomentumInterpolation::Implicit:
            MWInterpolationInteriorImplicitParallel<Axis::X>(fvCoeffs, mesh);
            MWInterpolationInteriorImplicitParallel<Axis::Y>(fvCoeffs, mesh);
            MWInterpolationInteriorImplicitParallel<Axis::Z>(fvCoeffs, mesh);
            break;

        case MomentumInterpolation::SemiExplicit:
            MWInterpolationInteriorSemiExplicitParallel<Axis::X>(fvCoeffs, P, mesh);
            MWInterpolationInteriorSemiExplicitParallel<Axis::Y>(fvCoeffs, P, mesh);
            MWInterpolationInteriorSemiExplicitParallel<Axis::Z>(fvCoeffs, P, mesh);
            break;
    }
}


// Set momentum interpolation coefficients
[[maybe_unused]]
void SetMomentumInterpolationCoefficients( FVCoefficients &fvCoeffs,
                                           const Mesh &mesh,
                                           const Tensor3D &P )
{
    switch ( fvCoeffs.momentumInterpolation ) {
        case MomentumInterpolation::Implicit:
            MWInterpolationInteriorImplicit<Axis::X>(fvCoeffs, mesh);
            MWInterpolationInteriorImplicit<Axis::Y>(fvCoeffs, mesh);
            MWInterpolationInteriorImplicit<Axis::Z>(fvCoeffs, mesh);
            break;

        case MomentumInterpolation::SemiExplicit:
            MWInterpolationInteriorSemiExplicit<Axis::X>(fvCoeffs, P, mesh);
            MWInterpolationInteriorSemiExplicit<Axis::Y>(fvCoeffs, P, mesh);
            MWInterpolationInteriorSemiExplicit<Axis::Z>(fvCoeffs, P, mesh);
            break;
    }
}


/*---------------------------------------------------------------------------------------------------------------*\
                                                Unsteady Terms
\*---------------------------------------------------------------------------------------------------------------*/


void AddBackwardsEuler( FVCoefficients &fvCoeffs, 
                        const Mesh &mesh,
                        const FieldData<Tensor3D> &fieldsPrevTime )
{
    using enum TransportCoefficients::ENUMDATA;
    const floatType dtInv = 1.0f / fvCoeffs.timeStep;

    for (intType k = 0; k != mesh.nCells(2); k++) {
        for (intType j = 0; j != mesh.nCells(1); j++) {
            for (intType i = 0; i != mesh.nCells(0); i++) {

                // Velocity coefficient is shared in momentum equations
                fvCoeffs.Mom[Axis::X].AU[p]( G(i, j, k) ) += dtInv;

                // Each equation has a different source term
                EnumFor<Axis>( [&] (Axis::ENUMDATA comp) {
                    fvCoeffs.Mom[comp].B( G(i, j, k) ) += - fieldsPrevTime.U[comp]( G(i, j, k) ) * dtInv;
                } );
                

            }
        }
    }
}



void AddBackwardsThreeLevel( FVCoefficients &fvCoeffs, 
                             const Mesh &mesh,
                             const FieldData<Tensor3D> &fieldsPrevTime,
                             const FieldData<Tensor3D> &fieldsPrevPrevTime )
{
    using enum TransportCoefficients::ENUMDATA;
    const floatType dtInv = 1.0f / fvCoeffs.timeStep;

    for (intType k = 0; k != mesh.nCells(2); k++) {
        for (intType j = 0; j != mesh.nCells(1); j++) {
            for (intType i = 0; i != mesh.nCells(0); i++) {

                // Velocity coefficient is shared in momentum equations
                fvCoeffs.Mom[Axis::X].AU[p]( G(i, j, k) ) += 1.5f * dtInv;

                // Each equation has a different source term
                EnumFor<Axis>( [&] (Axis::ENUMDATA comp) {
                    fvCoeffs.Mom[comp].B( G(i, j, k) )     += (
                                                                - 2.0f * fieldsPrevTime.U[comp]( G(i, j, k) )
                                                                + 0.5f * fieldsPrevPrevTime.U[comp]( G(i, j, k) )
                                                              ) * dtInv;
                } );

            }
        }
    }
}



void AddUnsteadyTerm( FVCoefficients &fvCoeffs,
                      const FieldData<Tensor3D> &fieldsPrevTime,
                      const FieldData<Tensor3D> &fieldsPrevPrevTime,
                      const Mesh &mesh )
{

    switch ( fvCoeffs.timeScheme ) {
        case TimeSchemes::BackwardsEuler:
            AddBackwardsEuler(fvCoeffs, mesh, fieldsPrevTime);
            break;

        case TimeSchemes::BackwardsThreeLevel:
            AddBackwardsThreeLevel(fvCoeffs, mesh, fieldsPrevTime, fieldsPrevPrevTime);
            break;

        case TimeSchemes::Steady:
            /* NULL */
            break;
    }

}



/*---------------------------------------------------------------------------------------------------------------*\
                                            Immersed Boundary Functions
\*---------------------------------------------------------------------------------------------------------------*/


// Set the IB source terms that come from the implicit stencil 
void MomentumIBSourceStencil( FVCoefficients &fvCoeffs,
                              const IBCell::SourceTermData &sourceTermData, 
                              const TensorIndex3D &cellIndex )
{
    using FVT::G;

    const Axis::ENUMDATA faceNormal = sourceTermData.direction;
    const TransportCoefficients::ENUMDATA coeff = ( sourceTermData.directionIndex == +1 ) ?  LUT::HiCoeff[faceNormal] : LUT::LoCoeff[faceNormal];

    EnumFor<Axis>( [&] (Axis::ENUMDATA momentumAxis) {

        // Velocity term
        floatType ibSource = fvCoeffs.Mom[momentumAxis].AU[coeff]( G(cellIndex) ) * sourceTermData.ghostCellValues.U[momentumAxis];

        // Pressure stencil
        if ( momentumAxis == faceNormal ) {
            ibSource += fvCoeffs.Mom[momentumAxis].AP[coeff]( G(cellIndex[faceNormal]) ) * sourceTermData.ghostCellValues.P;
        }

        fvCoeffs.Mom[momentumAxis].B( G(cellIndex) ) += ibSource;

    } );
}



template< AdvectionSchemes advectionScheme >
void MomentumIBSourceDeferredCorrection( FVCoefficients &fvCoeffs,
                                         const FieldData<Tensor3D> &fields,
                                         const EnumVector< Axis, Tensor3D > &faceFluxes,
                                         const Mesh &mesh,
                                         const IBCell::SourceTermData &sourceTermData, 
                                         const TensorIndex3D &cellIndex )
{
    using FVT::G; 
    
    const Axis::ENUMDATA faceNormal = sourceTermData.direction;
    TensorIndex3D faceIndex = cellIndex;
    faceIndex[faceNormal] += sourceTermData.faceDirectionIndex;

    EnumFor<Axis>( [&] (Axis::ENUMDATA momentumAxis) {

        // Advected velocities to remove
        floatType highOrderAdvectedVelocity{0.0f}, upwindAdvectedVelocity{0.0f};
        Tensor3D const &U = fields.U[ momentumAxis ];
        TensorIndex3D loIndex = faceIndex,
                    hiIndex = loIndex;
        hiIndex[faceNormal] += 1;
        if ( faceFluxes[faceNormal]( faceIndex ) >= 0.0f ) {
            highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme, +1>(U, fvCoeffs, mesh, faceNormal, hiIndex, loIndex);
            upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, +1>(U, fvCoeffs, mesh, faceNormal, hiIndex, loIndex);
        } else {
            highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme, -1>(U, fvCoeffs, mesh, faceNormal, hiIndex, loIndex);
            upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, -1>(U, fvCoeffs, mesh, faceNormal, hiIndex, loIndex);
        }

        // Remove the deferred correction term
        floatType ibSource = static_cast<floatType>( sourceTermData.directionIndex ) 
                           * fvCoeffs.advectionBlendingFactor 
                           * faceFluxes[faceNormal]( faceIndex )
                           * ( highOrderAdvectedVelocity - upwindAdvectedVelocity )
                           * mesh.cellLengthsInv[faceNormal]( cellIndex[faceNormal] );


        // Need to add effect of ghost cell for face one in from boundary
        constexpr bool hasWideAdvectionStencil = ( advectionScheme == AdvectionSchemes::SOU   ) ||
                                                 ( advectionScheme == AdvectionSchemes::QUICK );
        if constexpr ( hasWideAdvectionStencil ) {

            TensorIndex3D faceIndex_a = faceIndex;
            faceIndex_a[faceNormal] -= sourceTermData.directionIndex;

            TensorIndex3D cellIndex_a = cellIndex;
            cellIndex_a[faceNormal] -= sourceTermData.directionIndex;

            TensorIndex3D loIndex_a = faceIndex_a,
                        hiIndex_a = loIndex_a;
            hiIndex_a[faceNormal] += 1;

            const floatType ghostCellValue = sourceTermData.ghostCellValues.U[momentumAxis];

            floatType oldHighOrderAdvectedVelocity{0.0f}, correctedHighOrderAdvectedVelocity{0.0f};
            if ( faceFluxes[faceNormal]( faceIndex ) >= 0.0f ) {
                oldHighOrderAdvectedVelocity       = FaceInterpolatedVelocity<advectionScheme, +1>(U, fvCoeffs, mesh, faceNormal, hiIndex_a, loIndex_a);
                correctedHighOrderAdvectedVelocity = ( sourceTermData.directionIndex == +1 ) ? FaceInterpolatedVelocity<advectionScheme, +1, +1>(U, fvCoeffs, mesh, faceNormal, hiIndex_a, loIndex_a, ghostCellValue):
                                                                                               FaceInterpolatedVelocity<advectionScheme, +1, -1>(U, fvCoeffs, mesh, faceNormal, hiIndex_a, loIndex_a, ghostCellValue);
                
            } else {
                oldHighOrderAdvectedVelocity       = FaceInterpolatedVelocity<advectionScheme, -1>(U, fvCoeffs, mesh, faceNormal, hiIndex_a, loIndex_a);
                correctedHighOrderAdvectedVelocity = ( sourceTermData.directionIndex == +1 ) ? FaceInterpolatedVelocity<advectionScheme, -1, +1>(U, fvCoeffs, mesh, faceNormal, hiIndex_a, loIndex_a, ghostCellValue):
                                                                                               FaceInterpolatedVelocity<advectionScheme, -1, -1>(U, fvCoeffs, mesh, faceNormal, hiIndex_a, loIndex_a, ghostCellValue);
            }


            // Boundary cell
            ibSource -= - static_cast<floatType>( sourceTermData.directionIndex ) 
                    *   fvCoeffs.advectionBlendingFactor 
                    *   faceFluxes[faceNormal]( faceIndex_a )
                    *   ( correctedHighOrderAdvectedVelocity - oldHighOrderAdvectedVelocity )
                    *   mesh.cellLengthsInv[faceNormal]( cellIndex[faceNormal] );

            // Interior cell
            const floatType ibSource_a = - static_cast<floatType>( sourceTermData.directionIndex )
                                        *   fvCoeffs.advectionBlendingFactor 
                                        *   faceFluxes[faceNormal]( faceIndex_a )
                                        *   ( correctedHighOrderAdvectedVelocity - oldHighOrderAdvectedVelocity )
                                        *   mesh.cellLengthsInv[faceNormal]( cellIndex_a[faceNormal] );

            fvCoeffs.Mom[momentumAxis].B( G(cellIndex_a) ) += ibSource_a;

        }

        fvCoeffs.Mom[momentumAxis].B( G(cellIndex) ) += ibSource;

    } );
}



void ContinuityIBSourceImplicitMWI( FVCoefficients &fvCoeffs,
                                    const IBCell::SourceTermData &sourceTermData, 
                                    const TensorIndex3D &cellIndex ) 
{
    using FVT::G;

    const Axis::ENUMDATA faceNormal = sourceTermData.direction;
    const TransportCoefficients::ENUMDATA coeff  = ( sourceTermData.directionIndex == +1 ) ? LUT::HiCoeff[faceNormal]   : LUT::LoCoeff[faceNormal];
    const TransportCoefficients::ENUMDATA ccoeff = ( sourceTermData.directionIndex == +1 ) ? LUT::HiHiCoeff[faceNormal] : LUT::LoLoCoeff[faceNormal];

    // Divergence term
    floatType ibSource = fvCoeffs.Cont.AU[faceNormal][coeff](G(cellIndex[faceNormal])) * sourceTermData.ghostCellValues.U[faceNormal];

    // Pressure terms
    ibSource += fvCoeffs.Cont.AP[coeff ](G(cellIndex)) * sourceTermData.ghostCellValues.P
              + fvCoeffs.Cont.AP[ccoeff](G(cellIndex)) * sourceTermData.farPressureGhostCellValue;

    fvCoeffs.Cont.B( G(cellIndex) ) += ibSource;
}


void ZeroInSolidStencilCoeffs( FVCoefficients &fvCoeffs,
                               const IBCell::SourceTermData &sourceTermData, 
                               const TensorIndex3D &cellIndex ) 
{
    using FVT::G;

    const Axis::ENUMDATA faceNormal = sourceTermData.direction;
    // const TransportCoefficients::ENUMDATA coeff  = ( sourceTermData.directionIndex == +1 ) ? LUT::HiCoeff[faceNormal]   : LUT::LoCoeff[faceNormal];
    
    // The immediate cell
    // fvCoeffs.Cont.AP[coeff ](G(cellIndex)) = 0.0f;
    if ( fvCoeffs.momentumInterpolation == MomentumInterpolation::Implicit ) {
        const TransportCoefficients::ENUMDATA ccoeff = ( sourceTermData.directionIndex == +1 ) ? LUT::HiHiCoeff[faceNormal] : LUT::LoLoCoeff[faceNormal];
        fvCoeffs.Cont.AP[ccoeff](G(cellIndex)) = 0.0f;

        // The interior cell
        // fvCoeffs.Cont.AP[ccoeff](G(sourceTermData.cellIndex_a)) = 0.0f;
    }

    

}



void InteriorContinuityIBSourceImplicitMWI( FVCoefficients &fvCoeffs,
                                            const IBCell::SourceTermData &sourceTermData ) 
{
    using FVT::G;
    
    const Axis::ENUMDATA faceNormal = sourceTermData.direction;
    const TransportCoefficients::ENUMDATA ccoeff = ( sourceTermData.directionIndex == +1 ) ? LUT::HiHiCoeff[faceNormal] : LUT::LoLoCoeff[faceNormal];

    // Far pressure term
    const floatType ibSource = fvCoeffs.Cont.AP[ccoeff](G(sourceTermData.cellIndex_a)) * sourceTermData.ghostCellValues.P;

    fvCoeffs.Cont.B( G(sourceTermData.cellIndex_a) ) += ibSource;
}



void ContinuityIBSourceSemiExplicitMWI( FVCoefficients &fvCoeffs,
                                        const IBCell::SourceTermData &sourceTermData, 
                                        const TensorIndex3D &cellIndex,
                                        const FieldData<Tensor3D> &fields,
                                        const Mesh &mesh ) 
{
    using FVT::G;

    const bool ghostIsHiSide = ( sourceTermData.directionIndex == +1 );

    const Axis::ENUMDATA faceNormal = sourceTermData.direction;
    const TransportCoefficients::ENUMDATA coeff  = ghostIsHiSide ? LUT::HiCoeff[faceNormal]   : LUT::LoCoeff[faceNormal];

    // Divergence term
    floatType ibSource = fvCoeffs.Cont.AU[faceNormal][coeff](G(cellIndex[faceNormal])) * sourceTermData.ghostCellValues.U[faceNormal];

    // Implicit Pressure terms
    ibSource += fvCoeffs.Cont.AP[coeff ](G(cellIndex)) * sourceTermData.ghostCellValues.P;

    // Explicit Pressure terms, face closest to IB
    const Tensor3D &momentumDiagCoeff = fvCoeffs.Mom[faceNormal].AU[TransportCoefficients::p];
    const std::array<Tensor1D, 4> &mwiSparseCoeffs = fvCoeffs.mwiSparseCoeffs[faceNormal];
    TensorIndex3D loIndex = ghostIsHiSide ? cellIndex : sourceTermData.cellIndex_g;
    TensorIndex3D hiIndex = ghostIsHiSide ? sourceTermData.cellIndex_g : cellIndex;
    intType idx = hiIndex[faceNormal];
    floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeff, mesh, faceNormal );

    floatType ghostSparseCoeff  = ghostIsHiSide ? d * mwiSparseCoeffs[2](idx) : - d * mwiSparseCoeffs[1](idx);
    floatType ghostSparseCCoeff = ghostIsHiSide ? d * mwiSparseCoeffs[3](idx) : - d * mwiSparseCoeffs[0](idx);

    floatType explicitIBSource = ghostSparseCoeff * sourceTermData.ghostCellValues.P
                               + ghostSparseCCoeff * sourceTermData.farPressureGhostCellValue;

    // Remove the effect of the wide stencil term incase there is no ghost cell
    TensorIndex3D farGhostCellIndex = sourceTermData.cellIndex_g;
    farGhostCellIndex[ sourceTermData.direction ] += sourceTermData.directionIndex;
    explicitIBSource -= ghostSparseCCoeff * fields.P( G( farGhostCellIndex ) );

    // Explicit Pressure terms, face farthest from IB
    loIndex = ghostIsHiSide ? sourceTermData.cellIndex_a : cellIndex;
    hiIndex = ghostIsHiSide ? cellIndex : sourceTermData.cellIndex_a;
    idx = hiIndex[faceNormal];
    d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeff, mesh, faceNormal );

    ghostSparseCoeff  = ghostIsHiSide ? - d * mwiSparseCoeffs[3](idx) : d * mwiSparseCoeffs[0](idx);

    explicitIBSource += ghostSparseCoeff * sourceTermData.ghostCellValues.P;


    // Add to the source term, divide by cell length
    ibSource += explicitIBSource * mesh.cellLengthsInv[faceNormal]( cellIndex[faceNormal] );

    fvCoeffs.Cont.B( G(cellIndex) ) += ibSource;

}



void InteriorContinuityIBSourceSemiExplicitMWI( FVCoefficients &fvCoeffs,
                                                const IBCell::SourceTermData &sourceTermData, 
                                                const TensorIndex3D &cellIndex,
                                                const Mesh &mesh ) 
{
    bool ghostIsHiSide = ( sourceTermData.directionIndex == +1 );
    Axis::ENUMDATA faceNormal = sourceTermData.direction;

    // Far pressure term (explicit)
    const TensorIndex3D loIndex = ghostIsHiSide ? sourceTermData.cellIndex_a : cellIndex;
    const TensorIndex3D hiIndex = ghostIsHiSide ? cellIndex : sourceTermData.cellIndex_a;
    const intType idx = hiIndex[faceNormal];
    const Tensor3D &momentumDiagCoeff = fvCoeffs.Mom[faceNormal].AU[TransportCoefficients::p];
    const std::array<Tensor1D, 4> &mwiSparseCoeffs = fvCoeffs.mwiSparseCoeffs[faceNormal];
    const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeff, mesh, faceNormal );

    const floatType ghostSparseCoeff  = ghostIsHiSide ? d * mwiSparseCoeffs[3](idx) : - d * mwiSparseCoeffs[0](idx);

    const floatType ibSource = ghostSparseCoeff * sourceTermData.ghostCellValues.P * mesh.cellLengthsInv[faceNormal]( sourceTermData.cellIndex_a[faceNormal] );

    fvCoeffs.Cont.B( G(sourceTermData.cellIndex_a) ) += ibSource;
}



void ChangeStencilToCentralAtIB( FVCoefficients &fvCoeffs,
                                 const IBData &ibData,
                                 const EnumVector<Axis, Tensor3D> &faceFluxes, 
                                 const Mesh &mesh )
{
    using enum TransportCoefficients::ENUMDATA;

    // All equations share velocity coefficient
    auto &AU = fvCoeffs.Mom[Axis::X].AU;

    for ( auto &ibCellComponent : ibData.ibCells ) {

        for ( auto &ibCell : ibCellComponent ) { 

            TensorIndex3D cellIndex = ibCell.cellIndex;

            for ( auto &sourceTermData : ibCell.sourceTermsData ) {

                const Axis::ENUMDATA faceNormal = sourceTermData.direction;
                TensorIndex3D faceIndex = cellIndex;
                faceIndex[faceNormal] += sourceTermData.faceDirectionIndex;
                const intType fidx = faceIndex[faceNormal],
                              cidx = cellIndex[faceNormal];
                cellIndex = G(cellIndex);   // Coefficients have dummy cells

                const floatType faceFlux = faceFluxes[faceNormal](faceIndex);

                if ( sourceTermData.directionIndex == +1 ) {    // Face on Hi side

                    const TransportCoefficients::ENUMDATA hi = LUT::HiCoeff[faceNormal];

                    // Subtract upwinding term
                    if ( faceFlux >= 0.0f ) {
                        AU[p ](cellIndex) -= faceFlux * mesh.cellLengthsInv[faceNormal]( cidx );
                    } else {
                        AU[hi](cellIndex) -= faceFlux * mesh.cellLengthsInv[faceNormal]( cidx );
                    }

                    // Add in central differencing term
                    AU[p ](cellIndex) += faceFlux * ( 1.0f - mesh.interpFactors[faceNormal]( fidx ) ) * mesh.cellLengthsInv[faceNormal]( cidx );
                    AU[hi](cellIndex) += faceFlux * mesh.interpFactors[faceNormal]( fidx ) * mesh.cellLengthsInv[faceNormal]( cidx );

                } else {                                        // Face on Lo side    

                    const TransportCoefficients::ENUMDATA lo = LUT::LoCoeff[faceNormal];

                    // Subtract upwinding term
                    if ( faceFlux >= 0.0f ) {
                        AU[lo](cellIndex) -= - faceFlux * mesh.cellLengthsInv[faceNormal]( cidx );
                    } else {
                        AU[p ](cellIndex) -= - faceFlux * mesh.cellLengthsInv[faceNormal]( cidx );
                    }

                    // Add in central differencing term
                    AU[p ](cellIndex) += - faceFlux * mesh.interpFactors[faceNormal]( fidx ) * mesh.cellLengthsInv[faceNormal]( cidx );
                    AU[lo](cellIndex) += - faceFlux * ( 1.0f - mesh.interpFactors[faceNormal]( fidx ) ) * mesh.cellLengthsInv[faceNormal]( cidx );

                }
                
            }

        }

    }

}



void AddIBSourceTerms( FVCoefficients &fvCoeffs,
                       const EnumVector<Axis, Tensor3D> &faceFluxes, 
                       const IBData &ibData,
                       const FieldData<Tensor3D> &fields,
                       const Mesh &mesh )
{

    // Set source terms

    for ( auto &ibCellComponent : ibData.ibCells ) {
        for ( auto &ibCell : ibCellComponent ) { 

            const TensorIndex3D cellIndex = ibCell.cellIndex;

            // A source term is added for each forced face
            for ( auto &sourceTermData : ibCell.sourceTermsData ) {

                // Momentum equations
                MomentumIBSourceStencil( fvCoeffs, sourceTermData, cellIndex );

                switch( fvCoeffs.advectionScheme ) {
                    case AdvectionSchemes::Upwind:
                        /* NULL */
                        break;
                    case AdvectionSchemes::Central:
                        MomentumIBSourceDeferredCorrection<AdvectionSchemes::Central>( fvCoeffs, fields, faceFluxes, mesh, sourceTermData, cellIndex );
                        break;
                    case AdvectionSchemes::SOU:
                        MomentumIBSourceDeferredCorrection<AdvectionSchemes::SOU>( fvCoeffs, fields, faceFluxes, mesh, sourceTermData, cellIndex );
                        break;
                    case AdvectionSchemes::QUICK:
                        MomentumIBSourceDeferredCorrection<AdvectionSchemes::QUICK>( fvCoeffs, fields, faceFluxes, mesh, sourceTermData, cellIndex );
                        break;
                }


                // Continuity equation
                switch ( fvCoeffs.momentumInterpolation ) {
                    case MomentumInterpolation::Implicit:
                        ContinuityIBSourceImplicitMWI( fvCoeffs, sourceTermData, cellIndex );
                        InteriorContinuityIBSourceImplicitMWI( fvCoeffs, sourceTermData );
                        break;

                    case MomentumInterpolation::SemiExplicit:
                        ContinuityIBSourceSemiExplicitMWI( fvCoeffs, sourceTermData, cellIndex, fields, mesh );
                        InteriorContinuityIBSourceSemiExplicitMWI( fvCoeffs, sourceTermData, cellIndex, mesh );
                        break;
                }

            }

        }
    }


    // Set in boundary coefficients for the continity equation to zero
    for ( auto &ibCellComponent : ibData.ibCells ) {
        for ( auto &ibCell : ibCellComponent ) { 

            const TensorIndex3D cellIndex = ibCell.cellIndex;
            for ( auto &sourceTermData : ibCell.sourceTermsData ) {

                ZeroInSolidStencilCoeffs( fvCoeffs, sourceTermData, cellIndex );

            }

        }
    }

}



/*---------------------------------------------------------------------------------------------------------------*\
                                                General Functions
\*---------------------------------------------------------------------------------------------------------------*/


void SetGhostCellsToConstant( Tensor3D &array,
                              const floatType value )
{
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        for ( intType i = 0; i != nGhost; i++ ) {

            // Lo side
            array.chip(i, axis).setConstant(value);

            // Hi side
            intType iHi = array.dimension(axis) - 1 - i;
            array.chip(iHi, axis).setConstant(value);

        } 
    } );
}


// Set the coefficients that need to be relinearised to zero
[[maybe_unused]]
void ZeroNonlinearCoeffs( FVCoefficients &fvCoeffs )
{
    using enum TransportCoefficients::ENUMDATA;
    using enum Axis::ENUMDATA;
    
    // Momentum equations
    fvCoeffs.Mom[X].AU[p].setZero();

    // Set the ghost cell central coefficients to 1 to avoid divide by zero in solver
    SetGhostCellsToConstant(fvCoeffs.Mom[X].AU[p], 1.0f);

    EnumFor<Axis> ( [&] (Axis::ENUMDATA axis) {
        fvCoeffs.Mom[axis].B.setZero();
        fvCoeffs.Mom[axis].F.setZero();
    } );

    
    // Continuity equation
    fvCoeffs.Cont.AP[p].setZero();

    fvCoeffs.Cont.B.setZero();
    fvCoeffs.Cont.F.setZero();
}


// Set the coefficients that need to be relinearised to zero
void ZeroNonlinearCoeffsParallel( FVCoefficients &fvCoeffs )
{
    using enum TransportCoefficients::ENUMDATA;
    using enum Axis::ENUMDATA;
    
    // Each task will be done by a seperate thread
    std::vector< std::function<void()> > tasks;

    // Momentum equations 
    tasks.emplace_back( [&] () { 
        fvCoeffs.Mom[X].AU[p].setZero();
        SetGhostCellsToConstant(fvCoeffs.Mom[X].AU[p], 1.0f); // To avoid divide by zero in solver
    } );
    
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        tasks.emplace_back( [&, axis] () { fvCoeffs.Mom[axis].B.setZero(); } );
        tasks.emplace_back( [&, axis] () { fvCoeffs.Mom[axis].F.setZero(); } );
    } );
    
    // Continuity equation
    tasks.emplace_back( [&] () { fvCoeffs.Cont.AP[p].setZero(); } );
    tasks.emplace_back( [&] () { fvCoeffs.Cont.B.setZero(); } );
    tasks.emplace_back( [&] () { fvCoeffs.Cont.F.setZero(); } );

    // Parallel execution
    size_t nTasks = tasks.size();
    RAJA::forall<RAJA::omp_parallel_for_exec>( 
        RAJA::TypedRangeSegment<size_t>(0, nTasks), 
        [&] (size_t i) {
            tasks[i]();
        } 
    );
}



}   // end anonymous namespace





/*---------------------------------------------------------------------------------------------------------------*\
                                            Set and Update Functions
\*---------------------------------------------------------------------------------------------------------------*/

FVCoefficients InitialiseFVCoefficients( const Mesh &mesh,
                                         const InputData &inputData)
{
    // Default construct the coefficients class
    FVCoefficients fvCoeffs( mesh.nCells, inputData.schemes.momentumInterpolation );

    fvCoeffs.rho = inputData.rho;
    fvCoeffs.nu  = inputData.nu;

    // Parts of momentum/continuity equation that don't change with linearisation
    fvCoeffs.advectionScheme         = inputData.schemes.advectionScheme;
    fvCoeffs.timeScheme              = inputData.schemes.timeScheme;
    fvCoeffs.timeStep                = inputData.schemes.timeStep;
    fvCoeffs.advectionBlendingFactor = inputData.schemes.advectionBlendingFactor;

    SetDiffusionCoeffients(fvCoeffs, mesh);
    SetCellGradientCoefficients(fvCoeffs, mesh);
    SetHighOrderAdvectionCoefficients(fvCoeffs, mesh);

    SetMomentumInterpolationSparseConstants(fvCoeffs, mesh);
    SetMomentumInterpolationCompactConstants(fvCoeffs, mesh);

    return fvCoeffs;
}



// Update linearisation in momenum and continuity equations
void UpdateFVCoefficients( FVCoefficients &fvCoeffs, 
                           const Mesh &mesh,
                           const FieldData<Tensor3D> &fields,
                           const FieldData<Tensor3D> &fieldsPrevTime,
                           const FieldData<Tensor3D> &fieldsPrevPrevTime,
                           const EnumVector<Axis, Tensor3D> &faceFluxes,
                           const IBData &ibData )
{
    TIC("Zero nonlinear coefficients")
    ZeroNonlinearCoeffsParallel(fvCoeffs);
    TOC()

    TIC("Set advection")
    SetAdvectionCoefficientsParallel(fvCoeffs, fields, faceFluxes, mesh);
    TOC()

    TIC("Add diffusion")
    AddDiffusionParallel(fvCoeffs, mesh);
    TOC()

    TIC("Change stencil to central at IB")
    ChangeStencilToCentralAtIB( fvCoeffs, ibData, faceFluxes, mesh );
    TOC()

    TIC("Add unsteady term")
    AddUnsteadyTerm(fvCoeffs, fieldsPrevTime, fieldsPrevPrevTime, mesh);
    TOC()

    TIC("Set MWI coefficients")
    SetMomentumInterpolationCoefficientsParallel(fvCoeffs, mesh, fields.P);
    TOC()

    TIC("Add IB source terms")
    AddIBSourceTerms(fvCoeffs, faceFluxes, ibData, fields, mesh);
    TOC()
}


}   // end namespace CFD
