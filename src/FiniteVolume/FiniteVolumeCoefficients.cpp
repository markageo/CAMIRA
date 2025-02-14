#include "FiniteVolume.h"
#include "../Core/Macros.h"
#include "../Core/FVTools.h"
#include "../Core/FVLookups.h"
#include "FaceInterpolatedVelocity.h"

#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

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


// Return true if all values in array are the same value
bool IsConstantArray( const Tensor2D &array )
{
    const floatType testValue = array(0, 0);
    const Eigen::Tensor<bool, 0> isConstant = ( array == testValue ).all();
    return isConstant(0);
}



// Check if continuity equation implies a zero gradient boundary condition. This occurs if both orthogonal fields have a uniform BC
BoundaryConditions::ENUMDATA GetDiffusionBC( const EnumVector< Axis, BoundaryConditionData::Patches > &momBoundaryPatches, 
                                             const BoundaryPatches::ENUMDATA boundaryPatch, 
                                             const Axis::ENUMDATA velocityComponent )
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;

    const Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[boundaryPatch];

    // Set the field we need to check based on the axis
    const Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[ axis ],
                         axis2 = LUT::LoOrthogonalAxis[ axis ];

    // Only check the field that in the direction of the current axis
    if (velocityComponent == axis) {

        // Only possible if we have a fixed velocity BC
        if (momBoundaryPatches[axis1][boundaryPatch].type == BC::fixed && 
            momBoundaryPatches[axis2][boundaryPatch].type == BC::fixed) {

            // Gradient is only zero if boundary value is the same all over
            if ( IsConstantArray( momBoundaryPatches[axis1][boundaryPatch].value ) &&
                 IsConstantArray( momBoundaryPatches[axis2][boundaryPatch].value ) ) {
                return BC::zeroGradient;
            }     

        }
    }

    return momBoundaryPatches[velocityComponent][boundaryPatch].type;
}


// Apply boundary conditions for diffusion terms on axis positive boundary
void DiffusionPositiveBoundary( EnumVector< Axis,  EnumVector<TransportCoefficients, Tensor1D> > &diff, 
                                EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                const Mesh &mesh,  
                                const BoundaryConditionData::Patches &bcDataPatches,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = LUT::LoCoeff[axis],
                                          east = LUT::HiCoeff[axis];
    const intType iCellBound = mesh.nCells(axis) - 1;

    switch ( bcDataPatches[boundaryPatch].type ) {
        
        case BC::zeroGradient: 
            /* NULL */
            break;

        case BC::fixed:
            diff[axis][p   ](iCellBound)     +=   2*mesh.cellLengthsInv[axis](iCellBound);
            boundaryConstants[boundaryPatch] += - 2*mesh.cellLengthsInv[axis](iCellBound);
            break;

        case BC::extrapolated:
            diff[axis][p   ](iCellBound) += - 2*mesh.cellLengthsInv[axis](iCellBound) * (mesh.extrapFactors[boundaryPatch].p - 1);  
            diff[axis][west](iCellBound) += - 2*mesh.cellLengthsInv[axis](iCellBound) * mesh.extrapFactors[boundaryPatch].a;
            break;

        case BC::periodic:
            diff[axis][p   ](iCellBound) +=   mesh.cellCenterDiffInv[axis](iCellBound + 1); 
            diff[axis][east](iCellBound) += - mesh.cellCenterDiffInv[axis](iCellBound + 1);           
           break;

        default:
            break;
    }

}


// Apply boundary conditions for diffusion terms on axis negative boundary
void DiffusionNegativeBoundary( EnumVector< Axis, EnumVector<TransportCoefficients, Tensor1D> > &diff, 
                                EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                const Mesh &mesh,  
                                const BoundaryConditionData::Patches &bcDataPatches,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::NegativePatch[axis];
    const TransportCoefficients::ENUMDATA west = LUT::LoCoeff[axis],
                                          east = LUT::HiCoeff[axis];
    const intType iCellBound = 0;

    switch ( bcDataPatches[boundaryPatch].type ) {
        
        case BC::zeroGradient: 
            /* NULL */
            break;

        case BC::fixed:
            diff[axis][p   ](iCellBound)     +=   2*mesh.cellLengthsInv[axis](iCellBound);
            boundaryConstants[boundaryPatch] += - 2*mesh.cellLengthsInv[axis](iCellBound);
            break;

        case BC::extrapolated:
            diff[axis][p   ](iCellBound) += - 2*mesh.cellLengthsInv[axis](iCellBound) * (mesh.extrapFactors[boundaryPatch].p - 1);
            diff[axis][east](iCellBound) += - 2*mesh.cellLengthsInv[axis](iCellBound) * mesh.extrapFactors[boundaryPatch].a;
            break;

        case BC::periodic:
            diff[axis][p   ](iCellBound) +=   mesh.cellCenterDiffInv[axis](0); 
            diff[axis][west](iCellBound) += - mesh.cellCenterDiffInv[axis](0);   
            break;

        default:
            break;
    }

}


// Set diffusion coefficients for a given momentum equation
void SetDiffusionCoeffients( MomentumEquation &momentumEquation, 
                             const BoundaryConditionData &bcData, 
                             const floatType nu,
                             const Mesh &mesh )
{

    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const Axis::ENUMDATA velocityComponent = momentumEquation.component;
    
    TransportCoefficients::ENUMDATA east, west;     // These are just names, they can be north, south etc.
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    BoundaryConditions::ENUMDATA positivePatchBC, negativePatchBC;  // Store these since the continuity equation can override a BC to be zeroGradient

    auto &diff = momentumEquation.diff;
    auto &diffBoundary = momentumEquation.diffBoundary;

    // Diffusion in each axis is calculated in the same way
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        positivePatch = LUT::PositivePatch[axis];
        negativePatch = LUT::NegativePatch[axis];
        east = LUT::HiCoeff[axis];
        west = LUT::LoCoeff[axis];     

        // Internal faces
        for (intType i = 1; i != mesh.nCells(axis); i++) {
            
            // Cell on west side
            diff[axis][p   ](i-1) +=   mesh.cellCenterDiffInv[axis](i);
            diff[axis][east](i-1) += - mesh.cellCenterDiffInv[axis](i);

            // Cell on east side
            diff[axis][p   ](i)   +=   mesh.cellCenterDiffInv[axis](i);
            diff[axis][west](i)   += - mesh.cellCenterDiffInv[axis](i);

        }


        // Check boundary condition from continuity condition
        positivePatchBC = GetDiffusionBC(bcData.fields.U, positivePatch, velocityComponent);
        negativePatchBC = GetDiffusionBC(bcData.fields.U, negativePatch, velocityComponent); 


        // Boundary conditions only need to be set if it is not zero gradient
        if (positivePatchBC != BC::zeroGradient) {
            DiffusionPositiveBoundary(diff, diffBoundary, mesh, bcData.fields.U[velocityComponent], axis);
        }

        if (negativePatchBC != BC::zeroGradient) {
            DiffusionNegativeBoundary(diff, diffBoundary, mesh, bcData.fields.U[velocityComponent], axis);
        }


        // Multiply by inverse cell length
        for (intType i = 0; i != mesh.nCells(axis); i++) {
            diff[axis][p   ](i) *= mesh.cellLengthsInv[axis](i);
            diff[axis][east](i) *= mesh.cellLengthsInv[axis](i);
            diff[axis][west](i) *= mesh.cellLengthsInv[axis](i);
        }
        diffBoundary[ LUT::PositivePatch[axis] ] *= mesh.cellLengthsInv[axis]( mesh.nCells(axis)-1 );
        diffBoundary[ LUT::NegativePatch[axis] ] *= mesh.cellLengthsInv[axis]( 0 );

        // Multiply by viscosity
        for (intType i = 0; i != mesh.nCells(axis); i++) {
            diff[axis][p   ](i) *= nu;
            diff[axis][east](i) *= nu;
            diff[axis][west](i) *= nu;
        }
        diffBoundary[ LUT::PositivePatch[axis] ] *= nu;
        diffBoundary[ LUT::NegativePatch[axis] ] *= nu;

    } );
}





/*---------------------------------------------------------------------------------------------------------------*\
                                       Momentum Picard Advection Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


void SetHighOrderAdvectionCoefficients( MomentumEquation &momentumEquation,
                                        const Mesh &mesh )
{
    auto &negativeFluxHiOrderAdvectionCoeffs = momentumEquation.negativeFluxHiOrderAdvectionCoeffs;
    auto &positiveFluxHiOrderAdvectionCoeffs = momentumEquation.positiveFluxHiOrderAdvectionCoeffs;

    switch ( momentumEquation.advectionScheme ) {

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
void InteriorAdvectionTerms( MomentumEquation &momentumEquation, 
                             const FieldData<Tensor3D> &fields,
                             const EnumVector<Axis, Tensor3D> &faceFluxes, 
                             const Mesh &mesh,
                             const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    constexpr bool hasDeferredCorrection = ( advectionScheme != AdvectionSchemes::Upwind );

    auto &coeffs     = momentumEquation.AU[momentumEquation.component];
    auto &sourceTerm = momentumEquation.B;
    const auto &U          = fields.U[momentumEquation.component];

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

                floatType highOrderAdvectedVelocity{0.0f},
                          upwindAdvectedVelocity{0.0f};

                if ( faceFlux >= 0.0f ) {
                    coeffs[p   ](G(loIndex)) +=   faceFlux * mesh.cellLengthsInv[axis]( loIndex[axis] );
                    coeffs[east](G(loIndex))  =   0.0f;
                    coeffs[p   ](G(hiIndex)) +=   0.0f;
                    coeffs[west](G(hiIndex))  = - faceFlux * mesh.cellLengthsInv[axis]( hiIndex[axis] );

                    if constexpr ( hasDeferredCorrection ){
                        highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme, +1>(U, momentumEquation, mesh, axis, hiIndex, loIndex);
                        upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, +1>(U, momentumEquation, mesh, axis, hiIndex, loIndex);
                    }
                        

                } else {
                    coeffs[p   ](G(loIndex)) +=   0.0f;
                    coeffs[east](G(loIndex))  =   faceFlux * mesh.cellLengthsInv[axis]( loIndex[axis] );
                    coeffs[p   ](G(hiIndex)) += - faceFlux * mesh.cellLengthsInv[axis]( hiIndex[axis] );
                    coeffs[west](G(hiIndex))  =   0.0f;

                    if constexpr ( hasDeferredCorrection ) {
                        highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme, -1>(U, momentumEquation, mesh, axis, hiIndex, loIndex);
                        upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, -1>(U, momentumEquation, mesh, axis, hiIndex, loIndex);
                    }
                }


                if constexpr ( !hasDeferredCorrection )
                    continue;

                // Deferred correction term
                sourceTerm(G(loIndex)) +=   momentumEquation.advectionBlendingFactor 
                                       *   faceFlux
                                       *   ( highOrderAdvectedVelocity - upwindAdvectedVelocity ) 
                                       *   mesh.cellLengthsInv[axis]( loIndex[axis] );
                    
                sourceTerm(G(hiIndex)) += - momentumEquation.advectionBlendingFactor 
                                       *   faceFlux
                                       *   ( highOrderAdvectedVelocity - upwindAdvectedVelocity ) 
                                       *   mesh.cellLengthsInv[axis]( hiIndex[axis] );
            }
        }
    }
}



void AdvectionPositiveBoundary( EnumVector<TransportCoefficients, Tensor3D> &coeffs, 
                                EnumVector<BoundaryPatches, Tensor2D> &boundaryConstants,
                                const EnumVector<Axis, Tensor3D> &laggedVelocity, 
                                const Mesh &mesh,  
                                const BoundaryConditionData::Patches &bcDataPatches,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = LUT::LoCoeff[axis],
                                          east = LUT::HiCoeff[axis];
    const intType iCellBound = mesh.nCells(axis) - 1;   // Index of cell at the boundary
    const intType iFaceBound = iCellBound + 1;          // Index of face at the boundary
    const TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                        extents = {mesh.nCells[0], mesh.nCells[1], mesh.nCells[2]};

    switch ( bcDataPatches[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p].slice(offsets, extents).chip(iCellBound, axis) += laggedVelocity[axis].chip(iFaceBound, axis) 
                                                                      * laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::fixed:
            boundaryConstants[boundaryPatch]  += laggedVelocity[axis].chip(iFaceBound, axis)
                                               * bcDataPatches[boundaryPatch].value
                                               * laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::extrapolated:
            coeffs[p   ].slice(offsets, extents).chip(iCellBound, axis) += laggedVelocity[axis].chip(iFaceBound, axis) 
                                                                         * laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].p * mesh.cellLengthsInv[axis](iCellBound) );
            coeffs[west].slice(offsets, extents).chip(iCellBound, axis) += laggedVelocity[axis].chip(iFaceBound, axis) 
                                                                         * laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].a * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::periodic:
            coeffs[p   ].slice(offsets, extents).chip(iCellBound, axis) += laggedVelocity[axis].chip(iFaceBound, axis).constant( ( 1 - mesh.interpFactors[axis](iFaceBound) ) *  mesh.cellLengthsInv[axis](iCellBound) )
                                                                         * laggedVelocity[axis].chip(iFaceBound, axis);
            coeffs[east].slice(offsets, extents).chip(iCellBound, axis) += laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.interpFactors[axis](iFaceBound) * mesh.cellLengthsInv[axis](iCellBound) )
                                                                         * laggedVelocity[axis].chip(iFaceBound, axis);
            break;

        default:
            break;
    }

}


void AdvectionNegativeBoundary( EnumVector<TransportCoefficients, Tensor3D> &coeffs, 
                                EnumVector<BoundaryPatches, Tensor2D> &boundaryConstants,
                                const EnumVector<Axis, CFD::Tensor3D> &laggedVelocity, 
                                const Mesh &mesh,  
                                const BoundaryConditionData::Patches &bcDataPatches,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::NegativePatch[axis];
    const TransportCoefficients::ENUMDATA west = LUT::LoCoeff[axis],
                                          east = LUT::HiCoeff[axis];
    const intType iCellBound = 0;   // Index of cell at the boundary 
    const intType iFaceBound = 0;   // Index of face at the boundary
    const TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                        extents = {mesh.nCells[0], mesh.nCells[1], mesh.nCells[2]};

    switch ( bcDataPatches[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p].slice(offsets, extents).chip(iCellBound, axis) += - laggedVelocity[axis].chip(iFaceBound, axis) 
                                                                      *   laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::fixed:
            boundaryConstants[boundaryPatch]  += - laggedVelocity[axis].chip(iFaceBound, axis)
                                               *   bcDataPatches[boundaryPatch].value 
                                               *   laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::extrapolated:
            coeffs[p   ].slice(offsets, extents).chip(iCellBound, axis) += - laggedVelocity[axis].chip(iFaceBound, axis) 
                                                                         *   laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].p * mesh.cellLengthsInv[axis](iCellBound) );
            coeffs[east].slice(offsets, extents).chip(iCellBound, axis) += - laggedVelocity[axis].chip(iFaceBound, axis) 
                                                                         *   laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].a * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::periodic:
            coeffs[p   ].slice(offsets, extents).chip(iCellBound, axis) += - laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.interpFactors[axis](iFaceBound) *  mesh.cellLengthsInv[axis](iCellBound) )
                                                                         *   laggedVelocity[axis].chip(iFaceBound, axis);
            coeffs[west].slice(offsets, extents).chip(iCellBound, axis) += - laggedVelocity[axis].chip(iFaceBound, axis).constant( ( 1 - mesh.interpFactors[axis](iFaceBound) ) * mesh.cellLengthsInv[axis](iCellBound) )
                                                                         *   laggedVelocity[axis].chip(iFaceBound, axis);
            break;

        default:
            break;
    }

}



void SetInteriorAdvectionPicardCoefficients( MomentumEquation &momentumEquation,
                                             const FieldData<Tensor3D> &fields,
                                             const EnumVector<Axis, Tensor3D> &faceFluxes,
                                             const Mesh &mesh )
{
    using enum TransportCoefficients::ENUMDATA;
    
    switch ( momentumEquation.advectionScheme ) {

        case AdvectionSchemes::Upwind:
            EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {
                InteriorAdvectionTerms<AdvectionSchemes::Upwind>(momentumEquation, fields, faceFluxes, mesh, axis);
            } );
            break;

        case AdvectionSchemes::Central:
            EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {
                InteriorAdvectionTerms<AdvectionSchemes::Central>(momentumEquation, fields, faceFluxes, mesh, axis);
            } );
            break;

        case AdvectionSchemes::SOU:
            EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {
                InteriorAdvectionTerms<AdvectionSchemes::SOU>(momentumEquation, fields, faceFluxes, mesh, axis);
            } );
            break;

        case AdvectionSchemes::QUICK:
            EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {
                InteriorAdvectionTerms<AdvectionSchemes::QUICK>(momentumEquation, fields, faceFluxes, mesh, axis);
            } );
            break;

    }
}



void SetBoundaryAdvectionPicardCoefficients( MomentumEquation &momentumEquation,
                                             const EnumVector<Axis, Tensor3D> &faceFluxes, 
                                             const BoundaryConditionData::Patches &bcDataPatches,
                                             const Mesh &mesh )
{
    using enum TransportCoefficients::ENUMDATA;
    
    auto &coeffs            = momentumEquation.AU[ momentumEquation.component ];
    auto &boundaryConstants = momentumEquation.BUBoundary;

    // Upwind internal faces
    EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {

        // Boundary faces
        AdvectionPositiveBoundary(coeffs, boundaryConstants, faceFluxes, mesh, bcDataPatches, axis);
        AdvectionNegativeBoundary(coeffs, boundaryConstants, faceFluxes, mesh, bcDataPatches, axis);

    } );

}





/*---------------------------------------------------------------------------------------------------------------*\
                                       Momentum Newton Advection Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

[[maybe_unused]]
void NewtonInteriorImplicit( EnumVector< TransportCoefficients, Tensor3D > &coeffs, 
                             const EnumVector< Axis, Tensor3D > &faceAdvectedVelocities,  
                             const Mesh &mesh,
                             const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis], 
                                          west = LUT::LoCoeff[axis];

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                TensorIndex3D hiIndex = { i, j, k },
                              loIndex = { i, j, k };
                loIndex[axis] -= 1;
                const intType idx = hiIndex[axis];

                const floatType coeffLo = faceAdvectedVelocities[axis](i, j, k) * mesh.cellLengthsInv[axis]( loIndex[axis] );
                loIndex = G(loIndex);
                coeffs[p   ](loIndex) += coeffLo * ( 1 - mesh.interpFactors[axis]( idx ) );
                coeffs[east](loIndex) += coeffLo * mesh.interpFactors[axis]( idx );

                const floatType coeffHi = faceAdvectedVelocities[axis](i, j, k) * mesh.cellLengthsInv[axis]( hiIndex[axis] );
                hiIndex = G(hiIndex);
                coeffs[p   ](hiIndex) += - coeffHi * mesh.interpFactors[axis]( idx );
                coeffs[west](hiIndex) += - coeffHi * ( 1 - mesh.interpFactors[axis]( idx ) );

            }
        }
    }
}



// Templated to allow compiler autovectorisation
template< Axis::ENUMDATA axis > [[maybe_unused]]
void NewtonInteriorImplicit_autoVec( EnumVector< TransportCoefficients, Tensor3D > &coeffs, 
                                     const EnumVector< Axis, Tensor3D > &faceAdvectedVelocities,  
                                     const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis], 
                                          west = LUT::LoCoeff[axis];
    
    // These will get set at different levels of the nested loop depending on what axis is
    floatType LoCellLengthInv, HiCellLengthInv,
              LoInterpFactor , HiInterpFactor;


    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {

        if constexpr ( axis == Z ) {
            LoCellLengthInv = mesh.cellLengthsInv[axis]( k-1 );
            HiCellLengthInv = mesh.cellLengthsInv[axis]( k   );
            LoInterpFactor  = mesh.interpFactors[axis]( k-1 );
            HiInterpFactor  = mesh.interpFactors[axis]( k   );
        }

        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {

            if constexpr ( axis == Y ) {
                LoCellLengthInv = mesh.cellLengthsInv[axis]( j-1 );
                HiCellLengthInv = mesh.cellLengthsInv[axis]( j   );
                LoInterpFactor  = mesh.interpFactors[axis]( j-1 );
                HiInterpFactor  = mesh.interpFactors[axis]( j   );
            }

            // Left side cells
            // CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                if constexpr ( axis == X ) {
                    LoCellLengthInv = mesh.cellLengthsInv[axis]( i-1 );
                    LoInterpFactor  = mesh.interpFactors[axis]( i-1 );
                }

                TensorIndex3D loIndex = { i, j, k };
                loIndex[axis] -= 1; 
                
                const floatType coeffLo = faceAdvectedVelocities[axis](i, j, k) * LoCellLengthInv;
                coeffs[p   ](G(loIndex)) += coeffLo * ( 1 - LoInterpFactor );
                coeffs[east](G(loIndex)) += coeffLo * LoInterpFactor;

            }

            // Right side cells
            // CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                if constexpr ( axis == X ) {
                    HiCellLengthInv = mesh.cellLengthsInv[axis]( i );
                    HiInterpFactor  = mesh.interpFactors[axis]( i );
                }

                TensorIndex3D hiIndex = { i, j, k };

                const floatType coeffHi = faceAdvectedVelocities[axis](i, j, k) * HiCellLengthInv;
                coeffs[p   ](G(hiIndex)) += - coeffHi * HiInterpFactor;
                coeffs[west](G(hiIndex)) += - coeffHi * ( 1 - HiInterpFactor );

            }
        }
    }

}




void NewtonConstants( Tensor3D &B,
                      const EnumVector< Axis, Tensor3D > &faceAdvectedVelocities,
                      const EnumVector< Axis, Tensor3D > &faceFluxes,
                      const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {

            CFD_PRAGMA_VECTORIZE
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                const floatType xFluxDiff = - mesh.cellLengthsInv[X](i) 
                                            * ( faceFluxes[X](i+1, j, k) * faceAdvectedVelocities[X](i+1, j, k) 
                                            - faceFluxes[X](i  , j, k) * faceAdvectedVelocities[X](i  , j, k) );

                const floatType yFluxDiff = - mesh.cellLengthsInv[Y](j) 
                                            * ( faceFluxes[Y](i, j+1, k) * faceAdvectedVelocities[Y](i, j+1, k) 
                                            - faceFluxes[Y](i, j  , k) * faceAdvectedVelocities[Y](i, j  , k) );

                const floatType zFluxDiff = - mesh.cellLengthsInv[Z](k) 
                                            * ( faceFluxes[Z](i, j, k+1) * faceAdvectedVelocities[Z](i, j, k+1) 
                                            - faceFluxes[Z](i, j, k  ) * faceAdvectedVelocities[Z](i, j, k  ) );

                B( G(i, j, k) ) += xFluxDiff + yFluxDiff + zFluxDiff;

            }
        }
    }

}



void AddAdvectionNewtonCoefficients( MomentumEquation &momentumEquation,
                                     const EnumVector< Axis, EnumVector<Axis, Tensor3D> > &faceAdvectedVelocities, 
                                     const EnumVector< Axis, Tensor3D> &faceFluxes,
                                     const EnumVector< Axis, BoundaryConditionData::Patches > &momBoundaryPatches,
                                     const Mesh &mesh )
{
    const auto &faceVelComp = faceAdvectedVelocities[momentumEquation.component];
    auto &boundaryConstants = momentumEquation.BUBoundary;
    
    // Implicit terms
    #if defined( CFD_USE_AUTOVEC_FUNCTIONS )
        using enum Axis::ENUMDATA;
        NewtonInteriorImplicit_autoVec<Axis::X>(momentumEquation.AU[X], faceVelComp, mesh );
        NewtonInteriorImplicit_autoVec<Axis::Y>(momentumEquation.AU[Y], faceVelComp, mesh );
        NewtonInteriorImplicit_autoVec<Axis::Z>(momentumEquation.AU[Z], faceVelComp, mesh );
    #else
        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
            NewtonInteriorImplicit(momentumEquation.AU[axis], faceVelComp, mesh, axis );
        } );
    #endif

    // Implicit boundary terms
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        auto &coeffs            = momentumEquation.AU[ axis ];
        AdvectionPositiveBoundary(coeffs, boundaryConstants, faceVelComp, mesh, momBoundaryPatches[axis], axis);
        AdvectionNegativeBoundary(coeffs, boundaryConstants, faceVelComp, mesh, momBoundaryPatches[axis], axis);
    } );

    // Explicit terms
    NewtonConstants( momentumEquation.B, faceVelComp, faceFluxes, mesh );
}





/*---------------------------------------------------------------------------------------------------------------*\
                                           Add Diffusion Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


void AddDiffusion( MomentumEquation &momentumEquation,
                   const BoundaryConditionData::Patches &bcDataPatches,
                   const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    auto &velCoeffs   = momentumEquation.AU[ momentumEquation.component ];
    auto &boundaryVel = momentumEquation.BUBoundary;

    const auto &diffCoeffs         = momentumEquation.diff;
    const auto &diffCoeffsBoundary = momentumEquation.diffBoundary;

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

    // Constant terms
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA patch) {
        if ( bcDataPatches[patch].type == BoundaryConditions::fixed ) {
            boundaryVel[patch] += boundaryVel[patch].constant( diffCoeffsBoundary[patch] )
                                * bcDataPatches[patch].value;
        }
    } );

}





/*---------------------------------------------------------------------------------------------------------------*\
                                        Linear Interpolated Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


void FaceInterpolatedPositiveBoundaryStencil( EnumVector< TransportCoefficients, Tensor1D > &coeffs, 
                                              EnumVector< BoundaryPatches, Tensor2D > &boundaryConstants,
                                              const Mesh &mesh,  
                                              const BoundaryConditionData::Patches &bcDataPatches,
                                              const Axis::ENUMDATA axis)
{

    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::PositivePatch[axis];
    const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis],
                                          west = LUT::LoCoeff[axis];
    const intType iCellBound = mesh.nCells(axis) - 1;

    switch ( bcDataPatches[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p   ]( G(iCellBound) ) += 1;
            coeffs[west]( G(iCellBound) ) += 0;
            break;

        case BC::fixed:
            coeffs[p   ]( G(iCellBound) ) += 0;
            coeffs[west]( G(iCellBound) ) += 0;
            boundaryConstants[boundaryPatch] = bcDataPatches[boundaryPatch].value; 
            break;

        case BC::extrapolated:
            coeffs[p   ]( G(iCellBound) ) += mesh.extrapFactors[boundaryPatch].p;
            coeffs[west]( G(iCellBound) ) += mesh.extrapFactors[boundaryPatch].a;
            break;

        case BC::periodic:
            coeffs[p   ]( G(iCellBound) ) += 1 - mesh.interpFactors[axis](iCellBound + 1);
            coeffs[east]( G(iCellBound) ) += mesh.interpFactors[axis](iCellBound + 1);
            break;

        default:
            break;
    }

}


void FaceInterpolatedNegativeBoundaryStencil( EnumVector< TransportCoefficients, Tensor1D > &coeffs, 
                                              EnumVector< BoundaryPatches, Tensor2D > &boundaryConstants,
                                              const Mesh &mesh,  
                                              const BoundaryConditionData::Patches &bcDataPatches, 
                                              const Axis::ENUMDATA axis)
{

    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::NegativePatch[axis];
    const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis],
                                          west = LUT::LoCoeff[axis];
    const intType iCellBound = 0;

    switch ( bcDataPatches[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p   ]( G(iCellBound) ) += - 1;
            coeffs[east]( G(iCellBound) ) +=   0;
            break;

        case BC::fixed:
            coeffs[p   ]( G(iCellBound) ) += 0;
            coeffs[east]( G(iCellBound) ) += 0;
            boundaryConstants[boundaryPatch] = - bcDataPatches[boundaryPatch].value; 
            break;

        case BC::extrapolated:
            coeffs[p   ]( G(iCellBound) ) += - mesh.extrapFactors[boundaryPatch].p;
            coeffs[east]( G(iCellBound) ) += - mesh.extrapFactors[boundaryPatch].a;
            break;

        case BC::periodic:
            coeffs[p   ]( G(iCellBound) ) += - mesh.interpFactors[axis](0);
            coeffs[west]( G(iCellBound) ) += - ( 1 - mesh.interpFactors[axis](0) );
            break;

        default:
            break;
    }

}


void AddFaceInterpolatedBoundarySources( EnumVector< BoundaryPatches, Tensor2D > &boundaryConstants,
                                         const Mesh &mesh,  
                                         const BoundaryConditionData::Patches &bcDataPatches,
                                         const BoundaryPatches::ENUMDATA boundaryPatch )
{
    using BC = BoundaryConditions::ENUMDATA;

    const Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[boundaryPatch];
    const floatType sign = ( boundaryPatch == LUT::PositivePatch[axis] ) ? +1.0f : -1.0f;

    switch ( bcDataPatches[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            /* NULL */
            break;

        case BC::fixed:
            boundaryConstants[boundaryPatch] = boundaryConstants[boundaryPatch].constant(sign) * bcDataPatches[boundaryPatch].value; 
            break;

        case BC::extrapolated:
            /* NULL */
            break;

        case BC::periodic:
        {
            /* NULL */
        }
            
        default:
            break;
    }

    if ( boundaryPatch == LUT::PositivePatch[axis] ) {
        boundaryConstants[ boundaryPatch ] *= boundaryConstants[ boundaryPatch ].constant( mesh.cellLengthsInv[axis]( mesh.nCells(axis)-1 ) );
    } else {
        boundaryConstants[ boundaryPatch ] *= boundaryConstants[ boundaryPatch ].constant( mesh.cellLengthsInv[axis]( 0 ) );
    }
}


// Set coefficients for quantities that are intrpolated linearly onto faces.
void SetFaceInterpolatedCoefficients( EnumVector<TransportCoefficients, Tensor1D> &coeffs, 
                                      EnumVector< BoundaryPatches, Tensor2D > &boundaryConstants, 
                                      const Mesh &mesh, 
                                      const BoundaryConditionData::Patches &bcDataPatches, 
                                      const Axis::ENUMDATA axis)
{ 
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

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

    // Boundary faces
    FaceInterpolatedPositiveBoundaryStencil(coeffs, boundaryConstants, mesh, bcDataPatches, axis);
    FaceInterpolatedNegativeBoundaryStencil(coeffs, boundaryConstants, mesh, bcDataPatches, axis);
    
    // Multiply by inverse of cell length
    for (intType i = 0; i != mesh.nCells(axis); i++) {
        coeffs[p   ](G(i)) *= mesh.cellLengthsInv[axis](i);
        coeffs[east](G(i)) *= mesh.cellLengthsInv[axis](i);
        coeffs[west](G(i)) *= mesh.cellLengthsInv[axis](i); 
    }
    
}


void SetContinuityVelocityCoefficients( ContinuityEquation &continuityEquation, 
                                        const Mesh &mesh, 
                                        const EnumVector< Axis, BoundaryConditionData::Patches > &velBcData)
{
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        SetFaceInterpolatedCoefficients(continuityEquation.AU[axis], continuityEquation.BUBoundary, mesh, velBcData[axis], axis);
    } );   
}


void SetMomentumPressureCoefficients( MomentumEquation &momentumEquation, 
                                      const floatType rho,
                                      const Mesh &mesh, 
                                      const BoundaryConditionData::Patches &bcDataPatches)
{
    using enum TransportCoefficients::ENUMDATA;

    Axis::ENUMDATA axis     = momentumEquation.component;
    const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis],    // These are just names, they can be north, south etc.
                                          west = LUT::LoCoeff[axis];  

    SetFaceInterpolatedCoefficients(momentumEquation.AP, momentumEquation.BPBoundary, mesh, bcDataPatches, axis);

    // Divide coefficients by fluid density
    momentumEquation.AP[p   ] /= momentumEquation.AP[p   ].constant( rho );
    momentumEquation.AP[east] /= momentumEquation.AP[east].constant( rho );
    momentumEquation.AP[west] /= momentumEquation.AP[west].constant( rho );
}


void AddMomentumPressureBoundarySources( MomentumEquation &momentumEquation,
                                         const floatType &rho,
                                         const Mesh &mesh,  
                                         const BoundaryConditionData::Patches &bcDataPatches )
{
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {

        AddFaceInterpolatedBoundarySources(momentumEquation.BPBoundary, mesh, bcDataPatches, bp);

        momentumEquation.BPBoundary[ bp ] /= momentumEquation.BPBoundary[ bp ].constant( rho );
    } );
}


void AddContinuityVelocityBoundarySources( ContinuityEquation &continuityEquation,
                                           const Mesh &mesh,  
                                           const EnumVector<Axis, BoundaryConditionData::Patches > &velocityBoundaryConditions )
{
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        const Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[bp];
        AddFaceInterpolatedBoundarySources( continuityEquation.BUBoundary, mesh, velocityBoundaryConditions[axis], bp);
    } );
}


/*---------------------------------------------------------------------------------------------------------------*\
                        Momentum Weighted Interpolation (Rhie-Chow Interpolation) Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

// Unweighted constants that appear in the sparse gradient part of MWI
void SetMomentumInterpolationSparseConstants( std::array< Tensor1D, 4 > &mwiSparseCoeffs,
                                              const EnumVector<TransportCoefficients, Tensor1D> &momentumPressureCoeffs,
                                              const Mesh &mesh,
                                              const Axis::ENUMDATA axis )
{
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

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
}



// Unweighted constants that appear in the compact gradient part of MWI
void SetMomentumInterpolationCompactConstants( std::array< Tensor1D, 2 > &mwiCompactCoeffs,
                                               const floatType rho,
                                               const Mesh &mesh,
                                               const Axis::ENUMDATA axis )
{
    using enum TransportCoefficients::ENUMDATA;

    const floatType rhoInv = 1 / rho;

    // Internal faces
    for ( intType i = 1; i != mesh.cellFaces[axis].size()-1; i++ ) {

        mwiCompactCoeffs[0](i) =   mesh.cellCenterDiffInv[axis](i) * rhoInv;

        mwiCompactCoeffs[1](i) = - mesh.cellCenterDiffInv[axis](i) * rhoInv;
    }
}



// Cell weighting coefficient for MWI.
__attribute__((always_inline)) 
inline floatType MWIWeightingCoeff( const TensorIndex3D &loIndex,
                             const TensorIndex3D &hiIndex,
                             const Tensor3D &AUUpInv, 
                             const Mesh& mesh,
                             const Axis::ENUMDATA axis)
{
    const intType idx = hiIndex[axis];    // Axis index of the face
    const floatType interpFactor = mesh.interpFactors[axis]( idx );
    return  ( 1.0f - interpFactor ) *  AUUpInv( G(loIndex) )  +  interpFactor * AUUpInv( G(hiIndex) );
}



// Fully implicit momentum interpolation coefficient for internal faces
[[maybe_unused]]
void MWInterpolationInteriorImplicit( ContinuityEquation &continuityEquation,
                                      const Tensor3D &momentumDiagCoeffInv,
                                      const Mesh &mesh,
                                      const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    // Unpack
    EnumVector<TransportCoefficients, Tensor3D> &continuityPressureCoeffs = continuityEquation.AP;
    const std::array< Tensor1D, 4 > &mwiSparseCoeffs                      = continuityEquation.mwiSparseCoeffs[axis];
    const std::array< Tensor1D, 2 > &mwiCompactCoeffs                     = continuityEquation.mwiCompactCoeffs[axis];

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

                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeffInv, mesh, axis );

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
// Autovectorisation friendly version
template< Axis::ENUMDATA axis > [[maybe_unused]]
void MWInterpolationInteriorImplicit_autoVec( ContinuityEquation &continuityEquation,
                                              const Tensor3D &momentumDiagCoeffInv,
                                              const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    using FVT::G;

    // Unpack
    EnumVector<TransportCoefficients, Tensor3D> &continuityPressureCoeffs = continuityEquation.AP;
    const std::array< Tensor1D, 4 > &mwiSparseCoeffs                      = continuityEquation.mwiSparseCoeffs[axis];
    const std::array< Tensor1D, 2 > &mwiCompactCoeffs                     = continuityEquation.mwiCompactCoeffs[axis];

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

    floatType HiCellLengthInv,  // Depending on axis, these will change at different levels of loop nesting. 
              LoCellLengthInv;
    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {

        if constexpr ( axis == Z ) {
            HiCellLengthInv = mesh.cellLengthsInv[axis]( k   );
            LoCellLengthInv = mesh.cellLengthsInv[axis]( k-1 );
        }


        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {

            if constexpr ( axis == Y ) {
                HiCellLengthInv = mesh.cellLengthsInv[axis]( j   );
                LoCellLengthInv = mesh.cellLengthsInv[axis]( j-1 );
            }
            

            // Precalculate the coefficients on each face to save recalculation in the seperate vectorised loops
            std::vector< std::array< floatType, 4 > > mwiLineCoeffs( static_cast<size_t>( nFaces[X] ) );

            // CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D hiIndex = { i, j, k },
                              loIndex = { i, j, k };
                loIndex[axis] -= 1;

                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeffInv, mesh, axis ); 

                // Coefficients for westmost to eastmost cell
                const intType idx = hiIndex[axis];
                const size_t iv = static_cast< size_t >( i );
                mwiLineCoeffs[iv][0] = d * mwiSparseCoeffs[0](idx),
                mwiLineCoeffs[iv][1] = d * ( mwiSparseCoeffs[1](idx) + mwiCompactCoeffs[0](idx) ),
                mwiLineCoeffs[iv][2] = d * ( mwiSparseCoeffs[2](idx) + mwiCompactCoeffs[1](idx) ),
                mwiLineCoeffs[iv][3] = d * mwiSparseCoeffs[3](idx);
            }

            // CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                if constexpr ( axis == X ) {
                    HiCellLengthInv = mesh.cellLengthsInv[axis]( i );
                }

                TensorIndex3D hiIndex = { i, j, k };
                hiIndex = G(hiIndex);
                const size_t iv = static_cast< size_t >( i );

                // Cell on east side
                continuityPressureCoeffs[wwest](hiIndex)  = - mwiLineCoeffs[iv][0] * HiCellLengthInv;
                continuityPressureCoeffs[west ](hiIndex)  = - mwiLineCoeffs[iv][1] * HiCellLengthInv;
                continuityPressureCoeffs[p    ](hiIndex) += - mwiLineCoeffs[iv][2] * HiCellLengthInv;
                continuityPressureCoeffs[east ](hiIndex)  = - mwiLineCoeffs[iv][3] * HiCellLengthInv;
            }

            // CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                if constexpr ( axis == X ) {
                    LoCellLengthInv = mesh.cellLengthsInv[axis]( i-1 );
                }

                TensorIndex3D loIndex = { i, j, k };
                loIndex[axis] -= 1;
                loIndex = G(loIndex);
                const size_t iv = static_cast< size_t >( i );

                // Cell on west side 
                continuityPressureCoeffs[west ](loIndex) += mwiLineCoeffs[iv][0] * LoCellLengthInv;
                continuityPressureCoeffs[p    ](loIndex) += mwiLineCoeffs[iv][1] * LoCellLengthInv;
                continuityPressureCoeffs[east ](loIndex) += mwiLineCoeffs[iv][2] * LoCellLengthInv;
                continuityPressureCoeffs[eeast](loIndex)  = mwiLineCoeffs[iv][3] * LoCellLengthInv;
            }

        }
    }

}



// Semi explicit momentum interpolation coefficient for internal faces
[[maybe_unused]]
void MWInterpolationInteriorSemiExplicit( ContinuityEquation &continuityEquation, 
                                          const Tensor3D &P,
                                          const Tensor3D &momentumDiagCoeffInv,
                                          const Mesh &mesh,
                                          const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    // Unpack
    EnumVector<TransportCoefficients, Tensor3D> &continuityPressureCoeffs = continuityEquation.AP;
    Tensor3D &continuitySourceTerm                                        = continuityEquation.B;
    const std::array< Tensor1D, 4 > &mwiSparseCoeffs                      = continuityEquation.mwiSparseCoeffs[axis];
    const std::array< Tensor1D, 2 > &mwiCompactCoeffs                     = continuityEquation.mwiCompactCoeffs[axis];

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
                 
                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeffInv, mesh, axis );

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
// Autovectorisation friendly version
template< Axis::ENUMDATA axis > [[maybe_unused]]
void MWInterpolationInteriorSemiExplicit_autoVec( ContinuityEquation &continuityEquation, 
                                                  const Tensor3D &P,
                                                  const Tensor3D &momentumDiagCoeffInv,
                                                  const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    // Unpack
    EnumVector<TransportCoefficients, Tensor3D> &continuityPressureCoeffs = continuityEquation.AP;
    Tensor3D &continuitySourceTerm                                        = continuityEquation.B;
    const std::array< Tensor1D, 4 > &mwiSparseCoeffs                      = continuityEquation.mwiSparseCoeffs[axis];
    const std::array< Tensor1D, 2 > &mwiCompactCoeffs                     = continuityEquation.mwiCompactCoeffs[axis];

    // For getting the index of a neighbouring cell
    auto NeighbourIndex = [] ( TensorIndex3D index, intType shift, Axis::ENUMDATA shiftAxis ) { 
        index[shiftAxis] += shift; 
        return index; 
    };

    // Cell indexing
    const TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                          west  = LUT::LoCoeff[axis];

    const auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    floatType HiCellLengthInv,
              LoCellLengthInv;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {

        if constexpr ( axis == Z ) {
            HiCellLengthInv = mesh.cellLengthsInv[axis]( k   );
            LoCellLengthInv = mesh.cellLengthsInv[axis]( k-1 );
        }

        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {

            if constexpr ( axis == Y ) {
                HiCellLengthInv = mesh.cellLengthsInv[axis]( j   );
                LoCellLengthInv = mesh.cellLengthsInv[axis]( j-1 );
            }

            std::vector< std::array< floatType, 4 > > mwiLineSparseCoeffs( static_cast<size_t>( nFaces[X] ) );
            std::vector< std::array< floatType, 2 > > mwiLineCompactCoeffs( static_cast<size_t>( nFaces[X] ) );

            // CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D hiIndex = { i, j, k },
                             loIndex = { i, j, k };
                loIndex[axis] -= 1;
                
                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeffInv, mesh, axis );

                // Coefficients for westmost to eastmost cell
                const intType idx = hiIndex[axis];
                const size_t iv = static_cast< size_t >( i );
                mwiLineCompactCoeffs[iv][0] = d * mwiCompactCoeffs[0](idx),
                mwiLineCompactCoeffs[iv][1] = d * mwiCompactCoeffs[1](idx);

                mwiLineSparseCoeffs[iv][0]  = d * mwiSparseCoeffs[0](idx),
                mwiLineSparseCoeffs[iv][1]  = d * mwiSparseCoeffs[1](idx),
                mwiLineSparseCoeffs[iv][2]  = d * mwiSparseCoeffs[2](idx),
                mwiLineSparseCoeffs[iv][3]  = d * mwiSparseCoeffs[3](idx);

            }


            // Implicit compact difference --------------------------------------------------------------------------

            // CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                if constexpr ( axis == X ) {
                    HiCellLengthInv = mesh.cellLengthsInv[axis]( i );
                }

                TensorIndex3D hiIndex = { i, j, k };
                hiIndex = G(hiIndex);
                const size_t iv = static_cast< size_t >( i );

                // Cell on east side
                continuityPressureCoeffs[west ](hiIndex)  = - mwiLineCompactCoeffs[iv][0] * HiCellLengthInv;
                continuityPressureCoeffs[p    ](hiIndex) += - mwiLineCompactCoeffs[iv][1] * HiCellLengthInv;
            }


            // CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                if constexpr ( axis == X ) {
                    LoCellLengthInv = mesh.cellLengthsInv[axis]( i-1 );
                }

                TensorIndex3D loIndex = { i, j, k };
                loIndex[axis] -= 1;
                loIndex = G(loIndex);
                const size_t iv = static_cast< size_t >( i );

                // Cell on west side 
                continuityPressureCoeffs[p    ](loIndex) += mwiLineCompactCoeffs[iv][0] * LoCellLengthInv;
                continuityPressureCoeffs[east ](loIndex)  = mwiLineCompactCoeffs[iv][1] * LoCellLengthInv;
            }


            // Explicit sparse difference ---------------------------------------------------------------------------

            // CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                if constexpr ( axis == X ) {
                    HiCellLengthInv = mesh.cellLengthsInv[axis]( i );
                }

                TensorIndex3D hiIndex = { i, j, k };
                const size_t iv = static_cast< size_t >( i );

                const TensorIndex3D hiWWest = NeighbourIndex( hiIndex, -2, axis ),
                                    hiWest  = NeighbourIndex( hiIndex, -1, axis ),
                                    hiEast  = NeighbourIndex( hiIndex,  1, axis );

                // Cell on east side
                continuitySourceTerm(G(hiIndex)) -= ( mwiLineSparseCoeffs[iv][0] * P( G(hiWWest) )
                                                    + mwiLineSparseCoeffs[iv][1] * P( G(hiWest)  )
                                                    + mwiLineSparseCoeffs[iv][2] * P( G(hiIndex) )
                                                    + mwiLineSparseCoeffs[iv][3] * P( G(hiEast)  )
                                                    ) * HiCellLengthInv;
            }


            // CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                if constexpr ( axis == X ) {
                    LoCellLengthInv = mesh.cellLengthsInv[axis]( i-1 );
                }

                TensorIndex3D loIndex = { i, j, k };
                loIndex[axis] -= 1;
                const size_t iv = static_cast< size_t >( i );

                const TensorIndex3D loWest  = NeighbourIndex( loIndex, -1, axis ),
                                    loEast  = NeighbourIndex( loIndex,  1, axis ),
                                    loEEast = NeighbourIndex( loIndex,  2, axis );

                // Cell on west side 
                continuitySourceTerm(G(loIndex)) += ( mwiLineSparseCoeffs[iv][0] * P( G(loWest)  )
                                                    + mwiLineSparseCoeffs[iv][1] * P( G(loIndex) )
                                                    + mwiLineSparseCoeffs[iv][2] * P( G(loEast)  )
                                                    + mwiLineSparseCoeffs[iv][3] * P( G(loEEast) )
                                                    ) * LoCellLengthInv;
            }
        }
    }

}


// Boundary constants that come from MWI for negative side boundary
void MWInterpolationNegativeBoundary( EnumVector<BoundaryPatches, Tensor2D> &continuityBoundaryPressure,
                                      const EnumVector<BoundaryPatches, Tensor2D> &momentumBoundaryPressure,
                                      const Tensor3D &momentumDiagCoeffInv, 
                                      const Mesh &mesh,
                                      const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;

    const BoundaryPatches::ENUMDATA negativePatch = LUT::NegativePatch[ axis ];
    const Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[axis],
                         axis2 = LUT::HiOrthogonalAxis[axis];

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);
    startIndex[axis] = 1;
    nFaces[axis] = startIndex[axis] + 1;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D idx = { i, j, k },
                              hiIndex = { i, j, k },
                              loIndex = hiIndex;
                loIndex[axis] -= 1;

                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeffInv, mesh, axis );
                continuityBoundaryPressure[ negativePatch ]( idx[axis1], idx[axis2] ) = d * (1 - mesh.interpFactors[axis]( idx[axis] )) * momentumBoundaryPressure[ negativePatch ]( idx[axis1], idx[axis2] )
                                                                                          * mesh.cellLengthsInv[axis]( idx[axis] );
            }
        }
    }

}



// Boundary constants that come from MWI for positive side boundary
void MWInterpolationPositiveBoundary( EnumVector<BoundaryPatches, Tensor2D> &continuityBoundaryPressure,
                                      const EnumVector<BoundaryPatches, Tensor2D> &momentumBoundaryPressure,
                                      const Tensor3D &momentumDiagCoeffInv, 
                                      const Mesh &mesh,
                                      const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;

    const BoundaryPatches::ENUMDATA positivePatch = LUT::PositivePatch[ axis ];
    const Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[axis],
                         axis2 = LUT::HiOrthogonalAxis[axis];

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);
    startIndex[axis] = mesh.nCells[axis] - 1;
    nFaces[axis] = startIndex[axis] + 1;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D idx = { i, j, k },
                              hiIndex = { i, j, k },
                              loIndex = hiIndex;
                loIndex[axis] -= 1;

                const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeffInv, mesh, axis );
                continuityBoundaryPressure[ positivePatch ]( idx[axis1], idx[axis2] ) = d * mesh.interpFactors[axis]( idx[axis] ) * momentumBoundaryPressure[ positivePatch ]( idx[axis1], idx[axis2] )
                                                                                          * mesh.cellLengthsInv[axis]( idx[axis] );
            }
        }
    }

}




// Set momentum interpolation coefficients
void SetMomentumInterpolationCoefficients( FVCoefficients &fvCoeffs,
                                           const Mesh &mesh,
                                           const BoundaryConditionData &bcData,
                          [[maybe_unused]] const Tensor3D &P )
{
    // Assumes that 'p' coefficient is set to zero
    // Correction is zero at the boundary faces

    #if defined( CFD_USE_AUTOVEC_FUNCTIONS )
        using enum Axis::ENUMDATA;
        switch ( fvCoeffs.Cont.momentumInterpolation ) {

            case MomentumInterpolation::Implicit:
                MWInterpolationInteriorImplicit_autoVec<Axis::X>(fvCoeffs.Cont, fvCoeffs.Mom[X].diagCoeffInv, mesh);
                MWInterpolationInteriorImplicit_autoVec<Axis::Y>(fvCoeffs.Cont, fvCoeffs.Mom[Y].diagCoeffInv, mesh);
                MWInterpolationInteriorImplicit_autoVec<Axis::Z>(fvCoeffs.Cont, fvCoeffs.Mom[Z].diagCoeffInv, mesh);
                break;

            case MomentumInterpolation::SemiExplicit:
                MWInterpolationInteriorSemiExplicit_autoVec<Axis::X>(fvCoeffs.Cont, P, fvCoeffs.Mom[X].diagCoeffInv, mesh);
                MWInterpolationInteriorSemiExplicit_autoVec<Axis::Y>(fvCoeffs.Cont, P, fvCoeffs.Mom[Y].diagCoeffInv, mesh);
                MWInterpolationInteriorSemiExplicit_autoVec<Axis::Z>(fvCoeffs.Cont, P, fvCoeffs.Mom[Z].diagCoeffInv, mesh);
                break;
        }
    #else
        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
            switch ( fvCoeffs.Cont.momentumInterpolation ) {
                case MomentumInterpolation::Implicit:
                    MWInterpolationInteriorImplicit(fvCoeffs.Cont, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);
                    break;

                case MomentumInterpolation::SemiExplicit:
                     MWInterpolationInteriorSemiExplicit(fvCoeffs.Cont, P, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);
                     break;
            }
        } );
    #endif


    // Boundary constants
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        
        if ( bcData.fields.P[ LUT::PositivePatch[axis] ].type == BoundaryConditions::fixed ) 
            MWInterpolationPositiveBoundary(fvCoeffs.Cont.BPBoundary, fvCoeffs.Mom[axis].BPBoundary, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);

        if ( bcData.fields.P[ LUT::NegativePatch[axis] ].type == BoundaryConditions::fixed ) 
            MWInterpolationNegativeBoundary(fvCoeffs.Cont.BPBoundary, fvCoeffs.Mom[axis].BPBoundary, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);

    } );

}


/*---------------------------------------------------------------------------------------------------------------*\
                                                Unsteady Terms
\*---------------------------------------------------------------------------------------------------------------*/


void AddBackwardsEuler( MomentumEquation &momentumEquation, 
                        const Mesh &mesh,
                        const FieldData<Tensor3D> &fieldsPrevTime )
{
    using enum TransportCoefficients::ENUMDATA;
    const Axis::ENUMDATA axis = momentumEquation.component;
    const floatType dtInv = 1.0f / momentumEquation.timeStep;

    for (intType k = 0; k != mesh.nCells(2); k++) {
        for (intType j = 0; j != mesh.nCells(1); j++) {
            for (intType i = 0; i != mesh.nCells(0); i++) {

                momentumEquation.AU[axis][p]( G(i, j, k) ) += dtInv;
                momentumEquation.B( G(i, j, k) )           += - fieldsPrevTime.U[axis]( G(i, j, k) ) * dtInv;

            }
        }
    }
}



void AddBackwardsThreeLevel( MomentumEquation &momentumEquation, 
                             const Mesh &mesh,
                             const FieldData<Tensor3D> &fieldsPrevTime,
                             const FieldData<Tensor3D> &fieldsPrevPrevTime )
{
    using enum TransportCoefficients::ENUMDATA;
    const Axis::ENUMDATA axis = momentumEquation.component;
    const floatType dtInv = 1.0f / momentumEquation.timeStep;

    for (intType k = 0; k != mesh.nCells(2); k++) {
        for (intType j = 0; j != mesh.nCells(1); j++) {
            for (intType i = 0; i != mesh.nCells(0); i++) {

                momentumEquation.AU[axis][p]( G(i, j, k) ) += 1.5f * dtInv;
                momentumEquation.B( G(i, j, k) )           += (
                                                              - 2.0f * fieldsPrevTime.U[axis]( G(i, j, k) )
                                                              + 0.5f * fieldsPrevPrevTime.U[axis]( G(i, j, k) )
                                                              ) * dtInv;

            }
        }
    }
}




/*---------------------------------------------------------------------------------------------------------------*\
                                        Boundary Constants to Source Term
\*---------------------------------------------------------------------------------------------------------------*/

void AddMomentumBoundaryConstants( MomentumEquation &momCoeffs,
                                   const Mesh &mesh )
{
    const TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                        extents = {mesh.nCells[0], mesh.nCells[1], mesh.nCells[2]};

    // Each axis
    for (intType axis = 0; axis != Axis::count; axis++) {

        const BoundaryPatches::ENUMDATA positivePatch = LUT::PositivePatch[ static_cast<size_t>( axis ) ],
                                        negativePatch = LUT::NegativePatch[ static_cast<size_t>( axis ) ];
        const intType iEnd = mesh.nCells[axis] - 1;

        // Negative side boundary
        if ( momCoeffs.BUBoundary[negativePatch].size() != 0 )
            momCoeffs.B.slice(offsets, extents).chip( 0   , axis ) += momCoeffs.BUBoundary[negativePatch];

        if ( momCoeffs.BPBoundary[negativePatch].size() != 0 )
            momCoeffs.B.slice(offsets, extents).chip( 0   , axis ) += momCoeffs.BPBoundary[negativePatch];


        // Positive side boundary
        if ( momCoeffs.BUBoundary[positivePatch].size() != 0 )
            momCoeffs.B.slice(offsets, extents).chip( iEnd, axis ) += momCoeffs.BUBoundary[positivePatch];

        if ( momCoeffs.BPBoundary[positivePatch].size() != 0 )
            momCoeffs.B.slice(offsets, extents).chip( iEnd, axis ) += momCoeffs.BPBoundary[positivePatch];

    }
}


void AddContinuityBoundaryConstants( ContinuityEquation &contCoeffs,
                                     const Mesh &mesh )
{
    const TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                        extents = {mesh.nCells[0], mesh.nCells[1], mesh.nCells[2]};

    // Each axis
    for (intType axis = 0; axis != Axis::count; axis++) {

        const BoundaryPatches::ENUMDATA positivePatch = LUT::PositivePatch[ static_cast<size_t>( axis ) ],
                                        negativePatch = LUT::NegativePatch[ static_cast<size_t>( axis ) ];
        const intType iEnd = mesh.nCells[axis] - 1;

        // Negative side boundary
        if ( contCoeffs.BUBoundary[negativePatch].size() != 0 )
            contCoeffs.B.slice(offsets, extents).chip( 0   , axis ) += contCoeffs.BUBoundary[negativePatch];

        if ( contCoeffs.BPBoundary[negativePatch].size() != 0 )
            contCoeffs.B.slice(offsets, extents).chip( 0   , axis ) += contCoeffs.BPBoundary[negativePatch];


        // Positive side boundary
        if ( contCoeffs.BUBoundary[positivePatch].size() != 0 )
            contCoeffs.B.slice(offsets, extents).chip( iEnd, axis ) += contCoeffs.BUBoundary[positivePatch];

        if ( contCoeffs.BPBoundary[positivePatch].size() != 0 )
            contCoeffs.B.slice(offsets, extents).chip( iEnd, axis ) += contCoeffs.BPBoundary[positivePatch];

    }
}



/*---------------------------------------------------------------------------------------------------------------*\
                                            Immersed Boundary Functions
\*---------------------------------------------------------------------------------------------------------------*/


// Set the IB source terms that come from the implicit stencil 
void MomentumIBSourceStencilPicard( MomentumEquation &momentumEquation,
                                    const IBCell::SourceTermData &sourceTermData, 
                                    const TensorIndex3D &cellIndex )
{
    using FVT::G;

    const Axis::ENUMDATA faceNormal = sourceTermData.direction;
    const Axis::ENUMDATA momentumAxis = momentumEquation.component;
    const TransportCoefficients::ENUMDATA coeff = ( sourceTermData.directionIndex == +1 ) ?  LUT::HiCoeff[faceNormal] : LUT::LoCoeff[faceNormal];

    // Velocity term
    floatType ibSource = momentumEquation.AU[momentumAxis][coeff]( G(cellIndex) ) * sourceTermData.ghostCellValues.U[momentumAxis];

    // Pressure stencil
    if ( momentumAxis == faceNormal ) {
        ibSource += momentumEquation.AP[coeff]( G(cellIndex[faceNormal]) ) * sourceTermData.ghostCellValues.P;
    }

    momentumEquation.B( G(cellIndex) ) += ibSource;
}



void MomentumIBSourceStencilNewton( MomentumEquation &momentumEquation,
                                    const IBCell::SourceTermData &sourceTermData, 
                                    const TensorIndex3D &cellIndex )
{
    using FVT::G;

    const Axis::ENUMDATA faceNormal = sourceTermData.direction;
    const Axis::ENUMDATA momentumAxis = momentumEquation.component;
    const TransportCoefficients::ENUMDATA coeff = ( sourceTermData.directionIndex == +1 ) ?  LUT::HiCoeff[faceNormal] : LUT::LoCoeff[faceNormal];
    
    // Velocity term
    floatType ibSource = momentumEquation.AU[momentumAxis][coeff]( G(cellIndex) ) * sourceTermData.ghostCellValues.U[momentumAxis];

    const Axis::ENUMDATA loAxis = LUT::LoOrthogonalAxis[ momentumAxis ];
    if ( faceNormal == loAxis ) {
        ibSource += momentumEquation.AU[loAxis][coeff]( G(cellIndex) ) * sourceTermData.ghostCellValues.U[loAxis];
    }

    const Axis::ENUMDATA hiAxis = LUT::HiOrthogonalAxis[ momentumAxis ];
    if ( faceNormal == hiAxis ) {
        ibSource += momentumEquation.AU[hiAxis][coeff]( G(cellIndex) ) * sourceTermData.ghostCellValues.U[hiAxis];
    }


    // Pressure term
    if ( momentumAxis == faceNormal ) {
        ibSource += momentumEquation.AP[coeff]( G(cellIndex[faceNormal]) ) * sourceTermData.ghostCellValues.P;
    }

    momentumEquation.B( G(cellIndex) ) += ibSource;
}



template< AdvectionSchemes advectionScheme >
void MomentumIBSourceDeferredCorrection( MomentumEquation &momentumEquation,
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

    // Advected velocities to remove
    floatType highOrderAdvectedVelocity{0.0f}, upwindAdvectedVelocity{0.0f};
    Tensor3D const &U = fields.U[ momentumEquation.component ];
    TensorIndex3D loIndex = faceIndex,
                  hiIndex = loIndex;
    hiIndex[faceNormal] += 1;
    if ( faceFluxes[faceNormal]( faceIndex ) >= 0.0f ) {
        highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme, +1>(U, momentumEquation, mesh, faceNormal, hiIndex, loIndex);
        upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, +1>(U, momentumEquation, mesh, faceNormal, hiIndex, loIndex);
    } else {
        highOrderAdvectedVelocity = FaceInterpolatedVelocity<advectionScheme, -1>(U, momentumEquation, mesh, faceNormal, hiIndex, loIndex);
        upwindAdvectedVelocity    = FaceInterpolatedVelocity<AdvectionSchemes::Upwind, -1>(U, momentumEquation, mesh, faceNormal, hiIndex, loIndex);
    }

    // Remove the deferred correction term
    floatType ibSource = static_cast<floatType>( sourceTermData.directionIndex ) 
                       * momentumEquation.advectionBlendingFactor 
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

        const floatType ghostCellValue = sourceTermData.ghostCellValues.U[momentumEquation.component];

        floatType oldHighOrderAdvectedVelocity{0.0f}, correctedHighOrderAdvectedVelocity{0.0f};
        if ( faceFluxes[faceNormal]( faceIndex ) >= 0.0f ) {
            oldHighOrderAdvectedVelocity       = FaceInterpolatedVelocity<advectionScheme, +1>(U, momentumEquation, mesh, faceNormal, hiIndex_a, loIndex_a);
            correctedHighOrderAdvectedVelocity = ( sourceTermData.directionIndex == +1 ) ? FaceInterpolatedVelocity<advectionScheme, +1, +1>(U, momentumEquation, mesh, faceNormal, hiIndex_a, loIndex_a, ghostCellValue):
                                                                                           FaceInterpolatedVelocity<advectionScheme, +1, -1>(U, momentumEquation, mesh, faceNormal, hiIndex_a, loIndex_a, ghostCellValue);
            
        } else {
            oldHighOrderAdvectedVelocity       = FaceInterpolatedVelocity<advectionScheme, -1>(U, momentumEquation, mesh, faceNormal, hiIndex_a, loIndex_a);
            correctedHighOrderAdvectedVelocity = ( sourceTermData.directionIndex == +1 ) ? FaceInterpolatedVelocity<advectionScheme, -1, +1>(U, momentumEquation, mesh, faceNormal, hiIndex_a, loIndex_a, ghostCellValue):
                                                                                           FaceInterpolatedVelocity<advectionScheme, -1, -1>(U, momentumEquation, mesh, faceNormal, hiIndex_a, loIndex_a, ghostCellValue);
        }


        // Boundary cell
        ibSource -= - static_cast<floatType>( sourceTermData.directionIndex ) 
                  *   momentumEquation.advectionBlendingFactor 
                  *   faceFluxes[faceNormal]( faceIndex_a )
                  *   ( correctedHighOrderAdvectedVelocity - oldHighOrderAdvectedVelocity )
                  *   mesh.cellLengthsInv[faceNormal]( cellIndex[faceNormal] );

        // Interior cell
        const floatType ibSource_a = - static_cast<floatType>( sourceTermData.directionIndex )
                                     *   momentumEquation.advectionBlendingFactor 
                                     *   faceFluxes[faceNormal]( faceIndex_a )
                                     *   ( correctedHighOrderAdvectedVelocity - oldHighOrderAdvectedVelocity )
                                     *   mesh.cellLengthsInv[faceNormal]( cellIndex_a[faceNormal] );

        momentumEquation.B( G(cellIndex_a) ) += ibSource_a;

    }

    momentumEquation.B( G(cellIndex) ) += ibSource;
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
    // TransportCoefficients::ENUMDATA coeff  = ( sourceTermData.directionIndex == +1 ) ? LUT::HiCoeff[faceNormal]   : LUT::LoCoeff[faceNormal];
    const TransportCoefficients::ENUMDATA ccoeff = ( sourceTermData.directionIndex == +1 ) ? LUT::HiHiCoeff[faceNormal] : LUT::LoLoCoeff[faceNormal];

    // The immediate cell
    // fvCoeffs.Cont.AP[coeff ](G(cellIndex)) = 0.0f;
    fvCoeffs.Cont.AP[ccoeff](G(cellIndex)) = 0.0f;

    // The interior cell
    // fvCoeffs.Cont.AP[ccoeff](G(sourceTermData.cellIndex_a)) = 0.0f;

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
    const Tensor3D &momentumDiagCoeffInv = fvCoeffs.Mom[faceNormal].diagCoeffInv;
    const std::array<Tensor1D, 4> &mwiSparseCoeffs = fvCoeffs.Cont.mwiSparseCoeffs[faceNormal];
    TensorIndex3D loIndex = ghostIsHiSide ? cellIndex : sourceTermData.cellIndex_g;
    TensorIndex3D hiIndex = ghostIsHiSide ? sourceTermData.cellIndex_g : cellIndex;
    intType idx = hiIndex[faceNormal];
    floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeffInv, mesh, faceNormal );

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
    d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeffInv, mesh, faceNormal );

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
    const Tensor3D &momentumDiagCoeffInv = fvCoeffs.Mom[faceNormal].diagCoeffInv;
    const std::array<Tensor1D, 4> &mwiSparseCoeffs = fvCoeffs.Cont.mwiSparseCoeffs[faceNormal];
    const floatType d = MWIWeightingCoeff( loIndex, hiIndex, momentumDiagCoeffInv, mesh, faceNormal );

    const floatType ghostSparseCoeff  = ghostIsHiSide ? d * mwiSparseCoeffs[3](idx) : - d * mwiSparseCoeffs[0](idx);

    const floatType ibSource = ghostSparseCoeff * sourceTermData.ghostCellValues.P * mesh.cellLengthsInv[faceNormal]( cellIndex[faceNormal] );

    fvCoeffs.Cont.B( G(sourceTermData.cellIndex_a) ) += ibSource;
}



void ChangeStencilToCentralAtIB( FVCoefficients &fvCoeffs,
                                 const EnumVector<Axis, Tensor3D> &faceFluxes, 
                                 const Mesh &mesh,
                                 const IBData &ibData )
{
    using enum TransportCoefficients::ENUMDATA;

    for ( auto &ibCell : ibData.ibCells ) { 

        TensorIndex3D cellIndex = ibCell.cellIndex;

        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            const Axis::ENUMDATA faceNormal = sourceTermData.direction;
            TensorIndex3D faceIndex = cellIndex;
            faceIndex[faceNormal] += sourceTermData.faceDirectionIndex;
            const intType fidx = faceIndex[faceNormal],
                          cidx = cellIndex[faceNormal];
            cellIndex = G(cellIndex);   // Coefficients have dummy cells

            const floatType faceFlux = faceFluxes[faceNormal](faceIndex);

            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                if ( sourceTermData.directionIndex == +1 ) {    // Face on Hi side

                    const TransportCoefficients::ENUMDATA hi = LUT::HiCoeff[faceNormal];

                    // Subtract upwinding term
                    if ( faceFlux >= 0.0f ) {
                        fvCoeffs.Mom[axis].AU[axis][p ](cellIndex) -= faceFlux * mesh.cellLengthsInv[faceNormal]( cidx );
                    } else {
                        fvCoeffs.Mom[axis].AU[axis][hi](cellIndex) -= faceFlux * mesh.cellLengthsInv[faceNormal]( cidx );
                    }

                    // Add in central differencing term
                    fvCoeffs.Mom[axis].AU[axis][p ](cellIndex) += faceFlux * ( 1.0f - mesh.interpFactors[faceNormal]( fidx ) ) * mesh.cellLengthsInv[faceNormal]( cidx );
                    fvCoeffs.Mom[axis].AU[axis][hi](cellIndex) += faceFlux * mesh.interpFactors[faceNormal]( fidx ) * mesh.cellLengthsInv[faceNormal]( cidx );

                } else {                                        // Face on Lo side    

                    const TransportCoefficients::ENUMDATA lo = LUT::LoCoeff[faceNormal];

                    // Subtract upwinding term
                    if ( faceFlux >= 0.0f ) {
                        fvCoeffs.Mom[axis].AU[axis][lo](cellIndex) -= - faceFlux * mesh.cellLengthsInv[faceNormal]( cidx );
                    } else {
                        fvCoeffs.Mom[axis].AU[axis][p ](cellIndex) -= - faceFlux * mesh.cellLengthsInv[faceNormal]( cidx );
                    }

                    // Add in central differencing term
                    fvCoeffs.Mom[axis].AU[axis][p ](cellIndex) += - faceFlux * mesh.interpFactors[faceNormal]( fidx ) * mesh.cellLengthsInv[faceNormal]( cidx );
                    fvCoeffs.Mom[axis].AU[axis][lo](cellIndex) += - faceFlux * ( 1.0f - mesh.interpFactors[faceNormal]( fidx ) ) * mesh.cellLengthsInv[faceNormal]( cidx );

                }

            } );
            
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
    for ( auto &ibCell : ibData.ibCells ) { 

        const TensorIndex3D cellIndex = ibCell.cellIndex;

        // A source term is added for each forced face
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            // Momentum equations
            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                switch ( fvCoeffs.Mom[axis].linearisation ) {
                    case Linearisation::Picard:
                        MomentumIBSourceStencilPicard( fvCoeffs.Mom[axis], sourceTermData, cellIndex );
                        break;

                    case Linearisation::Newton:
                        MomentumIBSourceStencilNewton( fvCoeffs.Mom[axis], sourceTermData, cellIndex );
                        break;
                }

                switch( fvCoeffs.Mom[axis].advectionScheme ) {
                    case AdvectionSchemes::Upwind:
                        /* NULL */
                        break;
                    case AdvectionSchemes::Central:
                        MomentumIBSourceDeferredCorrection<AdvectionSchemes::Central>( fvCoeffs.Mom[axis], fields, faceFluxes, mesh, sourceTermData, cellIndex );
                        break;
                    case AdvectionSchemes::SOU:
                        MomentumIBSourceDeferredCorrection<AdvectionSchemes::SOU>( fvCoeffs.Mom[axis], fields, faceFluxes, mesh, sourceTermData, cellIndex );
                        break;
                    case AdvectionSchemes::QUICK:
                        MomentumIBSourceDeferredCorrection<AdvectionSchemes::QUICK>( fvCoeffs.Mom[axis], fields, faceFluxes, mesh, sourceTermData, cellIndex );
                        break;
                }
               
            } );

            // Continuity equation
            switch ( fvCoeffs.Cont.momentumInterpolation ) {
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


    // Set in boundary coefficients for the continity equation to zero
    for ( auto &ibCell : ibData.ibCells ) { 

        const TensorIndex3D cellIndex = ibCell.cellIndex;
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            ZeroInSolidStencilCoeffs( fvCoeffs, sourceTermData, cellIndex );

        }

    }

}



/*---------------------------------------------------------------------------------------------------------------*\
                                                General Functions
\*---------------------------------------------------------------------------------------------------------------*/

// Allocate the 2D arrays which store constant values of boundary conditions if they are needed
void AllocateBoundaryConstants( FVCoefficients &fvCoeffs, 
                                const BoundaryConditionData &bcData )
{
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {

        // Dimensions of the patch
        const Axis::ENUMDATA patchAxis = LUT::BoundaryPatchAxis[ bp ];
        const intType patchDimLo = fvCoeffs.nCells( LUT::LoOrthogonalAxis[ patchAxis ] ),
                      patchDimHi = fvCoeffs.nCells( LUT::HiOrthogonalAxis[ patchAxis ] );


        // Velocity boundary conditions
        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

            if ( bcData.fields.U[axis][bp].type == BoundaryConditions::fixed ) {

                if ( patchAxis == axis ) {
                        if ( fvCoeffs.Mom[ LUT::LoOrthogonalAxis[patchAxis] ].linearisation == Linearisation::Newton ) 
                            fvCoeffs.Mom[ LUT::LoOrthogonalAxis[patchAxis] ].BUBoundary[bp] = Tensor2D(patchDimLo, patchDimHi).setZero();

                        if ( fvCoeffs.Mom[ LUT::HiOrthogonalAxis[patchAxis] ].linearisation == Linearisation::Newton ) 
                            fvCoeffs.Mom[ LUT::HiOrthogonalAxis[patchAxis] ].BUBoundary[bp] = Tensor2D(patchDimLo, patchDimHi).setZero();
                }

                fvCoeffs.Mom[axis].BUBoundary[bp] = Tensor2D(patchDimLo, patchDimHi).setZero();
                fvCoeffs.Cont.BUBoundary[bp] = Tensor2D(patchDimLo, patchDimHi).setZero();

            }

        } );


        // Pressure boundary conditions
        if ( bcData.fields.P[bp].type == BoundaryConditions::fixed ) {

            EnumFor<Axis>( [&] (Axis::ENUMDATA momAxis) {
                fvCoeffs.Mom[momAxis].BPBoundary[bp] = Tensor2D(patchDimLo, patchDimHi).setZero();
            } );
            fvCoeffs.Cont.BPBoundary[bp] = Tensor2D(patchDimLo, patchDimHi).setZero();

        }

    } );

}



// Set the coefficients that need to be relinearised to zero
void ZeroNonlinearCoeffs( FVCoefficients &fvCoeffs )
{
    using enum TransportCoefficients::ENUMDATA;
    using enum Axis::ENUMDATA;
    
    // Zero momentum equations
    EnumFor<Axis> ( [&] (Axis::ENUMDATA axis) {

        EnumFor<TransportCoefficients>( [&] (TransportCoefficients::ENUMDATA tc) {
            fvCoeffs.Mom[axis].AU[axis][tc].setZero();
        } );
        // fvCoeffs.Mom[axis].AU[axis][p].setZero();

        if ( fvCoeffs.Mom[axis].linearisation == Linearisation::Newton ) {

            if ( axis != X ) {
                fvCoeffs.Mom[axis].AU[X][e].setZero();
                fvCoeffs.Mom[axis].AU[X][p].setZero();
                fvCoeffs.Mom[axis].AU[X][w].setZero();
            }
            
            if ( axis != Y ) {
                fvCoeffs.Mom[axis].AU[Y][n].setZero();
                fvCoeffs.Mom[axis].AU[Y][p].setZero();
                fvCoeffs.Mom[axis].AU[Y][s].setZero();
            }

            if ( axis != Z ) {
                fvCoeffs.Mom[axis].AU[Z][t].setZero();
                fvCoeffs.Mom[axis].AU[Z][p].setZero();
                fvCoeffs.Mom[axis].AU[Z][b].setZero();
            }
            
        } 

        EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
            fvCoeffs.Mom[axis].BUBoundary[bp].setZero();
            fvCoeffs.Mom[axis].BPBoundary[bp].setZero();
        } );

        fvCoeffs.Mom[axis].B.setZero();
        fvCoeffs.Mom[axis].F.setZero();

    } );

    
    // Zero continuity equation
    EnumFor<TransportCoefficients>( [&] (TransportCoefficients::ENUMDATA tc) {
            fvCoeffs.Cont.AP[tc].setZero();
        } );
    // fvCoeffs.Cont.AP[p].setZero();

    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        fvCoeffs.Cont.BUBoundary[bp].setZero();
        fvCoeffs.Cont.BPBoundary[bp].setZero();
    } );

    fvCoeffs.Cont.B.setZero();
    fvCoeffs.Cont.F.setZero();
}


void SetDiagCoeffInverse( MomentumEquation &momentumEquation,
                          const Mesh &mesh )
{
    using TC = TransportCoefficients;
    const TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                        extents = {mesh.nCells[0], mesh.nCells[1], mesh.nCells[2]};

    const Axis::ENUMDATA axis = momentumEquation.component;
    momentumEquation.diagCoeffInv.slice(offsets, extents) = momentumEquation.AU[axis][TC::p].slice(offsets, extents).inverse();

}


}   // end anonymous namespace





/*---------------------------------------------------------------------------------------------------------------*\
                                            Set and Update Functions
\*---------------------------------------------------------------------------------------------------------------*/

FVCoefficients InitialiseFVCoefficients( const Mesh &mesh,
                                         const FieldData<Tensor3D> &fields,
                                         const FieldData<Tensor3D> &fieldsPrevTime,
                                         const FieldData<Tensor3D> &fieldsPrevPrevTime,
                                         const EnumVector< Axis, EnumVector< Axis, Tensor3D> > &faceAdvectedVelocities,
                                         const EnumVector<Axis, Tensor3D> &faceFluxes, 
                                         const IBData &ibData,
                                         const BoundaryConditionData &bcData,
                                         const InputData &inputData)
{
    // Default construct the coefficients class
    FVCoefficients fvCoeffs(mesh.nCells, inputData.schemes.linearisation, inputData.schemes.momentumInterpolation);

    fvCoeffs.rho = inputData.rho;
    fvCoeffs.nu  = inputData.nu;

    // Allocate boundary constant arrays for fixed boundary conditions
    AllocateBoundaryConstants( fvCoeffs, bcData );

    // Parts of momentum equation that don't change with linearisation
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        fvCoeffs.Mom[axis].advectionScheme         = inputData.schemes.advectionScheme;
        fvCoeffs.Mom[axis].timeScheme              = inputData.schemes.timeScheme;
        fvCoeffs.Mom[axis].timeStep                = inputData.schemes.timeStep;
        fvCoeffs.Mom[axis].advectionBlendingFactor = inputData.schemes.advectionBlendingFactor;

        SetDiffusionCoeffients(fvCoeffs.Mom[axis], bcData, inputData.nu, mesh);
        SetMomentumPressureCoefficients(fvCoeffs.Mom[axis], inputData.rho, mesh, bcData.fields.P);
        SetHighOrderAdvectionCoefficients(fvCoeffs.Mom[axis], mesh);
    } );


    // Parts of continuity equation that don't change with linearisation
    SetContinuityVelocityCoefficients(fvCoeffs.Cont, mesh, bcData.fields.U);
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        SetMomentumInterpolationSparseConstants(fvCoeffs.Cont.mwiSparseCoeffs[axis], fvCoeffs.Mom[axis].AP, mesh, axis);
        SetMomentumInterpolationCompactConstants(fvCoeffs.Cont.mwiCompactCoeffs[axis], inputData.rho, mesh, axis);
    } );

    // Set the coefficients that depend on linearisation
    UpdateFVCoefficients( fvCoeffs, mesh, fields, fieldsPrevTime, fieldsPrevPrevTime, faceAdvectedVelocities, faceFluxes, ibData, bcData );

    return fvCoeffs;
}



// Update linearisation in momenum and continuity equations
void UpdateFVCoefficients( FVCoefficients &fvCoeffs, 
                           const Mesh &mesh,
                           const FieldData<Tensor3D> &fields,
                           const FieldData<Tensor3D> &fieldsPrevTime,
                           const FieldData<Tensor3D> &fieldsPrevPrevTime,
                           const EnumVector< Axis, EnumVector< Axis, Tensor3D> > &faceAdvectedVelocities,
                           const EnumVector<Axis, Tensor3D> &faceFluxes,
                           const IBData &ibData,
                           const BoundaryConditionData &bcData )
{
    ZeroNonlinearCoeffs( fvCoeffs );

    // The implicit Picard coefficients for all momentum equations are the same, so just use the ones from the U momentum 
    // equation after its been set
    if ( fvCoeffs.Mom[Axis::X].advectionScheme == AdvectionSchemes::Upwind ) {
        SetInteriorAdvectionPicardCoefficients(fvCoeffs.Mom[Axis::X], fields, faceFluxes, mesh);
        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
            if ( axis != Axis::X ) {
                EnumFor<TransportCoefficients>( [&] (TransportCoefficients::ENUMDATA tc) {
                    fvCoeffs.Mom[axis].AU[axis][tc] = fvCoeffs.Mom[Axis::X].AU[Axis::X][tc];      
                } );
            }
        } );
    } else {    // Source terms are different due to deferred correction
        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
            SetInteriorAdvectionPicardCoefficients(fvCoeffs.Mom[axis], fields, faceFluxes, mesh);
        } );
    }

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        // Boundaries need to be done after since they can affect the internal coefficients
        SetBoundaryAdvectionPicardCoefficients(fvCoeffs.Mom[axis], faceFluxes, bcData.fields.U[axis], mesh);
        AddMomentumPressureBoundarySources(fvCoeffs.Mom[axis], fvCoeffs.rho, mesh, bcData.fields.P);
        AddDiffusion(fvCoeffs.Mom[axis], bcData.fields.U[axis], mesh);
    } );

    // Change stencil in momentum equations to have central differencing at the boundaries
    ChangeStencilToCentralAtIB( fvCoeffs, faceFluxes, mesh, ibData );

    // Add unsteady term
    switch ( fvCoeffs.Mom[Axis::X].timeScheme ) {
        case TimeSchemes::BackwardsEuler:
            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                AddBackwardsEuler(fvCoeffs.Mom[axis], mesh, fieldsPrevTime);
            } );
            break;

        case TimeSchemes::BackwardsThreeLevel:
            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                AddBackwardsThreeLevel(fvCoeffs.Mom[axis], mesh, fieldsPrevTime, fieldsPrevPrevTime);
            } );
            break;

        case TimeSchemes::Steady:
            /* NULL */
            break;
    }

    // Inverse of AP coefficient (Picard)
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {     
        SetDiagCoeffInverse(fvCoeffs.Mom[axis], mesh);   
    } );

    SetMomentumInterpolationCoefficients( fvCoeffs, mesh, bcData, fields.P );
    AddContinuityVelocityBoundarySources( fvCoeffs.Cont, mesh, bcData.fields.U );

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        // Add Newton Linearisation terms if selected
        if ( fvCoeffs.Mom[axis].linearisation == Linearisation::Newton ) {
            AddAdvectionNewtonCoefficients(fvCoeffs.Mom[axis], faceAdvectedVelocities, faceFluxes, bcData.fields.U, mesh);
            SetDiagCoeffInverse(fvCoeffs.Mom[axis], mesh);
        }

        // Add boundary constants to source terms
        AddMomentumBoundaryConstants(fvCoeffs.Mom[axis], mesh);

    } );
    AddContinuityBoundaryConstants(fvCoeffs.Cont, mesh);


    // Add effect of immersed boundary
    AddIBSourceTerms( fvCoeffs, faceFluxes, ibData, fields, mesh );
}


}   // end namespace CFD
