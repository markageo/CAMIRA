#include "FiniteVolume.h"
#include "../Macros.h"
#include "../Tools/FVTools.h"
#include "../Tools/FVLookups.h"

#include <algorithm>
#include <iostream>
#include <vector>

namespace CFD
{

using namespace FVT;
 
namespace
{


/*---------------------------------------------------------------------------------------------------------------*\
                                                    Diffusion
\*---------------------------------------------------------------------------------------------------------------*/

// Check if continuity equation implies a zero gradient boundary condition. This occurs if both orthogonal fields have a uniform BC
BoundaryConditions::ENUMDATA GetDiffusionBC( const EnumVector< Axis, EnumVector< BoundaryPatches, BoundaryConditionConfig > > &MomBoundaryConditions, 
                                             const BoundaryPatches::ENUMDATA boundaryPatch, 
                                             const Axis::ENUMDATA velocityComponent )
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;

    const Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[boundaryPatch];

    // Set the field we need to check based on the axis
    Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[ axis ];
    Axis::ENUMDATA axis2 = LUT::LoOrthogonalAxis[ axis ];

    // Only check the field that in the direction of the current axis
    if (velocityComponent == axis) {

        if (MomBoundaryConditions[axis1][boundaryPatch].type == BC::fixed && 
            MomBoundaryConditions[axis2][boundaryPatch].type == BC::fixed) {
            return BC::zeroGradient;
        }
    }

    return MomBoundaryConditions[velocityComponent][boundaryPatch].type;
}


// Apply boundary conditions for diffusion terms on axis positive boundary
void DiffusionPositiveBoundary( EnumVector< Axis,  EnumVector<TransportCoefficients, array1D> > &diff, 
                                EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, BoundaryConditionConfig > &boundaryConditionStructs,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = LUT::LoCoeff[axis];
    const intType iCellBound = mesh.nCells(axis) - 1;

    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
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

        default:
            break;
    }

}


// Apply boundary conditions for diffusion terms on axis negative boundary
void DiffusionNegativeBoundary( EnumVector< Axis, EnumVector<TransportCoefficients, array1D> > &diff, 
                                EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, BoundaryConditionConfig > &boundaryConditionStructs,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::NegativePatch[axis];
    const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis];
    const intType iCellBound = 0;

    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
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

        default:
            break;
    }

}


// Set diffusion coefficients for a given momentum equation
void SetDiffusionCoeffients( MomentumEquation &momentumEquation, 
                             const FieldData< EnumVector< BoundaryPatches, BoundaryConditionConfig > > &bcData, 
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
        positivePatchBC = GetDiffusionBC(bcData.U, positivePatch, velocityComponent);
        negativePatchBC = GetDiffusionBC(bcData.U, negativePatch, velocityComponent); 


        // Boundary conditions only need to be set if it is not zero gradient
        if (positivePatchBC != BC::zeroGradient) {
            DiffusionPositiveBoundary(diff, diffBoundary, mesh, bcData.U[velocityComponent], axis);
        }

        if (negativePatchBC != BC::zeroGradient) {
            DiffusionNegativeBoundary(diff, diffBoundary, mesh, bcData.U[velocityComponent], axis);
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


// Upwind coefficients
// Assumes that the 'p' coefficient has been set to zero
[[maybe_unused]]
void UpwindInteriorPicard( EnumVector<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
                           const EnumVector<Axis, array3D> &faceFluxes, 
                           const Mesh &mesh,
                           const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis], 
                                    west = LUT::LoCoeff[axis];

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                arrayIndex3D HiIndex = { i, j, k },
                             LoIndex = { i, j, k };
                LoIndex[axis] -= 1;

                floatType uf = faceFluxes[ axis ](i, j, k);
                coeffs[p   ](LoIndex) += std::max(   uf * mesh.cellLengthsInv[axis]( LoIndex[axis] ), static_cast<floatType>(0.0f) );
                coeffs[east](LoIndex)  = std::min(   uf * mesh.cellLengthsInv[axis]( LoIndex[axis] ), static_cast<floatType>(0.0f) );

                coeffs[p   ](HiIndex) += std::max( - uf * mesh.cellLengthsInv[axis]( HiIndex[axis] ), static_cast<floatType>(0.0f) );
                coeffs[west](HiIndex)  = std::min( - uf * mesh.cellLengthsInv[axis]( HiIndex[axis] ), static_cast<floatType>(0.0f) );

            }
        }
    }
}



// Upwind coefficients
// Assumes that the 'p' coefficient has been set to zero
// Helps the compiler autovectorise
template< Axis::ENUMDATA axis> [[maybe_unused]]
void UpwindInteriorPicard_autoVec( EnumVector<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
                                   const EnumVector<Axis, array3D> &faceFluxes, 
                                   const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis], 
                                    west = LUT::LoCoeff[axis];

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {

            // Left side cells
            CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                arrayIndex3D  LoIndex = { i, j, k };
                LoIndex[axis] -= 1;

                floatType uf = faceFluxes[ axis ](i, j, k);
                coeffs[p   ](LoIndex) += std::max(   uf * mesh.cellLengthsInv[axis]( LoIndex[axis] ), static_cast<floatType>(0.0f) );
                coeffs[east](LoIndex)  = std::min(   uf * mesh.cellLengthsInv[axis]( LoIndex[axis] ), static_cast<floatType>(0.0f) );

            }


            // Right side cells
            CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                arrayIndex3D HiIndex = { i, j, k };

                floatType uf = faceFluxes[ axis ](i, j, k);
                coeffs[p   ](HiIndex) += std::max( - uf * mesh.cellLengthsInv[axis]( HiIndex[axis] ), static_cast<floatType>(0.0f) );
                coeffs[west](HiIndex)  = std::min( - uf * mesh.cellLengthsInv[axis]( HiIndex[axis] ), static_cast<floatType>(0.0f) );

            }
        }
    }
}



void AdvectionPositiveBoundary( EnumVector<TransportCoefficients, array3D> &coeffs, 
                                EnumVector<BoundaryPatches, array2D> &boundaryConstants,
                                const EnumVector<Axis, array3D> &laggedVelocity, 
                                const Mesh &mesh,  
                                const BoundaryConditionData &boundaryConditionData,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = LUT::LoCoeff[axis];
    const intType iCellBound = mesh.nCells(axis) - 1;   // Index of cell at the boundary
    const intType iFaceBound = iCellBound + 1;          // Index of face at the boundary

    switch ( boundaryConditionData[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p].chip(iCellBound, axis) += laggedVelocity[axis].chip(iFaceBound, axis) 
                                              * laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::fixed:
            boundaryConstants[boundaryPatch]  += laggedVelocity[axis].chip(iFaceBound, axis)
                                               * boundaryConditionData[boundaryPatch].value
                                               * laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::extrapolated:
            coeffs[p   ].chip(iCellBound, axis) += laggedVelocity[axis].chip(iFaceBound, axis) 
                                                 * laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].p * mesh.cellLengthsInv[axis](iCellBound) );
            coeffs[west].chip(iCellBound, axis) += laggedVelocity[axis].chip(iFaceBound, axis) 
                                                 * laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].a * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        default:
            break;
    }

}


void AdvectionNegativeBoundary( EnumVector<TransportCoefficients, array3D> &coeffs, 
                                EnumVector<BoundaryPatches, array2D> &boundaryConstants,
                                const EnumVector<Axis, CFD::array3D> &laggedVelocity, 
                                const Mesh &mesh,  
                                const BoundaryConditionData &boundaryConditionData,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::NegativePatch[axis];
    const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis];
    const intType iCellBound = 0;   // Index of cell at the boundary 
    const intType iFaceBound = 0;   // Index of face at the boundary

    switch ( boundaryConditionData[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p].chip(iCellBound, axis) += - laggedVelocity[axis].chip(iFaceBound, axis) 
                                              *   laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::fixed:
            boundaryConstants[boundaryPatch]  += - laggedVelocity[axis].chip(iFaceBound, axis)
                                               *   boundaryConditionData[boundaryPatch].value 
                                               *   laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::extrapolated:
            coeffs[p   ].chip(iCellBound, axis) += - laggedVelocity[axis].chip(iFaceBound, axis) 
                                                 *   laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].p * mesh.cellLengthsInv[axis](iCellBound) );
            coeffs[east].chip(iCellBound, axis) += - laggedVelocity[axis].chip(iFaceBound, axis) 
                                                 *   laggedVelocity[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].a * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        default:
            break;
    }

}



void SetInteriorAdvectionPicardCoefficients( MomentumEquation &momentumEquation,
                                             const EnumVector<Axis, array3D> &faceFluxes,
                                             const Mesh &mesh )
{
    using enum TransportCoefficients::ENUMDATA;
    
    auto &coeffs            = momentumEquation.AU[ momentumEquation.component ];

    #if defined( CFD_USE_AUTOVEC_FUNCTIONS )
        using enum Axis::ENUMDATA;
        UpwindInteriorPicard_autoVec<X>(coeffs, faceFluxes, mesh);
        UpwindInteriorPicard_autoVec<Y>(coeffs, faceFluxes, mesh);
        UpwindInteriorPicard_autoVec<Z>(coeffs, faceFluxes, mesh);
    #else
        EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {
            UpwindInteriorPicard(coeffs, faceFluxes, mesh, axis);
        } );
    #endif

}



void SetBoundaryAdvectionPicardCoefficients( MomentumEquation &momentumEquation,
                                             const EnumVector<Axis, array3D> &faceFluxes, 
                                             const BoundaryConditionData &boundaryConditions,
                                             const Mesh &mesh )
{
    using enum TransportCoefficients::ENUMDATA;
    
    auto &coeffs            = momentumEquation.AU[ momentumEquation.component ];
    auto &boundaryConstants = momentumEquation.BUBoundary;

    // Upwind internal faces
    EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {

        // Boundary faces
        AdvectionPositiveBoundary(coeffs, boundaryConstants, faceFluxes, mesh, boundaryConditions, axis);
        AdvectionNegativeBoundary(coeffs, boundaryConstants, faceFluxes, mesh, boundaryConditions, axis);

    } );

}





/*---------------------------------------------------------------------------------------------------------------*\
                                       Momentum Newton Advection Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

[[maybe_unused]]
void NewtonInteriorImplicit( EnumVector< TransportCoefficients, array3D > &coeffs, 
                             const EnumVector< Axis, array3D > &faceAdvectedVelocities,  
                             const Mesh &mesh,
                             const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis], 
                                    west = LUT::LoCoeff[axis];

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                arrayIndex3D HiIndex = { i, j, k },
                             LoIndex = { i, j, k };
                LoIndex[axis] -= 1;
                intType idx = HiIndex[axis];

                floatType coeffLo = faceAdvectedVelocities[axis](i, j, k) * mesh.cellLengthsInv[axis]( LoIndex[axis] );
                coeffs[p   ](LoIndex) += coeffLo * ( 1 - mesh.interpFactors[axis]( idx ) );
                coeffs[east](LoIndex) += coeffLo * mesh.interpFactors[axis]( idx );

                floatType coeffHi = faceAdvectedVelocities[axis](i, j, k) * mesh.cellLengthsInv[axis]( HiIndex[axis] );
                coeffs[p   ](HiIndex) += - coeffHi * mesh.interpFactors[axis]( idx );
                coeffs[west](HiIndex) += - coeffHi * ( 1 - mesh.interpFactors[axis]( idx ) );

            }
        }
    }
}



// Templated to allow compiler autovectorisation
template< Axis::ENUMDATA axis > [[maybe_unused]]
void NewtonInteriorImplicit_autoVec( EnumVector< TransportCoefficients, array3D > &coeffs, 
                                     const EnumVector< Axis, array3D > &faceAdvectedVelocities,  
                                     const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis], 
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
            CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                if constexpr ( axis == X ) {
                    LoCellLengthInv = mesh.cellLengthsInv[axis]( i-1 );
                    LoInterpFactor  = mesh.interpFactors[axis]( i-1 );
                }

                arrayIndex3D LoIndex = { i, j, k };
                LoIndex[axis] -= 1;
                
                floatType coeffLo = faceAdvectedVelocities[axis](i, j, k) * LoCellLengthInv;
                coeffs[p   ](LoIndex) += coeffLo * ( 1 - LoInterpFactor );
                coeffs[east](LoIndex) += coeffLo * LoInterpFactor;

            }

            // Right side cells
            CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                if constexpr ( axis == X ) {
                    HiCellLengthInv = mesh.cellLengthsInv[axis]( i );
                    HiInterpFactor  = mesh.interpFactors[axis]( i );
                }

                arrayIndex3D HiIndex = { i, j, k };

                floatType coeffHi = faceAdvectedVelocities[axis](i, j, k) * HiCellLengthInv;
                coeffs[p   ](HiIndex) += - coeffHi * HiInterpFactor;
                coeffs[west](HiIndex) += - coeffHi * ( 1 - HiInterpFactor );

            }
        }
    }

}




void NewtonConstants( array3D &B,
                      const EnumVector< Axis, array3D > &faceAdvectedVelocities,
                      const EnumVector< Axis, array3D > &faceFluxes,
                      const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {

            CFD_PRAGMA_VECTORIZE
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                floatType xFluxDiff = mesh.cellLengthsInv[X](i) 
                                    * ( faceFluxes[X](i+1, j, k) * faceAdvectedVelocities[X](i+1, j, k) 
                                      - faceFluxes[X](i  , j, k) * faceAdvectedVelocities[X](i  , j, k) );

                floatType yFluxDiff = mesh.cellLengthsInv[Y](j) 
                                    * ( faceFluxes[Y](i, j+1, k) * faceAdvectedVelocities[Y](i, j+1, k) 
                                      - faceFluxes[Y](i, j  , k) * faceAdvectedVelocities[Y](i, j  , k) );

                floatType zFluxDiff = mesh.cellLengthsInv[Z](k) 
                                    * ( faceFluxes[Z](i, j, k+1) * faceAdvectedVelocities[Z](i, j, k+1) 
                                      - faceFluxes[Z](i, j, k  ) * faceAdvectedVelocities[Z](i, j, k  ) );

                B( i, j, k ) += xFluxDiff + yFluxDiff + zFluxDiff;

            }
        }
    }

}



void AddAdvectionNewtonCoefficients( MomentumEquation &momentumEquation,
                                     const EnumVector< Axis, EnumVector<Axis, array3D> > &faceAdvectedVelocities, 
                                     const EnumVector< Axis, array3D> &faceFluxes,
                                     const EnumVector< Axis, BoundaryConditionData> &boundaryConditions,
                                     const Mesh &mesh )
{
    const auto &faceVelComp = faceAdvectedVelocities[momentumEquation.component];
    auto &boundaryConstants = momentumEquation.BUBoundary;
    
    // Implicit terms
    #if defined( CFD_USE_AUTOVEC_FUNCTIONS )
        using enum Axis::ENUMDATA;
        NewtonInteriorImplicit_autoVec<X>(momentumEquation.AU[X], faceVelComp, mesh );
        NewtonInteriorImplicit_autoVec<Y>(momentumEquation.AU[Y], faceVelComp, mesh );
        NewtonInteriorImplicit_autoVec<Z>(momentumEquation.AU[Z], faceVelComp, mesh );
    #else
        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
            NewtonInteriorImplicit(momentumEquation.AU[axis], faceVelComp, mesh, axis );
        } );
    #endif

    // Implicit boundary terms
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        auto &coeffs            = momentumEquation.AU[ axis ];
        AdvectionPositiveBoundary(coeffs, boundaryConstants, faceVelComp, mesh, boundaryConditions[axis], axis);
        AdvectionNegativeBoundary(coeffs, boundaryConstants, faceVelComp, mesh, boundaryConditions[axis], axis);
    } );

    // Explicit terms
    NewtonConstants( momentumEquation.B, faceVelComp, faceFluxes, mesh );
}





/*---------------------------------------------------------------------------------------------------------------*\
                                           Add Diffusion Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


void AddDiffusion( MomentumEquation &momentumEquation,
                   const EnumVector< BoundaryPatches, BoundaryConditionConfig > &bcStructs,
                   const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    auto &velCoeffs   = momentumEquation.AU[ momentumEquation.component ];
    auto &boundaryVel = momentumEquation.BUBoundary;

    const auto &diffCoeffs         = momentumEquation.diff;
    const auto &diffCoeffsBoundary = momentumEquation.diffBoundary;

    for (intType k = 0; k != mesh.nCells(Z); k++) {

        floatType zpk = diffCoeffs[Z][p](k),
                  ztk = diffCoeffs[Z][t](k),
                  zbk = diffCoeffs[Z][b](k);

        for (intType j = 0; j != mesh.nCells(Y); j++) {

            floatType ypj = diffCoeffs[Y][p](j),
                      ynj = diffCoeffs[Y][n](j),
                      ysj = diffCoeffs[Y][s](j);

            CFD_PRAGMA_VECTORIZE
            for (intType i = 0; i != mesh.nCells(X); i++) {

                velCoeffs[p](i, j, k) += diffCoeffs[X][p](i) + ypj + zpk;

                velCoeffs[e](i, j, k) += diffCoeffs[X][e](i);
                velCoeffs[w](i, j, k) += diffCoeffs[X][w](i);
                
                velCoeffs[n](i, j, k) += ynj;
                velCoeffs[s](i, j, k) += ysj;

                velCoeffs[t](i, j, k) += ztk;
                velCoeffs[b](i, j, k) += zbk;
            }
        }
    }


    // Constant terms
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA patch) {

        boundaryVel[patch] += boundaryVel[patch].constant( diffCoeffsBoundary[patch] )
                            * bcStructs[patch].value;
        
    } );

}





/*---------------------------------------------------------------------------------------------------------------*\
                                        Linear Interpolated Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


void InterpolationPositiveBoundary( EnumVector< TransportCoefficients, array1D > &coeffs, 
                                    EnumVector< BoundaryPatches, array2D > &boundaryConstants,
                                    const Mesh &mesh,  
                                    const EnumVector< BoundaryPatches, BoundaryConditionConfig > &bcStructs,
                                    const Axis::ENUMDATA axis)
{

    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = LUT::LoCoeff[axis];
    const intType iCellBound = mesh.nCells(axis) - 1;

    switch ( bcStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p   ]( iCellBound ) += 1;
            coeffs[west]( iCellBound ) += 0;
            break;

        case BC::fixed:
            coeffs[p   ]( iCellBound ) += 0;
            coeffs[west]( iCellBound ) += 0;
            boundaryConstants[boundaryPatch] = bcStructs[boundaryPatch].value;
            break;

        case BC::extrapolated:
            coeffs[p   ]( iCellBound ) += mesh.extrapFactors[boundaryPatch].p;
            coeffs[west]( iCellBound ) += mesh.extrapFactors[boundaryPatch].a;
            break;

        default:
            break;
    }

}


void InterpolationNegativeBoundary( EnumVector< TransportCoefficients, array1D > &coeffs, 
                                    EnumVector< BoundaryPatches, array2D > &boundaryConstants,
                                    const Mesh &mesh,  
                                    const EnumVector< BoundaryPatches, BoundaryConditionConfig > &bcStructs, 
                                    const Axis::ENUMDATA axis)
{

    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::NegativePatch[axis];
    const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis];
    const intType iCellBound = 0;

    switch ( bcStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p   ]( iCellBound ) += - 1;
            coeffs[east]( iCellBound ) +=   0;
            break;

        case BC::fixed:
            coeffs[p   ]( iCellBound ) += 0;
            coeffs[east]( iCellBound ) += 0;
            boundaryConstants[boundaryPatch] = - bcStructs[boundaryPatch].value;
            break;

        case BC::extrapolated:
            coeffs[p   ]( iCellBound ) += - mesh.extrapFactors[boundaryPatch].p;
            coeffs[east]( iCellBound ) += - mesh.extrapFactors[boundaryPatch].a;
            break;

        default:
            break;
    }

}


// Set coefficients for quantities that are intrpolated linearly onto faces.
void SetFaceInterpolatedCoefficients( EnumVector<TransportCoefficients, array1D> &coeffs, 
                                      EnumVector< BoundaryPatches, array2D > &boundaryConstants, 
                                      const Mesh &mesh, 
                                      const EnumVector< BoundaryPatches, BoundaryConditionConfig> &bcStructs, 
                                      const Axis::ENUMDATA axis)
{ 
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis],    // These are just names, they can be north, south etc.
                                    west = LUT::LoCoeff[axis];  

    // Internal faces
    for (intType i = 1; i != mesh.nCells(axis); i++) {
        
        // Cell on west side
        coeffs[p   ](i-1) += 1 - mesh.interpFactors[axis](i);
        coeffs[east](i-1) += mesh.interpFactors[axis](i);

        // Cell on east side
        coeffs[p   ](i) += - mesh.interpFactors[axis](i);
        coeffs[west](i) += - ( 1 - mesh.interpFactors[axis](i) ); 

    }

    // Boundary faces
    InterpolationPositiveBoundary(coeffs, boundaryConstants, mesh, bcStructs, axis);
    InterpolationNegativeBoundary(coeffs, boundaryConstants, mesh, bcStructs, axis);
    
    // Multiply by inverse of cell length
    for (intType i = 0; i != mesh.nCells(axis); i++) {
        coeffs[p   ](i) *= mesh.cellLengthsInv[axis](i);
        coeffs[east](i) *= mesh.cellLengthsInv[axis](i);
        coeffs[west](i) *= mesh.cellLengthsInv[axis](i); 
    }
    boundaryConstants[ LUT::PositivePatch[axis] ] *= boundaryConstants[ LUT::PositivePatch[axis] ].constant( mesh.cellLengthsInv[axis]( mesh.nCells(axis)-1 ) );
    boundaryConstants[ LUT::NegativePatch[axis] ] *= boundaryConstants[ LUT::NegativePatch[axis] ].constant( mesh.cellLengthsInv[axis]( 0 ) );
}



// Divide pressure terms in momentum equaqtions by density
void DivideMomentumPressureByDensity( MomentumEquation &momentumEquation, 
                                      const floatType rho)
{ 
    using enum TransportCoefficients::ENUMDATA;

    Axis::ENUMDATA axis = momentumEquation.component;

    auto &coeffs            = momentumEquation.AP;
    auto &boundaryConstants = momentumEquation.BPBoundary;

    TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis],    // These are just names, they can be north, south etc.
                                    west = LUT::LoCoeff[axis];  

    coeffs[p   ] /= coeffs[p   ].constant( rho );
    coeffs[east] /= coeffs[west].constant( rho );
    coeffs[west] /= coeffs[west].constant( rho );

    boundaryConstants[ LUT::PositivePatch[axis] ] /= boundaryConstants[ LUT::PositivePatch[axis] ].constant( rho );
    boundaryConstants[ LUT::NegativePatch[axis] ] /= boundaryConstants[ LUT::NegativePatch[axis] ].constant( rho );
}





/*---------------------------------------------------------------------------------------------------------------*\
                        Momentum Weighted Interpolation (Rhie-Chow Interpolation) Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

// Unweighted constants that appear in the sparse gradient part of MWI
void SetMomentumInterpolationSparseConstants( std::array< array1D, 4 > &mwiSparseCoeffs,
                                              const EnumVector<TransportCoefficients, array1D> &momentumPressureCoeffs,
                                              const Mesh &mesh,
                                              const Axis::ENUMDATA axis )
{
    using enum TransportCoefficients::ENUMDATA;

    TransportCoefficients::ENUMDATA west = LUT::LoCoeff[axis],
                                    east = LUT::HiCoeff[axis];

    // Internal faces
    for ( intType i = 1; i != mesh.cellFaces[axis].size()-1; i++ ) {

        mwiSparseCoeffs[0](i) = (1 - mesh.interpFactors[axis](i))   * momentumPressureCoeffs[west](i-1);

        mwiSparseCoeffs[1](i) = ( (1 - mesh.interpFactors[axis](i)) * momentumPressureCoeffs[p   ](i-1)
                              +    mesh.interpFactors[axis](i)      * momentumPressureCoeffs[west](i  ) );

        mwiSparseCoeffs[2](i) = ( (1 - mesh.interpFactors[axis](i)) * momentumPressureCoeffs[east](i-1)
                              +   mesh.interpFactors[axis](i)       * momentumPressureCoeffs[p   ](i  ) );

        mwiSparseCoeffs[3](i) = mesh.interpFactors[axis](i) * momentumPressureCoeffs[east](i);
    }
}



// Unweighted constants that appear in the compact gradient part of MWI
void SetMomentumInterpolationCompactConstants( std::array< array1D, 2 > &mwiCompactCoeffs,
                                               const floatType rho,
                                               const Mesh &mesh,
                                               const Axis::ENUMDATA axis )
{
    using enum TransportCoefficients::ENUMDATA;

    floatType rhoInv = 1 / rho;

    // Internal faces
    for ( intType i = 1; i != mesh.cellFaces[axis].size()-1; i++ ) {

        mwiCompactCoeffs[0](i) =   mesh.cellCenterDiffInv[axis](i) * rhoInv;

        mwiCompactCoeffs[1](i) = - mesh.cellCenterDiffInv[axis](i) * rhoInv;
    }
}



// Cell weighting coefficient for MWI. The given idx is the face index
floatType MWIWeightingCoeff( const arrayIndex3D &idx, 
                             const array3D &AUUpInv, 
                             const Axis::ENUMDATA axis)
{
    arrayIndex3D HiIndex = idx,
                 LoIndex = idx;
    LoIndex[axis] -= 1;
    return 0.5f * ( AUUpInv( HiIndex )  +  AUUpInv( LoIndex ) ); 
}



// Fully implicit momentum interpolation coefficient for internal faces
[[maybe_unused]]
void MWInterpolationInteriorImplicit( ContinuityEquation &continuityEquation,
                                      const array3D &momentumDiagCoeffInv,
                                      const Mesh &mesh,
                                      const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    // Unpack
    EnumVector<TransportCoefficients, array3D> &continuityPressureCoeffs = continuityEquation.AP;
    const std::array< array1D, 4 > &mwiSparseCoeffs                      = continuityEquation.mwiSparseCoeffs[axis];
    const std::array< array1D, 2 > &mwiCompactCoeffs                     = continuityEquation.mwiCompactCoeffs[axis];

    // Cell indexing
    TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                    eeast = LUT::HiHiCoeff[axis],
                                    west  = LUT::LoCoeff[axis],
                                    wwest = LUT::LoLoCoeff[axis];

    // Set the first most plane to zero only. We do this so that the high coefficients can be set in place, and less coefficients have to 
    // be zeroed upon re-linearisation.
    continuityPressureCoeffs[east ].chip(0, axis) = continuityPressureCoeffs[west ].chip(0, axis).constant( 0.0f );
    continuityPressureCoeffs[west ].chip(0, axis) = continuityPressureCoeffs[west ].chip(0, axis).constant( 0.0f );

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D HiIndex = { i, j, k },
                             LoIndex = { i, j, k };
                LoIndex[axis] -= 1;

                floatType d =  0.5f * ( momentumDiagCoeffInv( HiIndex )  +  momentumDiagCoeffInv( LoIndex ) );   

                // Coefficients for westmost to eastmost cell
                intType idx = HiIndex[axis];
                floatType coeff0 = d * mwiSparseCoeffs[0](idx),
                          coeff1 = d * ( mwiSparseCoeffs[1](idx) + mwiCompactCoeffs[0](idx) ),
                          coeff2 = d * ( mwiSparseCoeffs[2](idx) + mwiCompactCoeffs[1](idx) ),
                          coeff3 = d * mwiSparseCoeffs[3](idx);

                // Cell on west side 
                floatType LoCellLengthInv = mesh.cellLengthsInv[axis]( LoIndex[axis] );
                continuityPressureCoeffs[west ](LoIndex) += coeff0 * LoCellLengthInv;
                continuityPressureCoeffs[p    ](LoIndex) += coeff1 * LoCellLengthInv;
                continuityPressureCoeffs[east ](LoIndex) += coeff2 * LoCellLengthInv;
                continuityPressureCoeffs[eeast](LoIndex)  = coeff3 * LoCellLengthInv;

                // Cell on east side
                floatType HiCellLengthInv = mesh.cellLengthsInv[axis]( HiIndex[axis] );
                continuityPressureCoeffs[wwest](HiIndex)  = - coeff0 * HiCellLengthInv;
                continuityPressureCoeffs[west ](HiIndex)  = - coeff1 * HiCellLengthInv;
                continuityPressureCoeffs[p    ](HiIndex) += - coeff2 * HiCellLengthInv;
                continuityPressureCoeffs[east ](HiIndex)  = - coeff3 * HiCellLengthInv;

            }
        }
    }

}



// Fully implicit momentum interpolation coefficient for internal faces
// Autovectorisation friendly version
template< Axis::ENUMDATA axis > [[maybe_unused]]
void MWInterpolationInteriorImplicit_autoVec( ContinuityEquation &continuityEquation,
                                              const array3D &momentumDiagCoeffInv,
                                              const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    // Unpack
    EnumVector<TransportCoefficients, array3D> &continuityPressureCoeffs = continuityEquation.AP;
    const std::array< array1D, 4 > &mwiSparseCoeffs                      = continuityEquation.mwiSparseCoeffs[axis];
    const std::array< array1D, 2 > &mwiCompactCoeffs                     = continuityEquation.mwiCompactCoeffs[axis];

    // Cell indexing
    TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                    eeast = LUT::HiHiCoeff[axis],
                                    west  = LUT::LoCoeff[axis],
                                    wwest = LUT::LoLoCoeff[axis];

    // Set the first most plane to zero only. We do this so that the high coefficients can be set in place, and less coefficients have to 
    // be zeroed upon re-linearisation.
    continuityPressureCoeffs[east ].chip(0, axis) = continuityPressureCoeffs[west ].chip(0, axis).constant( 0.0f );
    continuityPressureCoeffs[west ].chip(0, axis) = continuityPressureCoeffs[west ].chip(0, axis).constant( 0.0f );

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

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
            std::vector< std::array< floatType, 4 > > mwiLineCoeffs( nFaces[X] );

            CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D HiIndex = { i, j, k },
                             LoIndex = { i, j, k };
                LoIndex[axis] -= 1;

                floatType d =  0.5f * ( momentumDiagCoeffInv( HiIndex )  +  momentumDiagCoeffInv( LoIndex ) );   

                // Coefficients for westmost to eastmost cell
                intType idx = HiIndex[axis];
                mwiLineCoeffs[i][0] = d * mwiSparseCoeffs[0](idx),
                mwiLineCoeffs[i][1] = d * ( mwiSparseCoeffs[1](idx) + mwiCompactCoeffs[0](idx) ),
                mwiLineCoeffs[i][2] = d * ( mwiSparseCoeffs[2](idx) + mwiCompactCoeffs[1](idx) ),
                mwiLineCoeffs[i][3] = d * mwiSparseCoeffs[3](idx);
            }

            CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                if constexpr ( axis == X ) {
                    HiCellLengthInv = mesh.cellLengthsInv[axis]( i );
                }

                arrayIndex3D HiIndex = { i, j, k };

                // Cell on east side
                continuityPressureCoeffs[wwest](HiIndex)  = - mwiLineCoeffs[i][0] * HiCellLengthInv;
                continuityPressureCoeffs[west ](HiIndex)  = - mwiLineCoeffs[i][1] * HiCellLengthInv;
                continuityPressureCoeffs[p    ](HiIndex) += - mwiLineCoeffs[i][2] * HiCellLengthInv;
                continuityPressureCoeffs[east ](HiIndex)  = - mwiLineCoeffs[i][3] * HiCellLengthInv;
            }

            CFD_PRAGMA_VECTORIZE
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                if constexpr ( axis == X ) {
                    LoCellLengthInv = mesh.cellLengthsInv[axis]( i-1 );
                }

                arrayIndex3D LoIndex = { i, j, k };
                LoIndex[axis] -= 1;

                // Cell on west side 
                continuityPressureCoeffs[west ](LoIndex) += mwiLineCoeffs[i][0] * LoCellLengthInv;
                continuityPressureCoeffs[p    ](LoIndex) += mwiLineCoeffs[i][1] * LoCellLengthInv;
                continuityPressureCoeffs[east ](LoIndex) += mwiLineCoeffs[i][2] * LoCellLengthInv;
                continuityPressureCoeffs[eeast](LoIndex)  = mwiLineCoeffs[i][3] * LoCellLengthInv;
            }

        }
    }

}



// Semi explicit momentum interpolation coefficient for internal faces
void MWInterpolationInteriorSemiExplicit( ContinuityEquation &continuityEquation, 
                                          const array3D &P,
                                          const array3D &momentumDiagCoeffInv,
                                          const Mesh &mesh,
                                          const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    // Unpack
    EnumVector<TransportCoefficients, array3D> &continuityPressureCoeffs = continuityEquation.AP;
    array3D &continuitySourceTerm                                        = continuityEquation.B;
    const std::array< array1D, 4 > &mwiSparseCoeffs                      = continuityEquation.mwiSparseCoeffs[axis];
    const std::array< array1D, 2 > &mwiCompactCoeffs                     = continuityEquation.mwiCompactCoeffs[axis];

    // For getting the index of a neighbouring cell
    auto NeighbourIndex = [] ( arrayIndex3D index, intType shift, Axis::ENUMDATA shiftAxis ) { 
        index[shiftAxis] += shift; 
        return index; 
    };

    // Cell indexing
    TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                    west  = LUT::LoCoeff[axis];

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D HiIndex = { i, j, k },
                             LoIndex = { i, j, k };
                LoIndex[axis] -= 1;
                
                floatType d =  0.5f * ( momentumDiagCoeffInv( HiIndex )  +  momentumDiagCoeffInv( LoIndex ) );   

                floatType LoCellLengthInv = mesh.cellLengthsInv[axis]( LoIndex[axis] ),
                          HiCellLengthInv = mesh.cellLengthsInv[axis]( HiIndex[axis] );


                // Implicit compact difference --------------------------------------------------------------------------
                                

                // Coefficients for westmost to eastmost cell
                intType idx = HiIndex[axis];
                floatType coeffCompact0 = d * mwiCompactCoeffs[0](idx),
                          coeffCompact1 = d * mwiCompactCoeffs[1](idx);

                // Cell on west side 
                continuityPressureCoeffs[p    ](LoIndex) +=   coeffCompact0 * LoCellLengthInv;
                continuityPressureCoeffs[east ](LoIndex)  =   coeffCompact1 * LoCellLengthInv;

                // Cell on east side
                continuityPressureCoeffs[west ](HiIndex)  = - coeffCompact0 * HiCellLengthInv;
                continuityPressureCoeffs[p    ](HiIndex) += - coeffCompact1 * HiCellLengthInv;


                // Explicit sparse difference ---------------------------------------------------------------------------

                arrayIndex3D LoWest  = NeighbourIndex( LoIndex, -1, axis ),
                             LoEast  = NeighbourIndex( LoIndex,  1, axis ),
                             LoEEast = NeighbourIndex( LoIndex,  2, axis ),

                             HiWWest = NeighbourIndex( HiIndex, -2, axis ),
                             HiWest  = NeighbourIndex( HiIndex, -1, axis ),
                             HiEast  = NeighbourIndex( HiIndex,  1, axis );

                floatType coeffSparse0 = d * mwiSparseCoeffs[0](idx),
                          coeffSparse1 = d * mwiSparseCoeffs[1](idx),
                          coeffSparse2 = d * mwiSparseCoeffs[2](idx),
                          coeffSparse3 = d * mwiSparseCoeffs[3](idx);

                // Cell on west side 
                continuitySourceTerm(LoIndex) -= ( coeffSparse0 * P( G(LoWest)  )
                                                 + coeffSparse1 * P( G(LoIndex) )
                                                 + coeffSparse2 * P( G(LoEast)  )
                                                 + coeffSparse3 * P( G(LoEEast) )
                                                 ) * LoCellLengthInv;

                // Cell on east side
                continuitySourceTerm(HiIndex) += ( coeffSparse0 * P( G(HiWWest) )
                                                 + coeffSparse1 * P( G(HiWest)  )
                                                 + coeffSparse2 * P( G(HiIndex) )
                                                 + coeffSparse3 * P( G(HiEast)  )
                                                 ) * HiCellLengthInv;


                // ------------------------------------------------------------------------------------------------------

            }
        }
    }

}



// Boundary constants that come from MWI for negative side boundary
void MWInterpolationNegativeBoundary( EnumVector<BoundaryPatches, array2D> &continuityBoundaryPressure,
                                      const EnumVector<BoundaryPatches, array2D> &momentumBoundaryPressure,
                                      const array3D &momentumDiagCoeffInv, 
                                      const Mesh &mesh,
                                      const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;

    BoundaryPatches::ENUMDATA negativePatch = LUT::NegativePatch[ axis ];
    Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[axis],
                   axis2 = LUT::HiOrthogonalAxis[axis];

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);
    startIndex[axis] = 1;
    nFaces[axis] = startIndex[axis] + 1;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D idx = { i, j, k };
                floatType d = MWIWeightingCoeff( idx, momentumDiagCoeffInv, axis );
                continuityBoundaryPressure[ negativePatch ]( idx[axis1], idx[axis2] ) = d * (1 - mesh.interpFactors[axis]( idx[axis] )) * momentumBoundaryPressure[ negativePatch ]( idx[axis1], idx[axis2] )
                                                                                          * mesh.cellLengthsInv[axis]( idx[axis] );
            }
        }
    }

}



// Boundary constants that come from MWI for positive side boundary
void MWInterpolationPositiveBoundary( EnumVector<BoundaryPatches, array2D> &continuityBoundaryPressure,
                                      const EnumVector<BoundaryPatches, array2D> &momentumBoundaryPressure,
                                      const array3D &momentumDiagCoeffInv, 
                                      const Mesh &mesh,
                                      const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;

    BoundaryPatches::ENUMDATA positivePatch = LUT::PositivePatch[ axis ];
    Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[axis],
                   axis2 = LUT::HiOrthogonalAxis[axis];

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);
    startIndex[axis] = mesh.nCells[axis] - 1;
    nFaces[axis] = startIndex[axis] + 1;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D idx = { i, j, k };
                floatType d = MWIWeightingCoeff( idx, momentumDiagCoeffInv, axis );
                continuityBoundaryPressure[ positivePatch ]( idx[axis1], idx[axis2] ) = d * mesh.interpFactors[axis]( idx[axis] ) * momentumBoundaryPressure[ positivePatch ]( idx[axis1], idx[axis2] )
                                                                                          * mesh.cellLengthsInv[axis]( idx[axis] );
            }
        }
    }

}




// Set momentum interpolation coefficients
void SetMomentumInterpolationCoefficients( FVCoefficients &fvCoeffs,
                                           const Mesh &mesh,
                                           const FieldData< BoundaryConditionData > &bcData,
                          [[maybe_unused]] const array3D &P )
{
    // Assumes that 'p' coefficient is set to zero
    // Correction is zero at the boundary faces

    #if defined( CFD_USE_AUTOVEC_FUNCTIONS )
        using enum Axis::ENUMDATA;
        switch ( fvCoeffs.Cont.momentumInterpolation ) {

            case MomentumInterpolation::Implicit:
                MWInterpolationInteriorImplicit_autoVec<X>(fvCoeffs.Cont, fvCoeffs.Mom[X].diagCoeffInv, mesh);
                MWInterpolationInteriorImplicit_autoVec<Y>(fvCoeffs.Cont, fvCoeffs.Mom[Y].diagCoeffInv, mesh);
                MWInterpolationInteriorImplicit_autoVec<Z>(fvCoeffs.Cont, fvCoeffs.Mom[Z].diagCoeffInv, mesh);
                break;

            case MomentumInterpolation::SemiExplicit:
                    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                        MWInterpolationInteriorSemiExplicit(fvCoeffs.Cont, P, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);
                    } );
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
        
        if ( bcData.P[ LUT::PositivePatch[axis] ].type == BoundaryConditions::fixed ) 
            MWInterpolationPositiveBoundary(fvCoeffs.Cont.BPBoundary, fvCoeffs.Mom[axis].BPBoundary, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);

        if ( bcData.P[ LUT::NegativePatch[axis] ].type == BoundaryConditions::fixed ) 
            MWInterpolationNegativeBoundary(fvCoeffs.Cont.BPBoundary, fvCoeffs.Mom[axis].BPBoundary, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);

    } );

}

/*---------------------------------------------------------------------------------------------------------------*\
                                        Boundary Constants to Source Term
\*---------------------------------------------------------------------------------------------------------------*/

void AddMomentumBoundaryConstants( MomentumEquation &momCoeffs )
{
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    intType iEnd;

    // Each axis
    for (intType axis = 0; axis != Axis::count; axis++) {

        positivePatch = LUT::PositivePatch[ static_cast<size_t>( axis ) ];
        negativePatch = LUT::NegativePatch[ static_cast<size_t>( axis ) ];
        iEnd  = momCoeffs.B.dimension( static_cast<size_t>( axis ) ) - 1;

        // Negative side boundary
        if ( momCoeffs.BUBoundary[negativePatch].size() != 0 )
            momCoeffs.B.chip( 0   , axis ) -= momCoeffs.BUBoundary[negativePatch];

        if ( momCoeffs.BPBoundary[negativePatch].size() != 0 )
            momCoeffs.B.chip( 0   , axis ) -= momCoeffs.BPBoundary[negativePatch];


        // Positive side boundary
        if ( momCoeffs.BUBoundary[positivePatch].size() != 0 )
            momCoeffs.B.chip( iEnd, axis ) -= momCoeffs.BUBoundary[positivePatch];

        if ( momCoeffs.BPBoundary[positivePatch].size() != 0 )
            momCoeffs.B.chip( iEnd, axis ) -= momCoeffs.BPBoundary[positivePatch];

    }
}


void AddContinuityBoundaryConstants( ContinuityEquation &contCoeffs )
{
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    intType iEnd;

    // Each axis
    for (intType axis = 0; axis != Axis::count; axis++) {

        positivePatch = LUT::PositivePatch[ static_cast<size_t>( axis ) ];
        negativePatch = LUT::NegativePatch[ static_cast<size_t>( axis ) ];
        iEnd = contCoeffs.B.dimension( static_cast<size_t>( axis ) ) - 1;

        // Negative side boundary
        if ( contCoeffs.BUBoundary[negativePatch].size() != 0 )
            contCoeffs.B.chip( 0   , axis ) -= contCoeffs.BUBoundary[negativePatch];

        if ( contCoeffs.BPBoundary[negativePatch].size() != 0 )
            contCoeffs.B.chip( 0   , axis ) -= contCoeffs.BPBoundary[negativePatch];


        // Positive side boundary
        if ( contCoeffs.BUBoundary[positivePatch].size() != 0 )
            contCoeffs.B.chip( iEnd, axis ) -= contCoeffs.BUBoundary[positivePatch];

        if ( contCoeffs.BPBoundary[positivePatch].size() != 0 )
            contCoeffs.B.chip( iEnd, axis ) -= contCoeffs.BPBoundary[positivePatch];

    }
}



/*---------------------------------------------------------------------------------------------------------------*\
                                                General functions
\*---------------------------------------------------------------------------------------------------------------*/

// Allocate the 2D arrays which store constant values of boundary conditions if they are needed
void AllocateBoundaryConstants( FVCoefficients &fvCoeffs, 
                                const FieldData< BoundaryConditionData > &bcData )
{
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {

        // Dimensions of the patch
        Axis::ENUMDATA patchAxis = LUT::BoundaryPatchAxis[ bp ];
        intType patchDimLo = fvCoeffs.nCells( LUT::LoOrthogonalAxis[ patchAxis ] ),
                patchDimHi = fvCoeffs.nCells( LUT::HiOrthogonalAxis[ patchAxis ] );


        // Velocity boundary conditions
        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

            if ( bcData.U[axis][bp].type == BoundaryConditions::fixed ) {

                if ( patchAxis == axis ) {
                        if ( fvCoeffs.Mom[ LUT::LoOrthogonalAxis[patchAxis] ].linearisation == Linearisation::Newton ) 
                            fvCoeffs.Mom[ LUT::LoOrthogonalAxis[patchAxis] ].BUBoundary[bp] = array2D(patchDimLo, patchDimHi).setZero();

                        if ( fvCoeffs.Mom[ LUT::HiOrthogonalAxis[patchAxis] ].linearisation == Linearisation::Newton ) 
                            fvCoeffs.Mom[ LUT::HiOrthogonalAxis[patchAxis] ].BUBoundary[bp] = array2D(patchDimLo, patchDimHi).setZero();
                }


                fvCoeffs.Mom[axis].BUBoundary[bp] = array2D(patchDimLo, patchDimHi).setZero();
                fvCoeffs.Cont.BUBoundary[bp] = array2D(patchDimLo, patchDimHi).setZero();

            }

        } );


        // Pressure boundary conditions
        if ( bcData.P[bp].type == BoundaryConditions::fixed ) {

            EnumFor<Axis>( [&] (Axis::ENUMDATA momAxis) {
                fvCoeffs.Mom[momAxis].BPBoundary[bp] = array2D(patchDimLo, patchDimHi).setZero();
            } );
            fvCoeffs.Cont.BPBoundary[bp] = array2D(patchDimLo, patchDimHi).setZero();

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

        fvCoeffs.Mom[axis].AU[axis][p].setZero();

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
        } );

        fvCoeffs.Mom[axis].B.setZero();

    } );

    
    // Zero continuity equation
    fvCoeffs.Cont.AP[p].setZero();

    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        fvCoeffs.Cont.BPBoundary[bp].setZero();
    } );

    fvCoeffs.Cont.B.setZero();
}


}   // end anonymous namespace





/*---------------------------------------------------------------------------------------------------------------*\
                                            Set and Update Functions
\*---------------------------------------------------------------------------------------------------------------*/


// Allocate and initialise finite volume coefficients for momentum and continuity equations
FVCoefficients InitialiseFVCoefficients( const Mesh &mesh,
                                         const FieldData<array3D> &fields,
                                         const EnumVector< Axis, EnumVector< Axis, array3D> > &faceAdvectedVelocities,
                                         const EnumVector<Axis, array3D> &faceFluxes, 
                                         const FieldData< BoundaryConditionData > &bcData,
                                         const InputData &inputData)
{
    using TC = TransportCoefficients::ENUMDATA;

    // Default construct the coefficients class
    FVCoefficients fvCoeffs(mesh.nCells, inputData.schemes.linearisation, inputData.schemes.momentumInterpolation);

    // Allocate boundary constant arrays for fixed boundary conditions
    AllocateBoundaryConstants( fvCoeffs, bcData );

    // Momentum equations
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        // Diffusion coefficients
        SetDiffusionCoeffients(fvCoeffs.Mom[axis], bcData, inputData.nu, mesh);

        // Advection terms
        SetInteriorAdvectionPicardCoefficients(fvCoeffs.Mom[axis], faceFluxes, mesh);
        SetBoundaryAdvectionPicardCoefficients(fvCoeffs.Mom[axis], faceFluxes, bcData.U[axis], mesh);

        // Add diffusion to velocity terms
        AddDiffusion(fvCoeffs.Mom[axis], bcData.U[axis], mesh);

        // Inverse of AP coefficient (Picard)
        fvCoeffs.Mom[axis].diagCoeffInv = fvCoeffs.Mom[axis].AU[axis][TC::p].inverse();

        // Momentum pressure terms
        SetFaceInterpolatedCoefficients(fvCoeffs.Mom[axis].AP, fvCoeffs.Mom[axis].BPBoundary, mesh, bcData.P, axis);
        DivideMomentumPressureByDensity(fvCoeffs.Mom[axis], inputData.rho);

        // Relaxation factor
        fvCoeffs.Mom[axis].relaxation = inputData.schemes.implicitRelaxation.U[axis];

    } );


    // Continuity equation
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        // Velocity terms
        SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AU[axis], fvCoeffs.Cont.BUBoundary, mesh, bcData.U[axis], axis);

        // Momentum weighted interpolation constants in the coefficients
        SetMomentumInterpolationSparseConstants(fvCoeffs.Cont.mwiSparseCoeffs[axis], fvCoeffs.Mom[axis].AP, mesh, axis);
        SetMomentumInterpolationCompactConstants(fvCoeffs.Cont.mwiCompactCoeffs[axis], inputData.rho, mesh, axis);

    } );


    // Momentum Weighted interpolation
    SetMomentumInterpolationCoefficients(fvCoeffs, mesh, bcData, fields.P);
    
    // Add Newton Linearisation terms if selected
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        if ( fvCoeffs.Mom[axis].linearisation == Linearisation::Newton ) {
            AddAdvectionNewtonCoefficients(fvCoeffs.Mom[axis], faceAdvectedVelocities, faceFluxes, bcData.U, mesh);
            fvCoeffs.Mom[axis].diagCoeffInv = fvCoeffs.Mom[axis].AU[axis][TC::p].inverse();
        }
            
        // Add boundary constants to source terms
        AddMomentumBoundaryConstants(fvCoeffs.Mom[axis]);
    } );

    AddContinuityBoundaryConstants(fvCoeffs.Cont);

    // Relaxation factor
    fvCoeffs.Cont.relaxation = inputData.schemes.implicitRelaxation.P;

    return fvCoeffs;
}





// Update linearisation in momenum and continuity equations
void UpdateFVCoefficients( FVCoefficients &fvCoeffs, 
                           const Mesh &mesh,
                           const FieldData<array3D> &fields,
                           const EnumVector< Axis, EnumVector< Axis, array3D> > &faceAdvectedVelocities,
                           const EnumVector<Axis, array3D> &faceFluxes,
                           const FieldData< BoundaryConditionData > &bcData )
{
    using TC = TransportCoefficients::ENUMDATA;

    TIC("Coefficient Update")

    TIC("Zeroing")
    ZeroNonlinearCoeffs( fvCoeffs );
    TOC()

    // The Picard coefficients for all momentum equations are the same, so just use the ones from the U momentum equation after 
    // its been set
    TIC("Picard advection")
    SetInteriorAdvectionPicardCoefficients(fvCoeffs.Mom[Axis::X], faceFluxes, mesh);
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        if ( axis != Axis::X ) {
            EnumFor<TransportCoefficients>( [&] (TransportCoefficients::ENUMDATA tc) {
                fvCoeffs.Mom[axis].AU[axis][tc] = fvCoeffs.Mom[Axis::X].AU[Axis::X][tc];
            } );
        }
    } );

    // Boundaries need to be done after since they can affect the internal coefficients
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        SetBoundaryAdvectionPicardCoefficients(fvCoeffs.Mom[axis], faceFluxes, bcData.U[axis], mesh);
    } );
    TOC()


    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        // Add diffusion to velocity terms
        TIC("Add diffusion")
        AddDiffusion(fvCoeffs.Mom[axis], bcData.U[axis], mesh);
        TOC()

        // Inverse of AP coefficient (Picard)
        TIC("Inverse AP (Picard)")
        fvCoeffs.Mom[axis].diagCoeffInv = fvCoeffs.Mom[axis].AU[axis][TC::p].inverse();
        TOC();

    } );

    // Set the momentum interpolation coefficients
    TIC("MWI")
    SetMomentumInterpolationCoefficients(fvCoeffs, mesh, bcData, fields.P);
    TOC()

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        // Add Newton Linearisation terms if selected
        TIC("Newton advection")
        if ( fvCoeffs.Mom[axis].linearisation == Linearisation::Newton ) {
            AddAdvectionNewtonCoefficients(fvCoeffs.Mom[axis], faceAdvectedVelocities, faceFluxes, bcData.U, mesh);
            fvCoeffs.Mom[axis].diagCoeffInv = fvCoeffs.Mom[axis].AU[axis][TC::p].inverse();
        }
        TOC()

         // Add boundary constants to source terms
        TIC("Momentum boundary constants")
        AddMomentumBoundaryConstants(fvCoeffs.Mom[axis]);
        TOC()

    } );

    TIC("Continuity boundary constants")
    AddContinuityBoundaryConstants(fvCoeffs.Cont);
    TOC()
    

    TOC()
}


}   // end namespace CFD
