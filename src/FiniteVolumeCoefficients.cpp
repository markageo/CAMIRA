#include "FiniteVolume.h"

#include "Utils.h"

#include <algorithm>
#include <iostream>

namespace CFD
{
 
namespace
{

/*---------------------------------------------------------------------------------------------------------------*\
                                                    Diffusion
\*---------------------------------------------------------------------------------------------------------------*/

// Check if continuity equation implies a zero gradient boundary condition. This occurs if both orthogonal fields have a uniform BC
BoundaryConditions::ENUMDATA GetDiffusionBC( const EnumVector< Axis, EnumVector< BoundaryPatches, InputData::BoundaryConditionData > > &MomBoundaryConditions, 
                                             const BoundaryPatches::ENUMDATA boundaryPatch, 
                                             const Axis::ENUMDATA velocityComponent )
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;

    const Axis::ENUMDATA axis = BoundaryPatchAxis[boundaryPatch];

    // Set the field we need to check based on the axis
    Axis::ENUMDATA axis1 = ( axis == Axis::X ) ? Axis::Y : Axis::X;
    Axis::ENUMDATA axis2 = ( axis == Axis::Z ) ? Axis::Y : Axis::Z;


    // Only check the field that in the direction of the current axis
    if (velocityComponent == axis) {

        if (MomBoundaryConditions[axis1][boundaryPatch].type == BC::uniform && 
            MomBoundaryConditions[axis2][boundaryPatch].type == BC::uniform) {
            return BC::zeroGradient;
        }
    }

    return MomBoundaryConditions[velocityComponent][boundaryPatch].type;
}


// Apply boundary conditions for diffusion terms on axis positive boundary
void DiffusionPositiveBoundary( EnumVector< Axis,  ArrayAllocator<TransportCoefficients, array1D> > &diff, 
                                EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = LoCoeff[axis];
    const intType iCellBound = mesh.nCells(axis) - 1;

    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient: 
            /* NULL */
            break;

        case BC::uniform:
            diff[axis][p   ](iCellBound)     +=   2*mesh.cellLengthsInv[axis](iCellBound);
            boundaryConstants[boundaryPatch] += - 2*mesh.cellLengthsInv[axis](iCellBound) * boundaryConditionStructs[boundaryPatch].value;
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
void DiffusionNegativeBoundary( EnumVector< Axis, ArrayAllocator<TransportCoefficients, array1D> > &diff, 
                                EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = NegativePatch[axis];
    const TransportCoefficients::ENUMDATA east = HiCoeff[axis];
    const intType iCellBound = 0;

    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient: 
            /* NULL */
            break;

        case BC::uniform:
            diff[axis][p   ](iCellBound)     +=   2*mesh.cellLengthsInv[axis](iCellBound);
            boundaryConstants[boundaryPatch] += - 2*mesh.cellLengthsInv[axis](iCellBound) * boundaryConditionStructs[boundaryPatch].value;
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
void SetDiffusionCoeffients(EnumVector< Axis, ArrayAllocator<TransportCoefficients, array1D> > &diff, 
                            EnumVector< BoundaryPatches, floatType > &boundaryConstants, 
                            const Mesh &mesh, 
                            const InputData &inputData, 
                            const Axis::ENUMDATA velocityComponent)
{

    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    
    const auto &boundaryConditions = inputData.boundaryConditions;

    TransportCoefficients::ENUMDATA east, west;     // These are just names, they can be north, south etc.
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    BoundaryConditions::ENUMDATA positivePatchBC, negativePatchBC;  // Store these since the continuity equation can override a BC to be zeroGradient


    // Diffusion in each axis is calculated in the same way
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        positivePatch = PositivePatch[axis];
        negativePatch = NegativePatch[axis];
        east = HiCoeff[axis];
        west = LoCoeff[axis];     

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
        positivePatchBC = GetDiffusionBC(boundaryConditions.U, positivePatch, velocityComponent);
        negativePatchBC = GetDiffusionBC(boundaryConditions.U, negativePatch, velocityComponent); 


        // Boundary conditions only need to be set if it is not zero gradient
        if (positivePatchBC != BC::zeroGradient) {
            DiffusionPositiveBoundary(diff, boundaryConstants, mesh, boundaryConditions.U[velocityComponent], axis);
        }

        if (negativePatchBC != BC::zeroGradient) {
            DiffusionNegativeBoundary(diff, boundaryConstants, mesh, boundaryConditions.U[velocityComponent], axis);
        }


        // Multiply by inverse cell length
        for (intType i = 0; i != mesh.nCells(axis); i++) {
            diff[axis][p   ](i) *= mesh.cellLengthsInv[axis](i);
            diff[axis][east](i) *= mesh.cellLengthsInv[axis](i);
            diff[axis][west](i) *= mesh.cellLengthsInv[axis](i);
        }
        boundaryConstants[ PositivePatch[axis] ] *= mesh.cellLengthsInv[axis]( mesh.nCells(axis)-1 );
        boundaryConstants[ NegativePatch[axis] ] *= mesh.cellLengthsInv[axis]( 0 );

        // Multiply by viscosity
        for (intType i = 0; i != mesh.nCells(axis); i++) {
            diff[axis][p   ](i) *= inputData.nu;
            diff[axis][east](i) *= inputData.nu;
            diff[axis][west](i) *= inputData.nu;
        }
        boundaryConstants[ PositivePatch[axis] ] *= inputData.nu;
        boundaryConstants[ NegativePatch[axis] ] *= inputData.nu;

    } );
}





/*---------------------------------------------------------------------------------------------------------------*\
                                         Momentum Advection Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


// Upwind coefficients
void Upwind( ArrayAllocator<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
             const EnumVector<Axis, array3D> &faceFluxes, 
             const Mesh &mesh,
             const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    // Starting index and number of faces to iterate over
    iVector3 startIndex, nFaces;
    EnumFor<Axis>( [&] ( Axis::ENUMDATA a) {
        startIndex[a] = 0;
        nFaces[a] = faceFluxes[ axis ].dimension(a);
    } );
    startIndex[axis] += 1;
    nFaces[axis] -= 1;

    TransportCoefficients::ENUMDATA east = HiCoeff[axis], 
                                    west = LoCoeff[axis];

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {
                
                arrayIndex3D HiIndex = { i, j, k },
                             LoIndex = { i, j, k };
                LoIndex[axis] -= 1;

                floatType uf = faceFluxes[ axis ](i, j, k);
                if ( uf >= 0 ) {
                    
                    coeffs[p   ](LoIndex) +=   uf * mesh.cellLengthsInv[axis]( LoIndex[axis] );
                    coeffs[west](HiIndex)  = - uf * mesh.cellLengthsInv[axis]( HiIndex[axis] );

                } else {

                    coeffs[east](LoIndex)  =   uf * mesh.cellLengthsInv[axis]( LoIndex[axis] );
                    coeffs[p   ](HiIndex) += - uf * mesh.cellLengthsInv[axis]( HiIndex[axis] );

                }

            }
        }
    }
}


void AdvectionPositiveBoundary( ArrayAllocator<TransportCoefficients, array3D> &coeffs, 
                                EnumVector<BoundaryPatches, array2D> &boundaryConstants,
                                const EnumVector<Axis, array3D> &faceFluxes, 
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = LoCoeff[axis];
    const intType iCellBound = mesh.nCells(axis) - 1;   // Index of cell at the boundary
    const intType iFaceBound = iCellBound + 1;          // Index of face at the boundary

    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p].chip(iCellBound, axis) += faceFluxes[axis].chip(iFaceBound, axis) 
                                              * faceFluxes[axis].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::uniform:
            boundaryConstants[boundaryPatch]  = faceFluxes[axis].chip(iFaceBound, axis)
                                              * faceFluxes[axis].chip(iFaceBound, axis).constant( boundaryConditionStructs[boundaryPatch].value * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::extrapolated:
            coeffs[p   ].chip(iCellBound, axis) += faceFluxes[axis].chip(iFaceBound, axis) 
                                                 * faceFluxes[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].p * mesh.cellLengthsInv[axis](iCellBound) );
            coeffs[west].chip(iCellBound, axis) += faceFluxes[axis].chip(iFaceBound, axis) 
                                                 * faceFluxes[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].a * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        default:
            break;
    }

}


void AdvectionNegativeBoundary( ArrayAllocator<TransportCoefficients, array3D> &coeffs, 
                                EnumVector<BoundaryPatches, array2D> &boundaryConstants,
                                const EnumVector<Axis, CFD::array3D> &faceFluxes, 
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = NegativePatch[axis];
    const TransportCoefficients::ENUMDATA east = HiCoeff[axis];
    const intType iCellBound = 0;   // Index of cell at the boundary 
    const intType iFaceBound = 0;   // Index of face at the boundary

    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p].chip(iCellBound, axis) += - faceFluxes[axis].chip(iFaceBound, axis) 
                                              *   faceFluxes[axis].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::uniform:
            boundaryConstants[boundaryPatch]  = - faceFluxes[axis].chip(iFaceBound, axis)
                                              *   faceFluxes[axis].chip(iFaceBound, axis).constant( boundaryConditionStructs[boundaryPatch].value * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::extrapolated:
            coeffs[p   ].chip(iCellBound, axis) += - faceFluxes[axis].chip(iFaceBound, axis) 
                                                 *   faceFluxes[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].p * mesh.cellLengthsInv[axis](iCellBound) );
            coeffs[east].chip(iCellBound, axis) += - faceFluxes[axis].chip(iFaceBound, axis) 
                                                 *   faceFluxes[axis].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].a * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        default:
            break;
    }

}


void SetAdvectionCoefficients( ArrayAllocator<TransportCoefficients, array3D> &coeffs, 
                               EnumVector<BoundaryPatches, array2D> &boundaryConstants,
                               const EnumVector<Axis, array3D> &faceFluxes, 
                               const Mesh &mesh, 
                               const InputData &inputData, 
                               const Axis::ENUMDATA velocityComponent)
{
    using enum TransportCoefficients::ENUMDATA;

    const auto &boundaryConditions = inputData.boundaryConditions;
    
    // For now assumes that all coefficiencients are set to zero

    // Upwind internal faces
    EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {

        // Upwind internal faces
        Upwind(coeffs, faceFluxes, mesh, axis);

        // Boundary faces
        AdvectionPositiveBoundary(coeffs, boundaryConstants, faceFluxes, mesh, boundaryConditions.U[velocityComponent], axis);
        AdvectionNegativeBoundary(coeffs, boundaryConstants, faceFluxes, mesh, boundaryConditions.U[velocityComponent], axis);

    } );

}





/*---------------------------------------------------------------------------------------------------------------*\
                                           Add Diffusion Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

void AddDiffusion( ArrayAllocator< TransportCoefficients, array3D > &velCoeffs, 
                   EnumVector< BoundaryPatches, array2D > &boundaryVel,
                   const EnumVector< Axis, ArrayAllocator<TransportCoefficients, array1D> > &diffCoeffs, 
                   const EnumVector< BoundaryPatches, floatType > &boundaryDiff,
                   const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    // Velocity coefficients
    for (intType k = 0; k != mesh.nCells(Z); k++) {
        for (intType j = 0; j != mesh.nCells(Y); j++) {
            for (intType i = 0; i != mesh.nCells(X); i++) {

                velCoeffs[p](i, j, k) += diffCoeffs[X][p](i) + diffCoeffs[Y][p](j) + diffCoeffs[Z][p](k);

                velCoeffs[e](i, j, k) += diffCoeffs[X][e](i);
                velCoeffs[w](i, j, k) += diffCoeffs[X][w](i);
                
                velCoeffs[n](i, j, k) += diffCoeffs[Y][n](j);
                velCoeffs[s](i, j, k) += diffCoeffs[Y][s](j);

                velCoeffs[t](i, j, k) += diffCoeffs[Z][t](k);
                velCoeffs[b](i, j, k) += diffCoeffs[Z][b](k);
            }
        }
    }

    // Constant terms
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA patch) {

        boundaryVel[patch] += boundaryVel[patch].constant( boundaryDiff[patch] );
        
    } );

}





/*---------------------------------------------------------------------------------------------------------------*\
                                        Linear Interpolated Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


void InterpolationPositiveBoundary( ArrayAllocator< TransportCoefficients, array1D > &coeffs, 
                                    EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                    const Mesh &mesh,  
                                    const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs,
                                    const Axis::ENUMDATA axis)
{

    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = LoCoeff[axis];
    const intType iCellBound = mesh.nCells(axis) - 1;

    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p   ]( iCellBound ) += 1;
            coeffs[west]( iCellBound ) += 0;
            break;

        case BC::uniform:
            coeffs[p   ]( iCellBound ) += 0;
            coeffs[west]( iCellBound ) += 0;
            boundaryConstants[boundaryPatch] += boundaryConditionStructs[boundaryPatch].value;
            break;

        case BC::extrapolated:
            coeffs[p   ]( iCellBound ) += mesh.extrapFactors[boundaryPatch].p;
            coeffs[west]( iCellBound ) += mesh.extrapFactors[boundaryPatch].a;
            break;

        default:
            break;
    }

}


void InterpolationNegativeBoundary( ArrayAllocator< TransportCoefficients, array1D > &coeffs, 
                                    EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                    const Mesh &mesh,  
                                    const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs, 
                                    const Axis::ENUMDATA axis)
{

    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = NegativePatch[axis];
    const TransportCoefficients::ENUMDATA east = HiCoeff[axis];
    const intType iCellBound = 0;

    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p   ]( iCellBound ) += - 1;
            coeffs[east]( iCellBound ) +=   0;
            break;

        case BC::uniform:
            coeffs[p   ]( iCellBound ) += 0;
            coeffs[east]( iCellBound ) += 0;
            boundaryConstants[boundaryPatch] += - boundaryConditionStructs[boundaryPatch].value;
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
void SetFaceInterpolatedCoefficients( ArrayAllocator<CFD::TransportCoefficients, CFD::array1D> &coeffs, 
                                      EnumVector< BoundaryPatches, floatType > &boundaryConstants, 
                                      const Mesh &mesh, 
                                      const EnumVector< BoundaryPatches, InputData::BoundaryConditionData> &boundaryConditionStructs, 
                                      const Axis::ENUMDATA axis)
{ 
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    TransportCoefficients::ENUMDATA east = HiCoeff[axis],    // These are just names, they can be north, south etc.
                                    west = LoCoeff[axis];  

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
    InterpolationPositiveBoundary(coeffs, boundaryConstants, mesh, boundaryConditionStructs, axis);
    InterpolationNegativeBoundary(coeffs, boundaryConstants, mesh, boundaryConditionStructs, axis);
    
    // Multiply by inverse of cell length
    for (intType i = 0; i != mesh.nCells(axis); i++) {
        coeffs[p   ](i) *= mesh.cellLengthsInv[axis](i);
        coeffs[east](i) *= mesh.cellLengthsInv[axis](i);
        coeffs[west](i) *= mesh.cellLengthsInv[axis](i); 
    }
    boundaryConstants[ PositivePatch[axis] ] *= mesh.cellLengthsInv[axis]( mesh.nCells(axis)-1 );
    boundaryConstants[ NegativePatch[axis] ] *= mesh.cellLengthsInv[axis]( 0 );
}



// Divide pressure terms in momentum equaqtions by density
void DivideMomentumPressureByDensity( ArrayAllocator<CFD::TransportCoefficients, CFD::array1D> &coeffs, 
                                      EnumVector< BoundaryPatches, floatType > &boundaryConstants, 
                                      const floatType rho, 
                                      const Axis::ENUMDATA axis)
{ 
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    TransportCoefficients::ENUMDATA east = HiCoeff[axis],    // These are just names, they can be north, south etc.
                                    west = LoCoeff[axis];  

    
    coeffs[p   ] /= coeffs[p   ].constant( rho );
    coeffs[east] /= coeffs[west].constant( rho );
    coeffs[west] /= coeffs[west].constant( rho );

    boundaryConstants[ PositivePatch[axis] ] /= rho;
    boundaryConstants[ NegativePatch[axis] ] /= rho;
}



/*---------------------------------------------------------------------------------------------------------------*\
                        Momentum Weighted Interpolation (Rhie-Chow Interpolation) Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

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


// Momentum interpolation coefficient for internal faces
void MWInterpolationFace( ArrayAllocator<TransportCoefficients, array3D> &continuityPressureCoeffs,
                          const array3D &momentumDiagCoeffInv,
                          const ArrayAllocator<TransportCoefficients, array1D> &momentumPressureCoeffs,
                          const Mesh &mesh, 
                          const floatType rho,
                          const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    
    // Starting index and number of faces to iterate over
    iVector3 startIndex, nFaces;
    EnumFor<Axis>( [&] ( Axis::ENUMDATA a) {
        startIndex[a] = 0;
        nFaces[a] = mesh.nCells[a];
    } );
    startIndex[axis] += 1;


    // Cell indexing
    TransportCoefficients::ENUMDATA east  = HiCoeff[axis], 
                                    eeast = HiHiCoeff[axis],
                                    west  = LoCoeff[axis],
                                    wwest = LoLoCoeff[axis];

    floatType rhoInv = 1.0f / rho;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D HiIndex = { i, j, k },
                             LoIndex = { i, j, k };
                LoIndex[axis] -= 1;

                floatType d = MWIWeightingCoeff( {i, j, k}, momentumDiagCoeffInv, axis );

                // Coefficients for westmost to eastmost cell. These coefficients assume that the momentum equations 
                // have been divided through by the cell volume
                intType idx = HiIndex[axis];
                floatType coeff1 = d * (1 - mesh.interpFactors[axis](idx))   * momentumPressureCoeffs[west](idx-1);

                floatType coeff2 = d * ( (1 - mesh.interpFactors[axis](idx)) * momentumPressureCoeffs[p   ](idx-1)
                                       +  mesh.interpFactors[axis](idx)      * momentumPressureCoeffs[west](idx  )
                                       +  mesh.cellCenterDiffInv[axis](idx)  * rhoInv );

                floatType coeff3 = d * ( (1 - mesh.interpFactors[axis](idx)) * momentumPressureCoeffs[east](idx-1)
                                       + mesh.interpFactors[axis](idx)       * momentumPressureCoeffs[p   ](idx  )
                                       - mesh.cellCenterDiffInv[axis](idx)   * rhoInv );

                floatType coeff4 = d * mesh.interpFactors[axis](i) * momentumPressureCoeffs[east](i);


                // Cell on west side 
                floatType LoCellLengthInv = mesh.cellLengthsInv[axis]( LoIndex[axis] );
                continuityPressureCoeffs[west ](LoIndex) += coeff1 * LoCellLengthInv;
                continuityPressureCoeffs[p    ](LoIndex) += coeff2 * LoCellLengthInv;
                continuityPressureCoeffs[east ](LoIndex) += coeff3 * LoCellLengthInv;
                continuityPressureCoeffs[eeast](LoIndex) += coeff4 * LoCellLengthInv;

                // Cell on east side
                floatType HiCellLengthInv = mesh.cellLengthsInv[axis]( HiIndex[axis] );
                continuityPressureCoeffs[wwest](HiIndex) -= coeff1 * HiCellLengthInv;
                continuityPressureCoeffs[west ](HiIndex) -= coeff2 * HiCellLengthInv;
                continuityPressureCoeffs[p    ](HiIndex) -= coeff3 * HiCellLengthInv;
                continuityPressureCoeffs[east ](HiIndex) -= coeff4 * HiCellLengthInv;

            }
        }
    }

}



void MWInterpolationBoundary( EnumVector<BoundaryPatches, array2D> &continuityBoundaryPressure,
                              const EnumVector<BoundaryPatches, floatType> &momentumBoundaryPressure,
                              const array3D &momentumDiagCoeffInv, 
                              const Mesh &mesh,
                              const Axis::ENUMDATA axis )
{
    // Boundary condition contribution comes from cell face one off the boundary, since MWI correction is taken as zero
    // at the boundary.

    using enum Axis::ENUMDATA;

    BoundaryPatches::ENUMDATA positivePatch = PositivePatch[ axis ],
                              negativePatch = NegativePatch[ axis ];

    // Other orthogonal axis directions
    Axis::ENUMDATA axis1 = ( axis == X ) ? Y : X,
                   axis2 = ( axis == Z ) ? Y : Z;

    iVector3 startIndex, nCells;
    EnumFor<Axis>( [&] ( Axis::ENUMDATA a) {
        startIndex[a] = 0;
        nCells[a] = mesh.nCells[a];
    } );


    // Negative patch
    startIndex[axis] = 1;
    nCells[axis] = startIndex[axis] + 1;

    for (intType k = startIndex[Z]; k != nCells[Z]; k++) {
        for (intType j = startIndex[Y]; j != nCells[Y]; j++) {
            for (intType i = startIndex[X]; i != nCells[X]; i++) {

                arrayIndex3D idx = { i, j, k };
                floatType d = MWIWeightingCoeff( idx, momentumDiagCoeffInv, axis );
                continuityBoundaryPressure[ negativePatch ]( idx[axis1], idx[axis2] ) = d * (1 - mesh.interpFactors[axis]( idx[axis] )) * momentumBoundaryPressure[ negativePatch ] 
                                                                                            * mesh.cellLengthsInv[axis]( idx[axis] );
            }
        }
    }


    // Positive patch
    startIndex[axis] = mesh.nCells[axis] - 1;
    nCells[axis] = startIndex[axis] + 1;

    for (intType k = startIndex[Z]; k != nCells[Z]; k++) {
        for (intType j = startIndex[Y]; j != nCells[Y]; j++) {
            for (intType i = startIndex[X]; i != nCells[X]; i++) {

                arrayIndex3D idx = { i, j, k };
                floatType d = MWIWeightingCoeff( idx, momentumDiagCoeffInv, axis );
                continuityBoundaryPressure[ positivePatch ]( idx[axis1], idx[axis2] ) = d * mesh.interpFactors[axis]( idx[axis] ) * momentumBoundaryPressure[ positivePatch ]
                                                                                            * mesh.cellLengthsInv[axis]( idx[axis] );
            }
        }
    }


}



void SetMomentumInterpolationCoefficients( FVCoefficients &fvCoeffs, 
                                           const Mesh &mesh, 
                                           const InputData &inputData)
{
    const floatType rho = inputData.rho;

    // Assumes that coefficients are set to zero

    // Only internal cells need to be done, since correction is zero at boundary, and
    // sparse pressure gradient is taken from momentum equations, which already contain
    // the boundary condition.
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        // Internal faces
        MWInterpolationFace(fvCoeffs.Cont.AP, fvCoeffs.Mom[axis].diagCoeffInv, fvCoeffs.Mom[axis].AP, mesh, rho, axis);

        // Boundary constants
        MWInterpolationBoundary(fvCoeffs.Cont.boundaryP, fvCoeffs.Mom[axis].boundaryP, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);

    } );

}
 

/*---------------------------------------------------------------------------------------------------------------*\
                                                Implicit Relaxation
\*---------------------------------------------------------------------------------------------------------------*/

[[ maybe_unused ]]
void AddRelaxation( array3D &diagonalCoeffs,
                    array3D &sourceTerms,
                    const array3D &oldField,
                    const floatType &relaxationFactor )
{
    using enum TransportCoefficients::ENUMDATA;

    if ( relaxationFactor == 1.0f )
        return;

    // Add to the diagonal coefficient
    diagonalCoeffs *= diagonalCoeffs.constant( 1.0f / relaxationFactor );

    // Coefficients do not have ghost cells, so need to slice the field to work with it
    arrayIndex3D offsets = {nGhost, nGhost, nGhost},
                 extents = {oldField.dimension(0) - 2*nGhost,
                            oldField.dimension(1) - 2*nGhost,
                            oldField.dimension(2) - 2*nGhost }; 

    // Add to the source term, note the diagonal coefficient is already relaxed
    sourceTerms += diagonalCoeffs * sourceTerms.constant( 1 - relaxationFactor  )
                 * oldField.slice( offsets, extents );

}





/*---------------------------------------------------------------------------------------------------------------*\
                                        Boundary Constants to Source Term
\*---------------------------------------------------------------------------------------------------------------*/

void AddMomentumBoundaryConstants( FVCoefficients::MomentumEquation &momCoeffs )
{
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    intType iEnd;

    // Each axis
    for (intType axis = 0; axis != Axis::count; axis++) {

        positivePatch = PositivePatch[ static_cast<size_t>( axis ) ];
        negativePatch = NegativePatch[ static_cast<size_t>( axis ) ];
        iEnd  = momCoeffs.B.dimension( static_cast<size_t>( axis ) ) - 1;

        // Negative side boundary
        momCoeffs.B.chip( 0   , axis ) -= momCoeffs.boundaryVel[negativePatch]
                                        + momCoeffs.B.chip( 0   , axis ).constant( momCoeffs.boundaryP[negativePatch] );

        // Positive side boundary
        momCoeffs.B.chip( iEnd, axis ) -= momCoeffs.boundaryVel[positivePatch]
                                        + momCoeffs.B.chip( iEnd, axis ).constant( momCoeffs.boundaryP[positivePatch] );
    }
}



void AddContinuityBoundaryConstants( FVCoefficients::ContinuityEquation &contCoeffs )
{
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    intType iEnd;

    // Each axis
    for (intType axis = 0; axis != Axis::count; axis++) {

        positivePatch = PositivePatch[ static_cast<size_t>( axis ) ];
        negativePatch = NegativePatch[ static_cast<size_t>( axis ) ];
        iEnd = contCoeffs.B.dimension( static_cast<size_t>( axis ) ) - 1;

        // Negative side boundary
        contCoeffs.B.chip( 0   , axis ) -= contCoeffs.boundaryP[negativePatch]
                                         + contCoeffs.B.chip( 0   , axis ).constant( contCoeffs.boundaryVel[negativePatch] );

        // Positive side boundary
        contCoeffs.B.chip( iEnd, axis ) -= contCoeffs.boundaryP[positivePatch]
                                         + contCoeffs.B.chip( iEnd, axis ).constant( contCoeffs.boundaryVel[positivePatch] );
    }
}


}   // end anonymous namespace





/*---------------------------------------------------------------------------------------------------------------*\
                                            Set and Update Functions
\*---------------------------------------------------------------------------------------------------------------*/


// Allocate and initialise finite volume coefficients for momentum and continuity equations
FVCoefficients InitialiseFVCoefficients( const Mesh &mesh,
                                         const EnumVector<Axis, array3D> &faceFluxes, 
                                         const InputData &inputData)
{
    using enum Axis::ENUMDATA;

    // Default construct the coefficients class
    FVCoefficients fvCoeffs(mesh.nCells);

    // Diffusion coefficients
    SetDiffusionCoeffients(fvCoeffs.Mom[X].diff, fvCoeffs.Mom[X].boundaryDiff, mesh, inputData, X);
    SetDiffusionCoeffients(fvCoeffs.Mom[Y].diff, fvCoeffs.Mom[Y].boundaryDiff, mesh, inputData, Y);
    SetDiffusionCoeffients(fvCoeffs.Mom[Z].diff, fvCoeffs.Mom[Z].boundaryDiff, mesh, inputData, Z);

    // Momentum advection terms
    SetAdvectionCoefficients(fvCoeffs.Mom[X].AU[X], fvCoeffs.Mom[X].boundaryVel, faceFluxes, mesh, inputData, X);
    SetAdvectionCoefficients(fvCoeffs.Mom[Y].AU[Y], fvCoeffs.Mom[Y].boundaryVel, faceFluxes, mesh, inputData, Y);
    SetAdvectionCoefficients(fvCoeffs.Mom[Z].AU[Z], fvCoeffs.Mom[Z].boundaryVel, faceFluxes, mesh, inputData, Z);

    // Add diffusion to the velocity coefficients in momentum equations
    AddDiffusion(fvCoeffs.Mom[X].AU[X], fvCoeffs.Mom[X].boundaryVel, fvCoeffs.Mom[X].diff, fvCoeffs.Mom[X].boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Mom[Y].AU[Y], fvCoeffs.Mom[Y].boundaryVel, fvCoeffs.Mom[Y].diff, fvCoeffs.Mom[Y].boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Mom[Z].AU[Z], fvCoeffs.Mom[Z].boundaryVel, fvCoeffs.Mom[Z].diff, fvCoeffs.Mom[Z].boundaryDiff, mesh);

    // Inverse of AP coefficient
    using TC = TransportCoefficients::ENUMDATA;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        fvCoeffs.Mom[axis].diagCoeffInv = fvCoeffs.Mom[axis].AU[axis][TC::p].inverse();
    } );

    // Momentum pressure terms
    SetFaceInterpolatedCoefficients(fvCoeffs.Mom[X].AP, fvCoeffs.Mom[X].boundaryP, mesh, inputData.boundaryConditions.P, X);
    SetFaceInterpolatedCoefficients(fvCoeffs.Mom[Y].AP, fvCoeffs.Mom[Y].boundaryP, mesh, inputData.boundaryConditions.P, Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.Mom[Z].AP, fvCoeffs.Mom[Z].boundaryP, mesh, inputData.boundaryConditions.P, Z);

    DivideMomentumPressureByDensity(fvCoeffs.Mom[X].AP, fvCoeffs.Mom[X].boundaryP, inputData.rho, X);
    DivideMomentumPressureByDensity(fvCoeffs.Mom[Y].AP, fvCoeffs.Mom[Y].boundaryP, inputData.rho, Y);
    DivideMomentumPressureByDensity(fvCoeffs.Mom[Z].AP, fvCoeffs.Mom[Z].boundaryP, inputData.rho, Z);

    // Continuity velocity terms
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AU[X], fvCoeffs.Cont.boundaryVel, mesh, inputData.boundaryConditions.U[X], X);
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AU[Y], fvCoeffs.Cont.boundaryVel, mesh, inputData.boundaryConditions.U[Y], Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AU[Z], fvCoeffs.Cont.boundaryVel, mesh, inputData.boundaryConditions.U[Z], Z);

    // Continuity pressure terms (from momentum weighted interpolation)
    SetMomentumInterpolationCoefficients(fvCoeffs, mesh, inputData);
    
    // Add boundary constants to source terms
    AddMomentumBoundaryConstants(fvCoeffs.Mom[X]);
    AddMomentumBoundaryConstants(fvCoeffs.Mom[Y]);
    AddMomentumBoundaryConstants(fvCoeffs.Mom[Z]);
    AddContinuityBoundaryConstants(fvCoeffs.Cont);

    // Set implicit under relaxation
    EnumFor<Axis> ( [&] (Axis::ENUMDATA axis) {
        fvCoeffs.Mom[axis].relaxation = inputData.schemes.implicitRelaxation.U[axis];
    } );
    fvCoeffs.Cont.relaxation = inputData.schemes.implicitRelaxation.P;

    return fvCoeffs;
}


// Update linearisation in momenum and continuity equations
void UpdateFVCoefficients(FVCoefficients &fvCoeffs, 
                          const Mesh &mesh,
                          const EnumVector<Axis, array3D> &faceFluxes,
                          const InputData &inputData)
{
    using enum Axis::ENUMDATA;

    TIC("Updating FV Coeffs")

    TIC("Zeroing coeffs")
    // Zero momentum equations
    EnumFor<Axis> ( [&] (Axis::ENUMDATA axis) {

        EnumFor<TransportCoefficients>( [&] (TransportCoefficients::ENUMDATA tc) {
            if ( fvCoeffs.Mom[axis].AU[axis].get(tc) )
                fvCoeffs.Mom[axis].AU[axis][tc].setZero();
        } );

        EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
            fvCoeffs.Mom[axis].boundaryVel[bp].setZero();
        } );

        fvCoeffs.Mom[axis].B.setZero();

    } );


    // Zero continuity equation
    EnumFor<TransportCoefficients>( [&] (TransportCoefficients::ENUMDATA tc) {
        if ( fvCoeffs.Cont.AP.get(tc) )
            fvCoeffs.Cont.AP[tc].setZero();
    } );

    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        fvCoeffs.Cont.boundaryP[bp].setZero();
    } );
    fvCoeffs.Cont.B.setZero();
    TOC()


    TIC("Advection")
    // Momentum advection terms
    SetAdvectionCoefficients(fvCoeffs.Mom[X].AU[X], fvCoeffs.Mom[X].boundaryVel, faceFluxes, mesh, inputData, X);
    SetAdvectionCoefficients(fvCoeffs.Mom[Y].AU[Y], fvCoeffs.Mom[Y].boundaryVel, faceFluxes, mesh, inputData, Y);
    SetAdvectionCoefficients(fvCoeffs.Mom[Z].AU[Z], fvCoeffs.Mom[Z].boundaryVel, faceFluxes, mesh, inputData, Z);
    TOC()

    TIC("Add diffusion")
    // Add diffusion to the velocity coefficients in momentum equations
    AddDiffusion(fvCoeffs.Mom[X].AU[X], fvCoeffs.Mom[X].boundaryVel, fvCoeffs.Mom[X].diff, fvCoeffs.Mom[X].boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Mom[Y].AU[Y], fvCoeffs.Mom[Y].boundaryVel, fvCoeffs.Mom[Y].diff, fvCoeffs.Mom[Y].boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Mom[Z].AU[Z], fvCoeffs.Mom[Z].boundaryVel, fvCoeffs.Mom[Z].diff, fvCoeffs.Mom[Z].boundaryDiff, mesh);
    TOC()

    TIC("Inverse AP")
    // Inverse of AP coefficient
    using TC = TransportCoefficients::ENUMDATA;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        fvCoeffs.Mom[axis].diagCoeffInv = fvCoeffs.Mom[axis].AU[axis][TC::p].inverse();
    } );
    TOC()

    TIC("MWI")
    // Set the momentum interpolation coefficients
    SetMomentumInterpolationCoefficients(fvCoeffs, mesh, inputData);
    TOC()

    TIC("Boundary constants")
    // Add boundary constants to source terms
    AddMomentumBoundaryConstants(fvCoeffs.Mom[X]);
    AddMomentumBoundaryConstants(fvCoeffs.Mom[Y]);
    AddMomentumBoundaryConstants(fvCoeffs.Mom[Z]);
    AddContinuityBoundaryConstants(fvCoeffs.Cont);
    TOC()

    TOC()

}


}   // end namespace CFD
