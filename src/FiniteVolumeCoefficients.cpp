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
BoundaryConditions::ENUMDATA GetDiffusionBC( const InputData::BoundaryConditionData &boundaryConditions, 
                                             const BoundaryPatches::ENUMDATA boundaryPatch, 
                                             const Fields::ENUMDATA field )
{
    using BC = BoundaryConditions::ENUMDATA;
    using F = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;

    Fields::ENUMDATA fieldToCheck = F::U;
    Fields::ENUMDATA orthogonalField1 = F::V; 
    Fields::ENUMDATA orthogonalField2 = F::W;
    const Axis::ENUMDATA axis = BoundaryPatchAxis[boundaryPatch];

    // Set the field we need to check based on the axis
    if        (axis == X) {
        fieldToCheck = F::U;
        orthogonalField1 = F::V;
        orthogonalField2 = F::W;

    } else if (axis == Y) {
        fieldToCheck = F::V;
        orthogonalField1 = F::W;
        orthogonalField2 = F::U;

    } else if (axis == Z) {
        fieldToCheck = F::W;
        orthogonalField1 = F::U;
        orthogonalField2 = F::V;

    } else {
        /* NULL */
        
    }

    // Only check the field that in the direction of the current axis
    if (field == fieldToCheck) {

        if (boundaryConditions[orthogonalField1][boundaryPatch].type == BC::uniform && 
            boundaryConditions[orthogonalField2][boundaryPatch].type == BC::uniform) {
            return BC::zeroGradient;
        }
    }

    return boundaryConditions[field][boundaryPatch].type;
}


// Apply boundary conditions for diffusion terms on axis positive boundary
void DiffusionPositiveBoundary( EnumVector< Axis,  ArrayAllocator<TransportCoefficients, array1D> > &diff, 
                                EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionStruct > &boundaryConditionStructs,
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
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionStruct > &boundaryConditionStructs,
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
                            const Fields::ENUMDATA field)
{

    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    
    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;

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
        positivePatchBC = GetDiffusionBC(boundaryConditions, positivePatch, field);
        negativePatchBC = GetDiffusionBC(boundaryConditions, negativePatch, field); 


        // Boundary conditions only need to be set if it is not zero gradient
        if (positivePatchBC != BC::zeroGradient) {
            DiffusionPositiveBoundary(diff, boundaryConstants, mesh, boundaryConditions[field], axis);
        }

        if (negativePatchBC != BC::zeroGradient) {
            DiffusionNegativeBoundary(diff, boundaryConstants, mesh, boundaryConditions[field], axis);
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
             const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
             const Mesh &mesh,
             const Axis::ENUMDATA axis )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    Fields::ENUMDATA field = AxisVelocity[axis];

    // Starting index and number of faces to iterate over
    iVector3 startIndex, nFaces;
    EnumFor<Axis>( [&] ( Axis::ENUMDATA a) {
        startIndex[a] = 0;
        nFaces[a] = faceVelocities[ field ].dimension(a);
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

                floatType uf = faceVelocities[ field ](i, j, k);
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
                                const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionStruct > &boundaryConditionStructs,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    static const std::array<Fields::ENUMDATA, 3> faceVelocityFields = {F::U, F::V, F::W}; // Used to get corresponding velocity field from axis
    
    const BoundaryPatches::ENUMDATA boundaryPatch = PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = LoCoeff[axis];
    const Fields::ENUMDATA axisVel = faceVelocityFields[axis];
    const intType iCellBound = mesh.nCells(axis) - 1;   // Index of cell at the boundary
    const intType iFaceBound = iCellBound + 1;          // Index of face at the boundary

    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p].chip(iCellBound, axis) += faceVelocities[axisVel].chip(iFaceBound, axis) 
                                              * faceVelocities[axisVel].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::uniform:
            boundaryConstants[boundaryPatch]  = faceVelocities[axisVel].chip(iFaceBound, axis)
                                              * faceVelocities[axisVel].chip(iFaceBound, axis).constant( boundaryConditionStructs[boundaryPatch].value * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::extrapolated:
            coeffs[p   ].chip(iCellBound, axis) += faceVelocities[axisVel].chip(iFaceBound, axis) 
                                                 * faceVelocities[axisVel].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].p * mesh.cellLengthsInv[axis](iCellBound) );
            coeffs[west].chip(iCellBound, axis) += faceVelocities[axisVel].chip(iFaceBound, axis) 
                                                 * faceVelocities[axisVel].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].a * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        default:
            break;
    }

}


void AdvectionNegativeBoundary( ArrayAllocator<TransportCoefficients, array3D> &coeffs, 
                                EnumVector<BoundaryPatches, array2D> &boundaryConstants,
                                const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionStruct > &boundaryConditionStructs,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    static const std::array<Fields::ENUMDATA, 3> faceVelocityFields = {F::U, F::V, F::W}; // Used to get corresponding velocity field from axis
    
    const BoundaryPatches::ENUMDATA boundaryPatch = NegativePatch[axis];
    const TransportCoefficients::ENUMDATA east = HiCoeff[axis];
    const Fields::ENUMDATA axisVel = faceVelocityFields[axis];
    const intType iCellBound = 0;   // Index of cell at the boundary 
    const intType iFaceBound = 0;   // Index of face at the boundary

    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p].chip(iCellBound, axis) += - faceVelocities[axisVel].chip(iFaceBound, axis) 
                                              *   faceVelocities[axisVel].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::uniform:
            boundaryConstants[boundaryPatch]  = - faceVelocities[axisVel].chip(iFaceBound, axis)
                                              *   faceVelocities[axisVel].chip(iFaceBound, axis).constant( boundaryConditionStructs[boundaryPatch].value * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::extrapolated:
            coeffs[p   ].chip(iCellBound, axis) += - faceVelocities[axisVel].chip(iFaceBound, axis) 
                                                 *   faceVelocities[axisVel].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].p * mesh.cellLengthsInv[axis](iCellBound) );
            coeffs[east].chip(iCellBound, axis) += - faceVelocities[axisVel].chip(iFaceBound, axis) 
                                                 *   faceVelocities[axisVel].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].a * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        default:
            break;
    }

}


void SetAdvectionCoefficients( ArrayAllocator<TransportCoefficients, array3D> &coeffs, 
                               EnumVector<BoundaryPatches, array2D> &boundaryConstants,
                               const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                               const Mesh &mesh, 
                               const InputData &inputData, 
                               const Fields::ENUMDATA field)
{
    using enum TransportCoefficients::ENUMDATA;

    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;
    
    // For now assumes that all coefficiencients are set to zero

    // Upwind internal faces
    EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {

        // Upwind internal faces
        Upwind(coeffs, faceVelocities, mesh, axis);

        // Boundary faces
        AdvectionPositiveBoundary(coeffs, boundaryConstants, faceVelocities, mesh, boundaryConditions[field], axis);
        AdvectionNegativeBoundary(coeffs, boundaryConstants, faceVelocities, mesh, boundaryConditions[field], axis);

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
                                    const EnumVector< BoundaryPatches, InputData::BoundaryConditionStruct > &boundaryConditionStructs,
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
                                    const EnumVector< BoundaryPatches, InputData::BoundaryConditionStruct > &boundaryConditionStructs, 
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
                                      const InputData &inputData, 
                                      const Fields::ENUMDATA field, 
                                      const Axis::ENUMDATA axis)
{ 
    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;

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
    InterpolationPositiveBoundary(coeffs, boundaryConstants, mesh, boundaryConditions[field], axis);
    InterpolationNegativeBoundary(coeffs, boundaryConstants, mesh, boundaryConditions[field], axis);
    
    // Multiply by inverse of cell length
    for (intType i = 0; i != mesh.nCells(axis); i++) {
        coeffs[p   ](i) *= mesh.cellLengthsInv[axis](i);
        coeffs[east](i) *= mesh.cellLengthsInv[axis](i);
        coeffs[west](i) *= mesh.cellLengthsInv[axis](i); 
    }
    boundaryConstants[ PositivePatch[axis] ] *= mesh.cellLengthsInv[axis]( mesh.nCells(axis)-1 );
    boundaryConstants[ NegativePatch[axis] ] *= mesh.cellLengthsInv[axis]( 0 );

    // Divide pressure terms by density
    if (field == F::P) {
        coeffs[p   ] /= coeffs[p   ].constant( inputData.rho );
        coeffs[east] /= coeffs[west].constant( inputData.rho );
        coeffs[west] /= coeffs[west].constant( inputData.rho );

        boundaryConstants[ PositivePatch[axis] ] /= inputData.rho;
        boundaryConstants[ NegativePatch[axis] ] /= inputData.rho;
    }
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


// Face velocity correction coefficients for a single face from Momentum Weighted Coefficients. In order of westmost to eastmost.
std::array<floatType, 4> MWICoeffs( const arrayIndex3D &idx, 
                                    const array3D &AUUpInv, 
                                    const ArrayAllocator<TransportCoefficients, array1D> &AUP,
                                    const Mesh &mesh, 
                                    const floatType rho, 
                                    const Axis::ENUMDATA axis)
{
    using enum TransportCoefficients::ENUMDATA;
    using enum Axis::ENUMDATA;
    std::array<floatType, 4> coeffs{ 0.0f };
    const TransportCoefficients::ENUMDATA east = HiCoeff[axis],
                                          west = LoCoeff[axis];

    const intType i = idx[ axis ];
    const floatType d = MWIWeightingCoeff(idx, AUUpInv, axis);

    // These coefficients assume that the momentum equations have been divided through by the cell volume
    coeffs[0] = d * (1 - mesh.interpFactors[axis](i))   * AUP[west](i-1);

    coeffs[1] = d * ( (1 - mesh.interpFactors[axis](i)) * AUP[p   ](i-1)
                    +  mesh.interpFactors[axis](i)      * AUP[west](i  )
                    +  mesh.cellCenterDiffInv[axis](i)  / rho );

    coeffs[2] = d * ( (1 - mesh.interpFactors[axis](i)) * AUP[east](i-1)
                    + mesh.interpFactors[axis](i)       * AUP[p   ](i  )
                    - mesh.cellCenterDiffInv[axis](i)   / rho );

    coeffs[3] = d * mesh.interpFactors[axis](i) * AUP[east](i);

    return coeffs;
}


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


    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D HiIndex = { i, j, k },
                             LoIndex = { i, j, k };
                LoIndex[axis] -= 1;

                // Coefficients vector, in order of westmost to east most
                std::array<floatType, 4> coeffs = MWICoeffs({i, j, k}, momentumDiagCoeffInv, momentumPressureCoeffs, mesh, rho, axis); 

                // Cell on west side 
                floatType LoCellLengthInv = mesh.cellLengthsInv[axis]( LoIndex[axis] );
                continuityPressureCoeffs[west ](LoIndex) += coeffs[0] * LoCellLengthInv;
                continuityPressureCoeffs[p    ](LoIndex) += coeffs[1] * LoCellLengthInv;
                continuityPressureCoeffs[east ](LoIndex) += coeffs[2] * LoCellLengthInv;
                continuityPressureCoeffs[eeast](LoIndex) += coeffs[3] * LoCellLengthInv;

                // Cell on east side
                floatType HiCellLengthInv = mesh.cellLengthsInv[axis]( HiIndex[axis] );
                continuityPressureCoeffs[wwest](HiIndex) -= coeffs[0] * HiCellLengthInv;
                continuityPressureCoeffs[west ](HiIndex) -= coeffs[1] * HiCellLengthInv;
                continuityPressureCoeffs[p    ](HiIndex) -= coeffs[2] * HiCellLengthInv;
                continuityPressureCoeffs[east ](HiIndex) -= coeffs[3] * HiCellLengthInv;

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

        array3D                                        *momentumDiagCoeffInv     = nullptr;
        ArrayAllocator<TransportCoefficients, array1D> *momentumPressureCoeffs   = nullptr;
        EnumVector<BoundaryPatches, floatType>         *momentumBoundaryPressure = nullptr;
        switch ( axis ) 
        {
            case Axis::X:
                momentumDiagCoeffInv     = &fvCoeffs.Umom.diagCoeffInv;
                momentumPressureCoeffs   = &fvCoeffs.Umom.AP;
                momentumBoundaryPressure = &fvCoeffs.Umom.boundaryP;
                break;

            case Axis::Y:
                momentumDiagCoeffInv     = &fvCoeffs.Vmom.diagCoeffInv;
                momentumPressureCoeffs   = &fvCoeffs.Vmom.AP;
                momentumBoundaryPressure = &fvCoeffs.Vmom.boundaryP;
                break;

            case Axis::Z:
                momentumDiagCoeffInv     = &fvCoeffs.Wmom.diagCoeffInv;
                momentumPressureCoeffs   = &fvCoeffs.Wmom.AP;
                momentumBoundaryPressure = &fvCoeffs.Wmom.boundaryP;
                break;
        }

        // Internal faces
        MWInterpolationFace(fvCoeffs.Cont.AP, *momentumDiagCoeffInv, *momentumPressureCoeffs, mesh, rho, axis);

        // Boundary constants
        MWInterpolationBoundary(fvCoeffs.Cont.boundaryP, *momentumBoundaryPressure, *momentumDiagCoeffInv, mesh, axis);

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
                                         const ArrayAllocator<Fields, array3D> &faceVelocities, 
                                         const InputData &inputData)
{
    using F = Fields::ENUMDATA;
    using A = Axis::ENUMDATA;

    // Default construct the coefficients class
    FVCoefficients fvCoeffs(mesh.nCells);

    // Diffusion coefficients
    SetDiffusionCoeffients(fvCoeffs.Umom.diff, fvCoeffs.Umom.boundaryDiff, mesh, inputData, F::U);
    SetDiffusionCoeffients(fvCoeffs.Vmom.diff, fvCoeffs.Vmom.boundaryDiff, mesh, inputData, F::V);
    SetDiffusionCoeffients(fvCoeffs.Wmom.diff, fvCoeffs.Wmom.boundaryDiff, mesh, inputData, F::W);

    // Momentum advection terms
    SetAdvectionCoefficients(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, faceVelocities, mesh, inputData, F::U);
    SetAdvectionCoefficients(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, faceVelocities, mesh, inputData, F::V);
    SetAdvectionCoefficients(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, faceVelocities, mesh, inputData, F::W);

    // Add diffusion to the velocity coefficients in momentum equations
    AddDiffusion(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, fvCoeffs.Umom.diff, fvCoeffs.Umom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, fvCoeffs.Vmom.diff, fvCoeffs.Vmom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, fvCoeffs.Wmom.diff, fvCoeffs.Wmom.boundaryDiff, mesh);

    // Inverse of AP coefficient
    using TC = TransportCoefficients::ENUMDATA;
    fvCoeffs.Umom.diagCoeffInv = fvCoeffs.Umom.AU[TC::p].inverse();
    fvCoeffs.Vmom.diagCoeffInv = fvCoeffs.Vmom.AV[TC::p].inverse();
    fvCoeffs.Wmom.diagCoeffInv = fvCoeffs.Wmom.AW[TC::p].inverse();

    // Momentum pressure terms
    SetFaceInterpolatedCoefficients(fvCoeffs.Umom.AP, fvCoeffs.Umom.boundaryP, mesh, inputData, F::P, A::X);
    SetFaceInterpolatedCoefficients(fvCoeffs.Vmom.AP, fvCoeffs.Vmom.boundaryP, mesh, inputData, F::P, A::Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.Wmom.AP, fvCoeffs.Wmom.boundaryP, mesh, inputData, F::P, A::Z);

    // Continuity velocity terms
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AU, fvCoeffs.Cont.boundaryVel, mesh, inputData, F::U, A::X);
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AV, fvCoeffs.Cont.boundaryVel, mesh, inputData, F::V, A::Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AW, fvCoeffs.Cont.boundaryVel, mesh, inputData, F::W, A::Z);

    // Continuity pressure terms (from momentum weighted interpolation)
    SetMomentumInterpolationCoefficients(fvCoeffs, mesh, inputData);
    
    // Add boundary constants to source terms
    AddMomentumBoundaryConstants(fvCoeffs.Umom);
    AddMomentumBoundaryConstants(fvCoeffs.Vmom);
    AddMomentumBoundaryConstants(fvCoeffs.Wmom);
    AddContinuityBoundaryConstants(fvCoeffs.Cont);

    // Set implicit under relaxation
    fvCoeffs.Umom.relaxation = inputData.schemes.implicitRelaxation[F::U];
    fvCoeffs.Vmom.relaxation = inputData.schemes.implicitRelaxation[F::V];
    fvCoeffs.Wmom.relaxation = inputData.schemes.implicitRelaxation[F::W];
    fvCoeffs.Cont.relaxation = inputData.schemes.implicitRelaxation[F::P];

    return fvCoeffs;
}


// Update linearisation in momenum and continuity equations
void UpdateFVCoefficients(FVCoefficients &fvCoeffs, 
                          const Mesh &mesh,
                          const ArrayAllocator<Fields, CFD::array3D> &faceVelocities,
                          const InputData &inputData)
{
    using F = Fields::ENUMDATA;


    // Zero the momentum coefficients and the boundary constants
    EnumFor<TransportCoefficients>( [&] (TransportCoefficients::ENUMDATA tc) {
        if ( fvCoeffs.Umom.AU.get(tc) )
            fvCoeffs.Umom.AU[tc].setZero();

        if ( fvCoeffs.Vmom.AV.get(tc) )
            fvCoeffs.Vmom.AV[tc].setZero();

        if ( fvCoeffs.Wmom.AW.get(tc) )
            fvCoeffs.Wmom.AW[tc].setZero();


        if ( fvCoeffs.Cont.AP.get(tc) )
            fvCoeffs.Cont.AP[tc].setZero();
    } );

    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {

        fvCoeffs.Umom.boundaryVel[bp].setZero();
        fvCoeffs.Vmom.boundaryVel[bp].setZero();
        fvCoeffs.Wmom.boundaryVel[bp].setZero();
        fvCoeffs.Cont.boundaryP[bp].setZero();

    } );

    fvCoeffs.Umom.B.setZero();
    fvCoeffs.Vmom.B.setZero();
    fvCoeffs.Wmom.B.setZero();
    fvCoeffs.Cont.B.setZero();


    // Set the advection terms
    SetAdvectionCoefficients(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, faceVelocities, mesh, inputData, F::U);
    SetAdvectionCoefficients(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, faceVelocities, mesh, inputData, F::V);
    SetAdvectionCoefficients(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, faceVelocities, mesh, inputData, F::W);

    // Add in the diffusion
    AddDiffusion(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, fvCoeffs.Umom.diff, fvCoeffs.Umom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, fvCoeffs.Vmom.diff, fvCoeffs.Vmom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, fvCoeffs.Wmom.diff, fvCoeffs.Wmom.boundaryDiff, mesh);

    // Inverse of AP coefficient
    using TC = TransportCoefficients::ENUMDATA;
    fvCoeffs.Umom.diagCoeffInv = fvCoeffs.Umom.AU[TC::p].inverse();
    fvCoeffs.Vmom.diagCoeffInv = fvCoeffs.Vmom.AV[TC::p].inverse();
    fvCoeffs.Wmom.diagCoeffInv = fvCoeffs.Wmom.AW[TC::p].inverse();

    // Set the momentum interpolation coefficients
    SetMomentumInterpolationCoefficients(fvCoeffs, mesh, inputData);

    // Add in the boundary constants to the source terms
    AddMomentumBoundaryConstants(fvCoeffs.Umom);
    AddMomentumBoundaryConstants(fvCoeffs.Vmom);
    AddMomentumBoundaryConstants(fvCoeffs.Wmom);
    AddContinuityBoundaryConstants(fvCoeffs.Cont);

}


}   // end namespace CFD
