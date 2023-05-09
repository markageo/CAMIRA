#include "FiniteVolume.h"

#include "Utils.h"

#include <algorithm>
#include <iostream>

// Implementation file for finite volume coefficient structure and update functions
 
namespace
{

using namespace CFD;

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

    Fields::ENUMDATA fieldToCheck;
    Fields::ENUMDATA orthogonalField1; 
    Fields::ENUMDATA orthogonalField2;
    const Axis::ENUMDATA axis = boundaryPatchAxis[boundaryPatch];

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

    const BoundaryPatches::ENUMDATA boundaryPatch = positivePatches[axis];
    const TransportCoefficients::ENUMDATA west = westCoefficients[axis];
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
            // These work out to zero, since for linear extrapolation, the gradients will be equal on both cell faces
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

    const BoundaryPatches::ENUMDATA boundaryPatch = negativePatches[axis];
    const TransportCoefficients::ENUMDATA east = eastCoefficients[axis];
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
            // These work out to zero, since for linear extrapolation, the gradients will be equal on both cell faces
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
    Axis::ENUMDATA axis;
    for (int a = 0; a != Axis::count; a++) {
        axis = static_cast<Axis::ENUMDATA>(a);

        positivePatch = positivePatches[axis];
        negativePatch = negativePatches[axis];
        east = eastCoefficients[axis];
        west = westCoefficients[axis];     

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
        // Axis positive boundary
        if (positivePatchBC != BC::zeroGradient) {
            DiffusionPositiveBoundary(diff, boundaryConstants, mesh, boundaryConditions[field], static_cast<Axis::ENUMDATA>(axis));
        }

        // Axis negative boundary
        if (negativePatchBC != BC::zeroGradient) {
            DiffusionNegativeBoundary(diff, boundaryConstants, mesh, boundaryConditions[field], static_cast<Axis::ENUMDATA>(axis));
        }


        // Divide by inverse cell length
        for (intType i = 0; i != mesh.nCells(axis); i++) {
            diff[axis][p   ](i) *= mesh.cellLengthsInv[axis](i);
            diff[axis][east](i) *= mesh.cellLengthsInv[axis](i);
            diff[axis][west](i) *= mesh.cellLengthsInv[axis](i);
        }


        // Multiply by viscosity
        for (intType i = 0; i != mesh.nCells(axis); i++) {
            diff[axis][p   ](i) *= inputData.nu;
            diff[axis][east](i) *= inputData.nu;
            diff[axis][west](i) *= inputData.nu;
        }

    }

}





/*---------------------------------------------------------------------------------------------------------------*\
                                         Momentum Advection Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

// Upwind coefficients for X normal faces
void UpwindXnormal( ArrayAllocator<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
                    const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                    const Mesh &mesh)
{

    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType uf, coeff_e, coeff_w;
    for (intType k = 0; k != faceVelocities[F::U].dimension(Z); k++) {
        for (intType j = 0; j != faceVelocities[F::U].dimension(Y); j++) {
            coeffs[p](0, j, k) = 0;     // This one doesn't get reset in the loop
            for (intType i = 1; i != faceVelocities[F::U].dimension(X)-1; i++) {    // Boundary condition in the x direction
                
                uf = faceVelocities[F::U](i, j, k);
                coeff_w =   uf * mesh.cellLengthsInv[X](i-1);
                coeff_e = - uf * mesh.cellLengthsInv[X](i);

                // Cell on west side
                coeffs[e](i-1, j, k) = std::min( coeff_w, static_cast<floatType>(0.0f) );
                coeffs[p](i-1, j, k) += coeff_w - coeffs[e](i, j, k);

                // Cell on east side
                coeffs[w](i, j, k)  = std::min( coeff_e, static_cast<floatType>(0.0f) );      // The sign of this coefficient is negative
                coeffs[p](i, j, k)  = coeff_e - coeffs[w](i, j, k);  // Shouldn't be += since this is the first time it is touched

            }
        }
    }

}


// Upwind coefficients Y normal faces
void UpwindYnormal( ArrayAllocator<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
                    const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                    const Mesh &mesh)
{

    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType uf, coeff_n, coeff_s;
    for (intType k = 0; k != faceVelocities[F::V].dimension(Z); k++) {
        for (intType j = 1; j != faceVelocities[F::V].dimension(Y)-1; j++) {    // Boundary condition in the y direction
            for (intType i = 0; i != faceVelocities[F::V].dimension(X); i++) {
                
                uf = faceVelocities[F::V](i, j, k);
                coeff_s =   uf * mesh.cellLengthsInv[Y](j-1);
                coeff_n = - uf * mesh.cellLengthsInv[Y](j);

                // Cell on south side
                coeffs[n](i, j-1, k) = std::min( coeff_s, static_cast<floatType>(0.0f) );
                coeffs[p](i, j-1, k) += coeff_s - coeffs[n](i, j, k);

                // Cell on north side
                coeffs[s](i, j, k)  = std::min( coeff_n, static_cast<floatType>(0.0f) );     // The sign of this coefficient is negative
                coeffs[p](i, j, k)  += coeff_n - coeffs[s](i, j, k); 

            }
        }
    }

}


// Upwind coefficients for Z normal faces
void UpwindZnormal( ArrayAllocator<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
                    const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                    const Mesh &mesh)
{

    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType uf, coeff_t, coeff_b;
    for (intType k = 1; k != faceVelocities[F::W].dimension(Z)-1; k++) {    // Boundary condition in the z direction
        for (intType j = 0; j != faceVelocities[F::W].dimension(Y); j++) {
            for (intType i = 0; i != faceVelocities[F::W].dimension(X); i++) {
                
                uf = faceVelocities[F::W](i, j, k);
                coeff_b =   uf * mesh.cellLengthsInv[Z](k-1);
                coeff_t = - uf * mesh.cellLengthsInv[Z](k);

                // Cell on bottom side 
                coeffs[t](i, j, k-1) = std::min( coeff_b, static_cast<floatType>(0.0f) );
                coeffs[p](i, j, k-1) += coeff_b - coeffs[t](i, j, k); 

                // Cell on top side
                coeffs[b](i, j, k)  = std::min( coeff_t, static_cast<floatType>(0.0f) );      // The sign of this coefficient is negative
                coeffs[p](i, j, k)  += coeff_t - coeffs[b](i, j, k);

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
    
    const BoundaryPatches::ENUMDATA boundaryPatch = positivePatches[axis];
    const TransportCoefficients::ENUMDATA west = westCoefficients[axis];
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
    
    const BoundaryPatches::ENUMDATA boundaryPatch = negativePatches[axis];
    const TransportCoefficients::ENUMDATA east = eastCoefficients[axis];
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
    
    // Upwind internal faces
    UpwindXnormal(coeffs, faceVelocities, mesh);
    UpwindYnormal(coeffs, faceVelocities, mesh);
    UpwindZnormal(coeffs, faceVelocities, mesh);

    // Boundary conditions by axis
    Axis::ENUMDATA axis;
    for (int a = 0; a != Axis::count; a++) {
        axis = static_cast<Axis::ENUMDATA>(a);
        AdvectionPositiveBoundary(coeffs, boundaryConstants, faceVelocities, mesh, boundaryConditions[field], axis);
        AdvectionNegativeBoundary(coeffs, boundaryConstants, faceVelocities, mesh, boundaryConditions[field], axis);
    }

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
    BoundaryPatches::ENUMDATA patch;
    for (int p = 0; p != BoundaryPatches::count; p++) {
        patch = static_cast<BoundaryPatches::ENUMDATA>(p);
        boundaryVel[patch] += boundaryVel[patch].constant( boundaryDiff[patch] );
    }

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

    const BoundaryPatches::ENUMDATA boundaryPatch = positivePatches[axis];
    const TransportCoefficients::ENUMDATA west = westCoefficients[axis];
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

    const BoundaryPatches::ENUMDATA boundaryPatch = negativePatches[axis];
    const TransportCoefficients::ENUMDATA east = eastCoefficients[axis];
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

    TransportCoefficients::ENUMDATA east = eastCoefficients[axis],    // These are just names, they can be north, south etc.
                                    west = westCoefficients[axis];  

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

    // Divide pressure terms by density
    if (field == F::P) {
        coeffs[p   ] /= coeffs[p   ].constant( inputData.rho );
        coeffs[east] /= coeffs[west].constant( inputData.rho );
        coeffs[west] /= coeffs[west].constant( inputData.rho );
    }
}





/*---------------------------------------------------------------------------------------------------------------*\
                        Momentum Weighted Interpolation (Rhie-Chow Interpolation) Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

// Face velocity correction coefficients for a single face from Momentum Weighted Coefficients. In order of westmost to eastmost.
std::vector<floatType> MWICoeffs( const indexVector3 &idx, 
                                  const ArrayAllocator<TransportCoefficients, array3D> &AUU, 
                                  const ArrayAllocator<TransportCoefficients, array1D> &AUP,
                                  const Mesh &mesh, 
                                  const floatType rho, 
                                  const Axis::ENUMDATA axis)
{
    using enum TransportCoefficients::ENUMDATA;
    using enum Axis::ENUMDATA;
    std::vector<floatType> coeffs(4);
    const TransportCoefficients::ENUMDATA east = eastCoefficients[axis];
    const TransportCoefficients::ENUMDATA west = westCoefficients[axis];

    // Temporary index vector that has the correct indices for the neighbouring cell of the current axis
    indexVector3 idxn( idx );
    const intType i = idx( axis );
    idxn(axis) -= 1;
    const floatType d = 1.0f/ ( AUU[p]( idx(X), idx(Y), idx(Z) )  +  AUU[p]( idxn(X), idxn(Y), idxn(Z) ) ); 

    // These coefficients assume that the momentum equations have been divided through by the cell volume
    coeffs[0] = d * (1 - mesh.interpFactors[axis](i))   * AUP[west](i-1);

    coeffs[1] = d * ( (1 - mesh.interpFactors[axis](i)) * AUP[p   ](i-1)
                    +  mesh.interpFactors[axis](i)      * AUP[west](i  )
                    +  mesh.cellCenterDiffInv[axis](i)  * rho );

    coeffs[2] = d * ( (1 - mesh.interpFactors[axis](i)) * AUP[east](i-1)
                    + mesh.interpFactors[axis](i)       * AUP[p   ](i  )
                    - mesh.cellCenterDiffInv[axis](i)   * rho );

    coeffs[3] = d * mesh.interpFactors[axis](i) * AUP[east](i);

    return coeffs;
}


void MWInterpolationXnormal( FVCoefficients &fvCoeffs, 
                             const Mesh &mesh, 
                             const floatType rho)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    std::vector<floatType> coeffs;

    // These will get written to multiple times and need to be reset to zero
    fvCoeffs.Cont.AP[w ].chip(0, X).setZero();
    fvCoeffs.Cont.AP[p ].chip(0, X).setZero();
    fvCoeffs.Cont.AP[e ].chip(0, X).setZero();

    for (intType k = 0; k != mesh.nCells(Z); k++) {
        for (intType j = 0; j != mesh.nCells(Y); j++) {
            for (intType i = 1; i != mesh.nCells(X)-1; i++) {   // Boundary condition in the x direction

                // Coefficients vector, in order of westmost to east most
                coeffs = MWICoeffs({i, j, k}, fvCoeffs.Umom.AU, fvCoeffs.Umom.AP, mesh, rho, X); 

                // Cell on west side 
                fvCoeffs.Cont.AP[w ](i-1, j, k) += coeffs[0] * mesh.cellLengthsInv[X](i-1);
                fvCoeffs.Cont.AP[p ](i-1, j, k) += coeffs[1] * mesh.cellLengthsInv[X](i-1);
                fvCoeffs.Cont.AP[e ](i-1, j, k) += coeffs[2] * mesh.cellLengthsInv[X](i-1);
                fvCoeffs.Cont.AP[ee](i-1, j, k)  = coeffs[3] * mesh.cellLengthsInv[X](i-1);

                // Cell on east side
                fvCoeffs.Cont.AP[ww](i, j, k) = coeffs[0] * mesh.cellLengthsInv[X](i);
                fvCoeffs.Cont.AP[w ](i, j, k) = coeffs[1] * mesh.cellLengthsInv[X](i);
                fvCoeffs.Cont.AP[p ](i, j, k) = coeffs[2] * mesh.cellLengthsInv[X](i);  // Shouldn't be += since this is the first time it is set
                fvCoeffs.Cont.AP[e ](i, j, k) = coeffs[3] * mesh.cellLengthsInv[X](i);

            }
        }
    }

}


void MWInterpolationYnormal( FVCoefficients &fvCoeffs, 
                             const Mesh &mesh, 
                             const floatType rho)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    std::vector<floatType> coeffs(4);

    // These will get written to multiple times and need to be reset to zero
    fvCoeffs.Cont.AP[n ].chip(0, Y).setZero();
    fvCoeffs.Cont.AP[s ].chip(0, Y).setZero();

    for (intType k = 0; k != mesh.nCells(Z); k++) {
        for (intType j = 1; j != mesh.nCells(Y)-1; j++) {   // Boundary condition in the y direction
            for (intType i = 0; i != mesh.nCells(X); i++) {
                
                // Coefficients vector, in order of westmost to east most
                coeffs = MWICoeffs({i, j, k}, fvCoeffs.Vmom.AV, fvCoeffs.Vmom.AP, mesh, rho, Y); 

                // Cell on south
                fvCoeffs.Cont.AP[s ](i, j-1, k) += coeffs[0] * mesh.cellLengthsInv[Y](j-1);
                fvCoeffs.Cont.AP[p ](i, j-1, k) += coeffs[1] * mesh.cellLengthsInv[Y](j-1);
                fvCoeffs.Cont.AP[n ](i, j-1, k) += coeffs[2] * mesh.cellLengthsInv[Y](j-1);
                fvCoeffs.Cont.AP[nn](i, j-1, k)  = coeffs[3] * mesh.cellLengthsInv[Y](j-1);

                // Cell on north
                fvCoeffs.Cont.AP[ss](i, j, k)  = coeffs[0] * mesh.cellLengthsInv[Y](j);
                fvCoeffs.Cont.AP[s ](i, j, k)  = coeffs[1] * mesh.cellLengthsInv[Y](j);
                fvCoeffs.Cont.AP[p ](i, j, k) += coeffs[2] * mesh.cellLengthsInv[Y](j);
                fvCoeffs.Cont.AP[n ](i, j, k)  = coeffs[3] * mesh.cellLengthsInv[Y](j);

            }
        }
    }

}


void MWInterpolationZnormal( FVCoefficients &fvCoeffs, 
                             const Mesh &mesh, 
                             const floatType rho)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    std::vector<floatType> coeffs;

    // These will get written to multiple times and need to be reset to zero
    fvCoeffs.Cont.AP[t ].chip(0, Z).setZero();
    fvCoeffs.Cont.AP[b ].chip(0, Z).setZero();

    for (intType k = 1; k != mesh.nCells(Z)-1; k++) {   // Boundary condition in the z direction
        for (intType j = 0; j != mesh.nCells(Y); j++) {
            for (intType i = 0; i != mesh.nCells(X); i++) {

                // Coefficients vector, in order of westmost to east most
                coeffs = MWICoeffs({i, j, k}, fvCoeffs.Wmom.AW, fvCoeffs.Wmom.AP, mesh, rho, Z); 

                // Cell on bottom side
                fvCoeffs.Cont.AP[b ](i, j, k-1) += coeffs[0] * mesh.cellLengthsInv[Z](k-1);
                fvCoeffs.Cont.AP[p ](i, j, k-1) += coeffs[1] * mesh.cellLengthsInv[Z](k-1);
                fvCoeffs.Cont.AP[t ](i, j, k-1) += coeffs[2] * mesh.cellLengthsInv[Z](k-1);
                fvCoeffs.Cont.AP[tt](i, j, k-1)  = coeffs[3] * mesh.cellLengthsInv[Z](k-1);

                // Cell on top side
                fvCoeffs.Cont.AP[bb](i, j, k)  = coeffs[0] * mesh.cellLengthsInv[Z](k);
                fvCoeffs.Cont.AP[b ](i, j, k)  = coeffs[1] * mesh.cellLengthsInv[Z](k);
                fvCoeffs.Cont.AP[p ](i, j, k) += coeffs[2] * mesh.cellLengthsInv[Z](k);
                fvCoeffs.Cont.AP[t ](i, j, k)  = coeffs[3] * mesh.cellLengthsInv[Z](k);

            }
        }
    }

}

void MWInterpolationBoundaryConstants( FVCoefficients &fvCoeffs )
{
    using BP = BoundaryPatches::ENUMDATA;

    // X axis, from U momentum
    fvCoeffs.Cont.boundaryP[BP::xPositive] = fvCoeffs.Umom.boundaryP[BP::xPositive];
    fvCoeffs.Cont.boundaryP[BP::xNegative] = fvCoeffs.Umom.boundaryP[BP::xNegative];

    // Y axis, from V momentum 
    fvCoeffs.Cont.boundaryP[BP::yPositive] = fvCoeffs.Vmom.boundaryP[BP::yPositive];
    fvCoeffs.Cont.boundaryP[BP::yNegative] = fvCoeffs.Vmom.boundaryP[BP::yNegative];

    // Z axis, from W momentum
    fvCoeffs.Cont.boundaryP[BP::zPositive] = fvCoeffs.Wmom.boundaryP[BP::zPositive];
    fvCoeffs.Cont.boundaryP[BP::zNegative] = fvCoeffs.Wmom.boundaryP[BP::zNegative];

}


void SetMomentumInterpolationCoefficients( FVCoefficients &fvCoeffs, 
                                           const Mesh &mesh, 
                                           const InputData &inputData)
{
    const floatType rho = inputData.rho;

    // Only internal cells need to be done, since correction is zero at boundary, and
    // sparse pressure gradient is taken from momentum equations, which already contain
    // the boundary condition.
    MWInterpolationXnormal(fvCoeffs, mesh, rho);
    MWInterpolationYnormal(fvCoeffs, mesh, rho);
    MWInterpolationZnormal(fvCoeffs, mesh, rho);

    // Constants that arise from boundary conditions need to be added
    MWInterpolationBoundaryConstants(fvCoeffs);
}
 




/*---------------------------------------------------------------------------------------------------------------*\
                                        Boundary Constants to Source Term
\*---------------------------------------------------------------------------------------------------------------*/

void AddMomentumBoundaryConstants( FVCoefficients::MomentumEquation &momCoeffs )
{
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    intType iEnd;

    // Each axis
    for (int axis = 0; axis != Axis::count; axis++) {

        positivePatch = positivePatches[axis];
        negativePatch = negativePatches[axis];
        iEnd = momCoeffs.B.dimension(axis)-1;

        // Negative side boundary
        momCoeffs.B.chip( 0   , axis ) += momCoeffs.boundaryVel[negativePatch]
                                        + momCoeffs.B.chip( 0   , axis ).constant( momCoeffs.boundaryP[negativePatch] );

        // Positive side boundary
        momCoeffs.B.chip( iEnd, axis ) += momCoeffs.boundaryVel[positivePatch]
                                        + momCoeffs.B.chip( iEnd, axis ).constant( momCoeffs.boundaryP[positivePatch] );
    }
}



void AddContinuityBoundaryConstants( FVCoefficients::ContinuityEquation &contCoeffs )
{
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    intType iEnd;

    // Each axis
    for (int axis = 0; axis != Axis::count; axis++) {

        positivePatch = positivePatches[axis];
        negativePatch = negativePatches[axis];
        iEnd = contCoeffs.B.dimension(axis)-1;

        // Negative side boundary
        contCoeffs.B.chip( 0   , axis ) += - contCoeffs.B.chip( 0   , axis ).constant( contCoeffs.boundaryVel[negativePatch] + contCoeffs.boundaryP[negativePatch] );

        // Positive side boundary
        contCoeffs.B.chip( iEnd, axis ) += - contCoeffs.B.chip( iEnd, axis ).constant( contCoeffs.boundaryVel[positivePatch] + contCoeffs.boundaryP[positivePatch] );
    }
}


}   // end anonymous namespace



namespace CFD 
{

// Allocate and initialise finite volume coefficients for momentum and continuity equations
FVCoefficients InitialiseFVCoefficients( const Mesh &mesh, 
                                         const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                                         const InputData &inputData)
{
    // Default construct the coefficients class
    FVCoefficients fvCoeffs(mesh.nCells);

    // Diffusion coefficients
    SetDiffusionCoeffients(fvCoeffs.Umom.diff, fvCoeffs.Umom.boundaryDiff, mesh, inputData, Fields::U);
    SetDiffusionCoeffients(fvCoeffs.Vmom.diff, fvCoeffs.Vmom.boundaryDiff, mesh, inputData, Fields::V);
    SetDiffusionCoeffients(fvCoeffs.Wmom.diff, fvCoeffs.Wmom.boundaryDiff, mesh, inputData, Fields::W);

    // Momentum velocity terms
    SetAdvectionCoefficients(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, faceVelocities, mesh, inputData, Fields::U);
    SetAdvectionCoefficients(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, faceVelocities, mesh, inputData, Fields::V);
    SetAdvectionCoefficients(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, faceVelocities, mesh, inputData, Fields::W);

    // Add diffusion to the velocity coefficients in momentum equations
    AddDiffusion(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, fvCoeffs.Umom.diff, fvCoeffs.Umom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, fvCoeffs.Vmom.diff, fvCoeffs.Vmom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, fvCoeffs.Wmom.diff, fvCoeffs.Wmom.boundaryDiff, mesh);

    // Momentum pressure terms
    SetFaceInterpolatedCoefficients(fvCoeffs.Umom.AP, fvCoeffs.Umom.boundaryP, mesh, inputData, Fields::P, Axis::X);
    SetFaceInterpolatedCoefficients(fvCoeffs.Vmom.AP, fvCoeffs.Vmom.boundaryP, mesh, inputData, Fields::P, Axis::Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.Wmom.AP, fvCoeffs.Wmom.boundaryP, mesh, inputData, Fields::P, Axis::Z);

    // Continuity velocity terms
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AU, fvCoeffs.Cont.boundaryVel, mesh, inputData, Fields::U, Axis::X);
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AV, fvCoeffs.Cont.boundaryVel, mesh, inputData, Fields::V, Axis::Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AW, fvCoeffs.Cont.boundaryVel, mesh, inputData, Fields::W, Axis::Z);

    // Continuity pressure terms (from momentum weighted interpolation)
    SetMomentumInterpolationCoefficients(fvCoeffs, mesh, inputData);
    
    // Set source terms
    /* NULL */

    // Add boundary constants to source terms
    AddMomentumBoundaryConstants(fvCoeffs.Umom);
    AddMomentumBoundaryConstants(fvCoeffs.Vmom);
    AddMomentumBoundaryConstants(fvCoeffs.Wmom);
    AddContinuityBoundaryConstants(fvCoeffs.Cont);

    return fvCoeffs;
}


// Update linearisation in momenum and continuity equations
void UpdateFVCoefficients(FVCoefficients &fvCoeffs, 
                          const Mesh &mesh, 
                          const ArrayAllocator<Fields, CFD::array3D> &faceVelocities,
                          const InputData &inputData)
{

    // Set the advection terms
    SetAdvectionCoefficients(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, faceVelocities, mesh, inputData, Fields::U);
    SetAdvectionCoefficients(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, faceVelocities, mesh, inputData, Fields::V);
    SetAdvectionCoefficients(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, faceVelocities, mesh, inputData, Fields::W);

    // Add in the diffusion
    AddDiffusion(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, fvCoeffs.Umom.diff, fvCoeffs.Umom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, fvCoeffs.Vmom.diff, fvCoeffs.Vmom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, fvCoeffs.Wmom.diff, fvCoeffs.Wmom.boundaryDiff, mesh);

    // Set the momentum interpolation coefficients
    SetMomentumInterpolationCoefficients(fvCoeffs, mesh, inputData);

    // Set the source terms, just set them back to zero for now, since there are no source terms
    fvCoeffs.Umom.B.setZero();
    fvCoeffs.Vmom.B.setZero();
    fvCoeffs.Wmom.B.setZero();
    fvCoeffs.Cont.B.setZero();

    // Add in the boundary constants
    AddMomentumBoundaryConstants(fvCoeffs.Umom);
    AddMomentumBoundaryConstants(fvCoeffs.Vmom);
    AddMomentumBoundaryConstants(fvCoeffs.Wmom);
    AddContinuityBoundaryConstants(fvCoeffs.Cont);

}


}   // end namespace CFD
