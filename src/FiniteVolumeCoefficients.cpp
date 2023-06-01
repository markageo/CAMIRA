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
    const TransportCoefficients::ENUMDATA west = NegativeCoeff[axis];
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
    const TransportCoefficients::ENUMDATA east = PositiveCoeff[axis];
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
    Axis::ENUMDATA axis;
    for (intType a = 0; a != Axis::count; a++) {
        axis = static_cast<Axis::ENUMDATA>(a);

        positivePatch = PositivePatch[axis];
        negativePatch = NegativePatch[axis];
        east = PositiveCoeff[axis];
        west = NegativeCoeff[axis];     

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
    
    const BoundaryPatches::ENUMDATA boundaryPatch = PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = NegativeCoeff[axis];
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
    const TransportCoefficients::ENUMDATA east = PositiveCoeff[axis];
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
    const TransportCoefficients::ENUMDATA west = NegativeCoeff[axis];
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
    const TransportCoefficients::ENUMDATA east = PositiveCoeff[axis];
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

    TransportCoefficients::ENUMDATA east = PositiveCoeff[axis],    // These are just names, they can be north, south etc.
                                    west = NegativeCoeff[axis];  

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
floatType MWIWeightingCoeff( const indexVector3 &idx, 
                             const ArrayAllocator<TransportCoefficients, array3D> &AUU, 
                             const Axis::ENUMDATA axis)
{
    using enum TransportCoefficients::ENUMDATA;
    using enum Axis::ENUMDATA;

    // Temporary index vector that has the correct indices for the neighbouring cell of the current axis
    indexVector3 idxn( idx );
    idxn(axis) -= 1;

    floatType d = 1.0f / ( AUU[p]( idx(X), idx(Y), idx(Z) )  +  AUU[p]( idxn(X), idxn(Y), idxn(Z) ) ); 
    return d;
}


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
    const TransportCoefficients::ENUMDATA east = PositiveCoeff[axis];
    const TransportCoefficients::ENUMDATA west = NegativeCoeff[axis];

    const intType i = idx( axis );
    const floatType d = MWIWeightingCoeff(idx, AUU, axis);

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


void MWInterpolationBoundaryXnormal( FVCoefficients &fvCoeffs,
                                     const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using BP = BoundaryPatches::ENUMDATA;

    // Boundary condition contribution comes from cell face one off the boundary, since MWI correction is taken as zero
    // at the boundary.

    floatType d;

    // X axis, from U momentum
    intType i;
    for (intType k = 0; k != mesh.nCells(Z); k++) { 
        for (intType j = 0; j != mesh.nCells(Y); j++) {

            // Negative boundary
            i = 1;  
            d = MWIWeightingCoeff( {i, j, k}, fvCoeffs.Umom.AU, X );
            fvCoeffs.Cont.boundaryP[BP::xNegative](j, k) = d * (1 - mesh.interpFactors[X](i)) * fvCoeffs.Umom.boundaryP[BP::xNegative];

            // Positive boundary
            i = mesh.nCells(X) - 1;
            d = MWIWeightingCoeff( {i, j, k}, fvCoeffs.Umom.AU, X );
            fvCoeffs.Cont.boundaryP[BP::xPositive](j, k) = d * mesh.interpFactors[X](i) * fvCoeffs.Umom.boundaryP[BP::xPositive];

        }
    }

}


void MWInterpolationBoundaryYnormal( FVCoefficients &fvCoeffs,
                                     const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using BP = BoundaryPatches::ENUMDATA;

    // Boundary condition contribution comes from cell face one off the boundary, since MWI correction is taken as zero
    // at the boundary.

    floatType d;

    // Y axis, from V momentum
    intType j;
    for (intType k = 0; k != mesh.nCells(Z); k++) { 
        for (intType i = 0; i != mesh.nCells(X); i++) {

            // Negative boundary
            j = 1;  
            d = MWIWeightingCoeff( {i, j, k}, fvCoeffs.Vmom.AV, Y );
            fvCoeffs.Cont.boundaryP[BP::yNegative](i, k) = d * (1 - mesh.interpFactors[Y](j)) * fvCoeffs.Vmom.boundaryP[BP::yNegative];

            // Positive boundary
            j = mesh.nCells(Y) - 1;
            d = MWIWeightingCoeff( {i, j, k}, fvCoeffs.Vmom.AV, Y );
            fvCoeffs.Cont.boundaryP[BP::yPositive](i, k) = d * mesh.interpFactors[Y](j) * fvCoeffs.Vmom.boundaryP[BP::yPositive];

        }
    }

}


void MWInterpolationBoundaryZnormal( FVCoefficients &fvCoeffs,
                                     const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using BP = BoundaryPatches::ENUMDATA;

    // Boundary condition contribution comes from cell face one off the boundary, since MWI correction is taken as zero
    // at the boundary.

    floatType d;

    // Z axis, from W momentum
    intType k;
    for (intType j = 0; j != mesh.nCells(Y); j++) {
        for (intType i = 0; i != mesh.nCells(X); i++) {

            // Negative boundary
            k = 1;  
            d = MWIWeightingCoeff( {i, j, k}, fvCoeffs.Wmom.AW, Z );
            fvCoeffs.Cont.boundaryP[BP::zNegative](i, j) = d * (1 - mesh.interpFactors[Z](k)) * fvCoeffs.Wmom.boundaryP[BP::zNegative];

            // Positive boundary
            k = mesh.nCells(Z) - 1;
            d = MWIWeightingCoeff( {i, j, k}, fvCoeffs.Wmom.AW, Z );
            fvCoeffs.Cont.boundaryP[BP::zPositive](i, j) = d * mesh.interpFactors[Z](k) * fvCoeffs.Wmom.boundaryP[BP::zPositive];

        }
    }

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
    MWInterpolationBoundaryXnormal(fvCoeffs, mesh);
    MWInterpolationBoundaryYnormal(fvCoeffs, mesh);
    MWInterpolationBoundaryZnormal(fvCoeffs, mesh);
}
 

/*---------------------------------------------------------------------------------------------------------------*\
                                                Implicit Relaxation
\*---------------------------------------------------------------------------------------------------------------*/


void AddRelaxation( array3D &diagonalCoeffs,
                    array3D &sourceTerms,
                    const array3D &oldField,
                    const floatType &relaxationFactor )
{
    using enum TransportCoefficients::ENUMDATA;

    // Add to the diagonal coefficient
    diagonalCoeffs *= diagonalCoeffs.constant( 1.0f / relaxationFactor );

    // Coefficients do not have ghost cells, so need to slice the field to work with it
    Eigen::array< Eigen::Index, 3 > offsets = {nGhost, nGhost, nGhost},
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





// Allocate and initialise finite volume coefficients for momentum and continuity equations
FVCoefficients InitialiseFVCoefficients( const Mesh &mesh, 
                                         const ArrayAllocator<Fields, array3D> &fieldsInitial,
                                         const ArrayAllocator<Fields, array3D> &faceVelocities, 
                                         const InputData &inputData)
{
    using TC = TransportCoefficients::ENUMDATA;
    using F = Fields::ENUMDATA;
    using A = Axis::ENUMDATA;

    // Default construct the coefficients class
    FVCoefficients fvCoeffs(mesh.nCells);

    // Diffusion coefficients
    SetDiffusionCoeffients(fvCoeffs.Umom.diff, fvCoeffs.Umom.boundaryDiff, mesh, inputData, F::U);
    SetDiffusionCoeffients(fvCoeffs.Vmom.diff, fvCoeffs.Vmom.boundaryDiff, mesh, inputData, F::V);
    SetDiffusionCoeffients(fvCoeffs.Wmom.diff, fvCoeffs.Wmom.boundaryDiff, mesh, inputData, F::W);

    // Momentum velocity terms
    // SetAdvectionCoefficients(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, faceVelocities, mesh, inputData, F::U);
    // SetAdvectionCoefficients(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, faceVelocities, mesh, inputData, F::V);
    // SetAdvectionCoefficients(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, faceVelocities, mesh, inputData, F::W);

    // Add diffusion to the velocity coefficients in momentum equations
    AddDiffusion(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, fvCoeffs.Umom.diff, fvCoeffs.Umom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, fvCoeffs.Vmom.diff, fvCoeffs.Vmom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, fvCoeffs.Wmom.diff, fvCoeffs.Wmom.boundaryDiff, mesh);

    // Momentum pressure terms
    SetFaceInterpolatedCoefficients(fvCoeffs.Umom.AP, fvCoeffs.Umom.boundaryP, mesh, inputData, F::P, A::X);
    SetFaceInterpolatedCoefficients(fvCoeffs.Vmom.AP, fvCoeffs.Vmom.boundaryP, mesh, inputData, F::P, A::Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.Wmom.AP, fvCoeffs.Wmom.boundaryP, mesh, inputData, F::P, A::Z);

    // Continuity velocity terms
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AU, fvCoeffs.Cont.boundaryVel, mesh, inputData, F::U, A::X);
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AV, fvCoeffs.Cont.boundaryVel, mesh, inputData, F::V, A::Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AW, fvCoeffs.Cont.boundaryVel, mesh, inputData, F::W, A::Z);

    // Continuity pressure terms (from momentum weighted interpolation)
    // SetMomentumInterpolationCoefficients(fvCoeffs, mesh, inputData);
    
    // Add boundary constants to source terms
    AddMomentumBoundaryConstants(fvCoeffs.Umom);
    AddMomentumBoundaryConstants(fvCoeffs.Vmom);
    AddMomentumBoundaryConstants(fvCoeffs.Wmom);
    AddContinuityBoundaryConstants(fvCoeffs.Cont);

    // Add implicit relaxation to the equations
    AddRelaxation(fvCoeffs.Umom.AU[TC::p], fvCoeffs.Umom.B, fieldsInitial[F::U], inputData.schemes.implicitRelaxation[F::U]);
    AddRelaxation(fvCoeffs.Vmom.AV[TC::p], fvCoeffs.Vmom.B, fieldsInitial[F::V], inputData.schemes.implicitRelaxation[F::V]);
    AddRelaxation(fvCoeffs.Wmom.AW[TC::p], fvCoeffs.Wmom.B, fieldsInitial[F::W], inputData.schemes.implicitRelaxation[F::W]);
    AddRelaxation(fvCoeffs.Cont.AP[TC::p], fvCoeffs.Cont.B, fieldsInitial[F::P], inputData.schemes.implicitRelaxation[F::P]);

    return fvCoeffs;
}


// Update linearisation in momenum and continuity equations
void UpdateFVCoefficients(FVCoefficients &fvCoeffs, 
                          const Mesh &mesh, 
                          const ArrayAllocator<Fields, CFD::array3D> &fieldsOld,
                          const ArrayAllocator<Fields, CFD::array3D> &faceVelocities,
                          const InputData &inputData)
{
    using TC = TransportCoefficients::ENUMDATA;
    using F = Fields::ENUMDATA;


    // Zero the momentum coefficients and the boundary constants
    // *** remove this when advection is added back in.
    EnumFor<TransportCoefficients>( [&] (TransportCoefficients::ENUMDATA tc) {
        if ( fvCoeffs.Umom.AU.get(tc) )
            fvCoeffs.Umom.AU[tc].setZero();

        if ( fvCoeffs.Vmom.AV.get(tc) )
            fvCoeffs.Vmom.AV[tc].setZero();

        if ( fvCoeffs.Wmom.AW.get(tc) )
            fvCoeffs.Wmom.AW[tc].setZero();
    } );

    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        fvCoeffs.Umom.boundaryP[bp] = 0;
        fvCoeffs.Vmom.boundaryP[bp] = 0;
        fvCoeffs.Wmom.boundaryP[bp] = 0;

        fvCoeffs.Umom.boundaryVel[bp].setZero();
        fvCoeffs.Vmom.boundaryVel[bp].setZero();
        fvCoeffs.Wmom.boundaryVel[bp].setZero();
    } );



    // // Set the advection terms
    // SetAdvectionCoefficients(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, faceVelocities, mesh, inputData, F::U);
    // SetAdvectionCoefficients(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, faceVelocities, mesh, inputData, F::V);
    // SetAdvectionCoefficients(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, faceVelocities, mesh, inputData, F::W);

    // Add in the diffusion
    AddDiffusion(fvCoeffs.Umom.AU, fvCoeffs.Umom.boundaryVel, fvCoeffs.Umom.diff, fvCoeffs.Umom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Vmom.AV, fvCoeffs.Vmom.boundaryVel, fvCoeffs.Vmom.diff, fvCoeffs.Vmom.boundaryDiff, mesh);
    AddDiffusion(fvCoeffs.Wmom.AW, fvCoeffs.Wmom.boundaryVel, fvCoeffs.Wmom.diff, fvCoeffs.Wmom.boundaryDiff, mesh);

    // Set the momentum interpolation coefficients
    // SetMomentumInterpolationCoefficients(fvCoeffs, mesh, inputData);

    // Set the source terms to zero... there may be a more efficient way to do this
    fvCoeffs.Umom.B.setZero();
    fvCoeffs.Vmom.B.setZero();
    fvCoeffs.Wmom.B.setZero();
    fvCoeffs.Cont.B.setZero();

    // Add in the boundary constants to the source terms
    AddMomentumBoundaryConstants(fvCoeffs.Umom);
    AddMomentumBoundaryConstants(fvCoeffs.Vmom);
    AddMomentumBoundaryConstants(fvCoeffs.Wmom);
    AddContinuityBoundaryConstants(fvCoeffs.Cont);

    // Add implicit relaxation to the equations
    AddRelaxation(fvCoeffs.Umom.AU[TC::p], fvCoeffs.Umom.B, fieldsOld[F::U], inputData.schemes.implicitRelaxation[F::U]);
    AddRelaxation(fvCoeffs.Vmom.AV[TC::p], fvCoeffs.Vmom.B, fieldsOld[F::V], inputData.schemes.implicitRelaxation[F::V]);
    AddRelaxation(fvCoeffs.Wmom.AW[TC::p], fvCoeffs.Wmom.B, fieldsOld[F::W], inputData.schemes.implicitRelaxation[F::W]);
    AddRelaxation(fvCoeffs.Cont.AP[TC::p], fvCoeffs.Cont.B, fieldsOld[F::P], inputData.schemes.implicitRelaxation[F::P]);
}


}   // end namespace CFD
