#include "FiniteVolumeFunctions.h"

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
BoundaryConditions::ENUMDATA GetDiffusionBC(const InputData::BoundaryConditionData &boundaryConditions, const BoundaryPatches::ENUMDATA boundaryPatch, 
                                const Fields::ENUMDATA field, const Fields::ENUMDATA fieldToCheck, Fields::ENUMDATA orthogonalField1, const Fields::ENUMDATA orthogonalField2)
{
    using BC = BoundaryConditions::ENUMDATA;

    if (field == fieldToCheck) {

        if (boundaryConditions[orthogonalField1][boundaryPatch].type == BC::uniform && 
            boundaryConditions[orthogonalField2][boundaryPatch].type == BC::uniform) {
            return BC::zeroGradient;
        }
    }

    return boundaryConditions[field][boundaryPatch].type;
}


// Set diffusion coefficients for a given momentum equation
void SetDiffusionCoeffients(std::vector< ArrayAllocator<TransportCoefficients, array1D> > &diff, std::vector<floatType> &boundaryConstants, 
                            const Mesh &mesh, const InputData &inputData, const Fields::ENUMDATA field)
{

    using BC = BoundaryConditions::ENUMDATA;
    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    
    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;
    indexVector3 endIndexVector = { mesh.nCells(X) - 1, mesh.nCells(Y) - 1 , mesh.nCells(Z) - 1};
    intType iEnd;

    TransportCoefficients::ENUMDATA east, west;     // These are just names, they can be north, south etc.
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    BoundaryConditions::ENUMDATA positivePatchBC, negativePatchBC;  // Store these since the continuity equation can override a BC to be zeroGradient


    // Diffusion in each axis is calculated in the same way
    for (int axis = 0; axis != Axis::count; axis++) {

        iEnd = endIndexVector(axis);
        positivePatch = positivePatches[axis];
        negativePatch = negativePatches[axis];
        east = eastCoefficients[axis];
        west = westCoefficients[axis];     

        // Set variables unique to each axis
        if         (axis == X) {  
            positivePatchBC = GetDiffusionBC(boundaryConditions, positivePatch, field, F::U, F::V, F::W);
            negativePatchBC = GetDiffusionBC(boundaryConditions, negativePatch, field, F::U, F::V, F::W); 
              
        } else if (axis == Y) {
            positivePatchBC = GetDiffusionBC(boundaryConditions, positivePatch, field, F::V, F::W, F::U);
            negativePatchBC = GetDiffusionBC(boundaryConditions, negativePatch, field, F::V, F::W, F::U); 

        } else if (axis == Z) {
            positivePatchBC = GetDiffusionBC(boundaryConditions, positivePatch, field, F::W, F::U, F::V);
            negativePatchBC = GetDiffusionBC(boundaryConditions, negativePatch, field, F::W, F::U, F::V); 

        }


        // Internal cells
        for (iterType i = 1; i != iEnd; i++) {
            diff[axis][p   ](i) = ( mesh.cellCenterDiffInv[axis](i+1) + mesh.cellCenterDiffInv[axis](i) );
            diff[axis][east](i) = - mesh.cellCenterDiffInv[axis](i+1);
            diff[axis][west](i) = - mesh.cellCenterDiffInv[axis](i);
        }


        // Axis positive boundary
        switch ( positivePatchBC ) {
            
            case BC::zeroGradient: 
                diff[axis][p   ](iEnd) = mesh.cellCenterDiffInv[axis](iEnd);
                diff[axis][east](iEnd) = 0;
                diff[axis][west](iEnd) = - mesh.cellCenterDiffInv[axis](iEnd);
                break;

            case BC::uniform:
                diff[axis][p   ](iEnd) = ( 2*mesh.cellLengthsInv[axis](iEnd) + mesh.cellCenterDiffInv[axis](iEnd) );
                diff[axis][east](iEnd) = 0;
                diff[axis][west](iEnd) = - mesh.cellCenterDiffInv[axis](iEnd);
                boundaryConstants[positivePatch] += -2*mesh.cellLengthsInv[axis](iEnd) * boundaryConditions[field][positivePatch].value;
                break;

            case BC::extrapolated:
                diff[axis][p   ](iEnd) = - 2*mesh.cellLengthsInv[axis](iEnd) * (mesh.extrapFactors[positivePatch].p - 1)  -  mesh.cellCenterDiffInv[axis](iEnd);
                diff[axis][east](iEnd) = 0;
                diff[axis][west](iEnd) = - 2*mesh.cellLengthsInv[axis](iEnd) * mesh.extrapFactors[positivePatch].a  +  mesh.cellCenterDiffInv[axis](iEnd);
                break;

            default:
                break;
        }


        // Axis negative boundary
        switch ( negativePatchBC) {
            
            case BC::zeroGradient: 
                diff[axis][p   ](0) = mesh.cellCenterDiffInv[axis](1);
                diff[axis][east](0) = - mesh.cellCenterDiffInv[axis](1);
                diff[axis][west](0) = 0;
                break;

            case BC::uniform:
                diff[axis][p   ](0) = ( mesh.cellCenterDiffInv[axis](1) + 2*mesh.cellLengthsInv[axis](0) );
                diff[axis][east](0) = - mesh.cellCenterDiffInv[axis](1);
                diff[axis][west](0) = 0;
                boundaryConstants[negativePatch] += -2*mesh.cellLengthsInv[axis](0) * boundaryConditions[field][negativePatch].value;
                break;

            case BC::extrapolated:
                diff[axis][p   ](0) = - 2*mesh.cellLengthsInv[axis](0) * (mesh.extrapFactors[negativePatch].p - 1)  -  mesh.cellCenterDiffInv[axis](1);
                diff[axis][east](0) = - 2*mesh.cellLengthsInv[axis](0) * mesh.extrapFactors[negativePatch].a  +  mesh.cellCenterDiffInv[axis](1);
                diff[axis][west](0) = 0;
                break;

            default:
                break;
        }


        // Divide by inverse cell length
        for (iterType i = 0; i != iEnd+1; i++) {
            diff[axis][p   ](i) *= mesh.cellLengthsInv[axis](i);
            diff[axis][east](i) *= mesh.cellLengthsInv[axis](i);
            diff[axis][west](i) *= mesh.cellLengthsInv[axis](i);
        }

        // Multiply by viscosity
        for (iterType i = 0; i != iEnd+1; i++) {
            diff[axis][p   ](i) *= inputData.nu;
            diff[axis][east](i) *= inputData.nu;
            diff[axis][west](i) *= inputData.nu;
        }

    }

}


/*---------------------------------------------------------------------------------------------------------------*\
                                         Momentum Velocity Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

// Upwind coefficients for X normal faces
void UpwindXnormal(ArrayAllocator<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
                   const ArrayAllocator<Fields> &faceVelocities, const Mesh &mesh)
{

    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType uf, coeff_e, coeff_w;
    for (iterType k = 0; k != faceVelocities[F::U].dimension(Z); k++) {
        for (iterType j = 0; j != faceVelocities[F::U].dimension(Y); j++) {
            for (iterType i = 1; i != faceVelocities[F::U].dimension(X)-1; i++) {
                
                uf = faceVelocities[F::U](i, j, k);
                coeff_w = uf * mesh.cellLengthsInv[X](i-1);
                coeff_e = uf * mesh.cellLengthsInv[X](i);

                // Cell on west side
                coeffs[e](i-1, j, k)  = std::min( coeff_w, 0.0 );
                coeffs[p](i-1, j, k) += coeff_w - coeffs[e](i, j, k);

                // Cell on east side
                coeffs[w](i, j, k)  = std::max( coeff_e, 0.0 );
                coeffs[p](i, j, k) += coeff_e - coeffs[w](i, j, k);

            }
        }
    }
}


// Upwind coefficients Y normal faces
void UpwindYnormal(ArrayAllocator<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
                   const ArrayAllocator<Fields> &faceVelocities, const Mesh &mesh)
{

    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType uf, coeff_n, coeff_s;
    for (iterType k = 0; k != faceVelocities[F::U].dimension(Z); k++) {
        for (iterType j = 1; j != faceVelocities[F::U].dimension(Y)-1; j++) {
            for (iterType i = 0; i != faceVelocities[F::U].dimension(X); i++) {
                
                uf = faceVelocities[F::V](i, j, k);
                coeff_s = uf * mesh.cellLengthsInv[Y](j-1);
                coeff_n = uf * mesh.cellLengthsInv[Y](j);

                // Cell on south side
                coeffs[n](i, j-1, k)  = std::min( coeff_s, 0.0 );
                coeffs[p](i, j-1, k) += coeff_s - coeffs[n](i, j, k);

                // Cell on north side
                coeffs[s](i, j, k)  = std::max( coeff_n, 0.0 );
                coeffs[p](i, j, k) += coeff_n - coeffs[s](i, j, k);

            }
        }
    }
}


// Upwind coefficients for Z normal faces
void UpwindZnormal(ArrayAllocator<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
                   const ArrayAllocator<Fields> &faceVelocities, const Mesh &mesh)
{

    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType uf, coeff_t, coeff_b;
    for (iterType k = 1; k != faceVelocities[F::U].dimension(Z)-1; k++) {
        for (iterType j = 0; j != faceVelocities[F::U].dimension(Y); j++) {
            for (iterType i = 0; i != faceVelocities[F::U].dimension(X); i++) {
                
                uf = faceVelocities[F::W](i, j, k);
                coeff_b = uf * mesh.cellLengthsInv[Z](k-1);
                coeff_t = uf * mesh.cellLengthsInv[Z](k);

                // Cell on bottom side 
                coeffs[t](i, j-1, k)  = std::min( coeff_b, 0.0 );
                coeffs[p](i, j-1, k) += coeff_b - coeffs[t](i, j, k);

                // Cell on top side
                coeffs[b](i, j, k)  = std::max( coeff_t, 0.0 );
                coeffs[p](i, j, k) += coeff_t - coeffs[b](i, j, k);

            }
        }
    }
}


void SetAdvectionCoefficients(ArrayAllocator<TransportCoefficients, array3D> &coeffs, std::vector<array2D> &boundaryConstants,
                           const ArrayAllocator<Fields> &faceVelocities, const Mesh &mesh, const InputData &inputData, const Fields::ENUMDATA field)
{
    using BC = BoundaryConditions::ENUMDATA;
    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;

    
    // ALl directions can contribute to the p coefficient, so initialise it to zero
    coeffs[p].setConstant(0);
    
    // Coefficients are updated by looping through cell faces. This makes boundary conditions simpler and halves the number of upwind checks.
    UpwindXnormal(coeffs, faceVelocities, mesh);
    UpwindYnormal(coeffs, faceVelocities, mesh);
    UpwindZnormal(coeffs, faceVelocities, mesh);


    // Boundary conditions by axis
    constexpr std::array<Fields::ENUMDATA, 3> faceVelocityFields = {F::U, F::V, F::W};
    Fields::ENUMDATA axisVel;
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    TransportCoefficients::ENUMDATA east, west;
    iterType iEnd;
    for (int axis = 0; axis != Axis::count; axis++) {
        positivePatch = positivePatches[axis];
        negativePatch = negativePatches[axis];
        east = eastCoefficients[axis];
        west = westCoefficients[axis];
        axisVel = faceVelocityFields[axis];
        iEnd = mesh.nCells(axis) - 1;


        // Axis positive boundary
        switch ( boundaryConditions[field][positivePatch].type ) {
            
            case BC::zeroGradient:
                coeffs[p].chip(iEnd, axis) += faceVelocities[axisVel].chip(iEnd+1, axis) * faceVelocities[axisVel].chip(iEnd+1, axis).constant( mesh.cellLengthsInv[axis](iEnd) );
                break;

            case BC::uniform:
                boundaryConstants[positivePatch] += faceVelocities[axisVel].chip(iEnd+1, axis)
                                                  * boundaryConstants[positivePatch].constant( boundaryConditions[field][positivePatch].value * mesh.cellLengthsInv[axis](iEnd) );
                break;

            case BC::extrapolated:
                coeffs[p   ].chip(iEnd, axis) += faceVelocities[axisVel].chip(iEnd+1, axis) 
                                               * faceVelocities[axisVel].chip(iEnd+1, axis).constant( mesh.extrapFactors[positivePatch].p * mesh.cellLengthsInv[axis](iEnd) );
                coeffs[west].chip(iEnd, axis) += faceVelocities[axisVel].chip(iEnd+1, axis) 
                                               * faceVelocities[axisVel].chip(iEnd+1, axis).constant( mesh.extrapFactors[positivePatch].a * mesh.cellLengthsInv[axis](iEnd) );
                break;

            default:
                break;
        }


        // Axis negative boundary
        switch ( boundaryConditions[field][negativePatch].type ) {
            
            case BC::zeroGradient:
                coeffs[p].chip(0, axis) += faceVelocities[axisVel].chip(0, axis) * faceVelocities[axisVel].chip(0, axis).constant( mesh.cellLengthsInv[axis](0) );
                break;

            case BC::uniform:
                boundaryConstants[negativePatch] += faceVelocities[axisVel].chip(0, axis)
                                                  * boundaryConstants[negativePatch].constant( boundaryConditions[field][negativePatch].value * mesh.cellLengthsInv[axis](0) );
                break;

            case BC::extrapolated:
                coeffs[p   ].chip(0, axis) += faceVelocities[axisVel].chip(0, axis) 
                                            * faceVelocities[axisVel].chip(0, axis).constant( mesh.extrapFactors[negativePatch].p * mesh.cellLengthsInv[axis](0) );
                coeffs[east].chip(0, axis) += faceVelocities[axisVel].chip(0, axis) 
                                            * faceVelocities[axisVel].chip(0, axis).constant( mesh.extrapFactors[negativePatch].a * mesh.cellLengthsInv[axis](0) );
                break;

            default:
                break;
        }

    }

}



/*---------------------------------------------------------------------------------------------------------------*\
                                           Add Diffusion Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

void AddDiffusion(ArrayAllocator<TransportCoefficients, array3D> &velCoeffs, std::vector< array2D > &boundaryVel,
                  const std::vector< ArrayAllocator<TransportCoefficients, array1D> > &diffCoeffs, const std::vector< floatType > &boundaryDiff,
                  const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    // Velocity coefficients
    for (iterType k = 0; k != mesh.nCells(Z); k++) {
        for (iterType j = 0; j != mesh.nCells(Y); j++) {
            for (iterType i = 0; i != mesh.nCells(X); i++) {
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
    for (int patch = 0; patch != BoundaryPatches::count; patch++) {
        boundaryVel[patch] += boundaryVel[patch].constant( boundaryDiff[patch] );
    }

}



/*---------------------------------------------------------------------------------------------------------------*\
                                        Linear Interpolated Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

// Set coefficients for quantities that are intrpolated linearly onto faces.
void SetFaceInterpolatedCoefficients(ArrayAllocator<CFD::TransportCoefficients, CFD::array1D> &coeffs, std::vector<floatType> &boundaryConstants, const Mesh &mesh, 
                                     const InputData &inputData, const Fields::ENUMDATA field, Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;

    BoundaryPatches::ENUMDATA positivePatch = positivePatches[axis],
                              negativePatch = negativePatches[axis];
    TransportCoefficients::ENUMDATA east = eastCoefficients[axis],    // These are just names, they can be north, south etc.
                                    west = westCoefficients[axis];  
    iterType iEnd = mesh.nCells(axis) - 1;

    // Internal cells
    for (iterType i = 1; i != iEnd; i++) {
        coeffs[p   ](i) = 1 - mesh.interpFactors[axis](i+1) - mesh.interpFactors[axis](i);
        coeffs[east](i) = mesh.interpFactors[axis](i+1);
        coeffs[west](i) = - ( 1 - mesh.interpFactors[axis](i) ); 
    }

    // Axis positive boundary
    switch ( boundaryConditions[field][positivePatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p   ]( iEnd ) =  ( 1 - mesh.interpFactors[axis](iEnd) );
            coeffs[east]( iEnd ) = 0;
            coeffs[west]( iEnd ) = -( 1 - mesh.interpFactors[axis](iEnd) );
            break;

        case BC::uniform:
            coeffs[p   ]( iEnd ) = -mesh.interpFactors[axis](iEnd);
            coeffs[east]( iEnd ) = 0;
            coeffs[west]( iEnd ) = -( 1 - mesh.interpFactors[axis](iEnd) );
            boundaryConstants[positivePatch] += -boundaryConditions[field][positivePatch].value;
            break;

        case BC::extrapolated:
            coeffs[p   ]( iEnd ) = mesh.extrapFactors[positivePatch].p + mesh.interpFactors[axis](iEnd);
            coeffs[east]( iEnd ) = 0;
            coeffs[west]( iEnd ) = mesh.extrapFactors[positivePatch].a - ( 1 - mesh.interpFactors[axis](iEnd) );
            break;

        default:
            break;
    }


    // Axis negative boundary
    switch ( boundaryConditions[field][negativePatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p   ]( 0 ) = -mesh.interpFactors[axis](1);
            coeffs[east]( 0 ) =  mesh.interpFactors[axis](1);
            coeffs[west]( 0 ) = 0;
            break;

        case BC::uniform:
            coeffs[p   ]( 0 ) = 1 - mesh.interpFactors[axis](1);
            coeffs[east]( 0 ) = mesh.interpFactors[axis](1);
            coeffs[west]( 0 ) = 0;
            boundaryConstants[negativePatch] += boundaryConditions[field][negativePatch].value;
            break;

        case BC::extrapolated:
            coeffs[p   ]( 0 ) = 1 - mesh.interpFactors[axis](1) - mesh.extrapFactors[negativePatch].p;
            coeffs[east]( 0 ) = mesh.interpFactors[axis](1) - mesh.extrapFactors[negativePatch].a;
            coeffs[west]( 0 ) = 0;
            break;

        default:
            break;
    } 


    // Multiply by inverse of cell length
    for (iterType i = 0; i != iEnd+1; i++) {
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

void SetMomentumInterpolationCoefficients()
{

}
 
 
}   // end anonymous namespace



namespace CFD 
{


void InitialiseFVCoefficients(FVCoefficients &fvCoeffs, const Mesh &mesh, const ArrayAllocator<Fields> &faceVelocities, 
                              const InputData &inputData)
{

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

    // Continuity pressure terms (Rhie-Chow interpolation)
    SetMomentumInterpolationCoefficients();
    

    // Source terms
    /* NULL */

    // Add boundary constants to source terms

}


void UpdateFVCoefficients(FVCoefficients &fvCoeffs, const Mesh &mesh, const ArrayAllocator<Fields> &faceVelocities)
{

}


}   // end namespace CFD



