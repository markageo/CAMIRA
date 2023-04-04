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

BoundaryConditions::ENUMDATA GetDiffusionBC(const InputData::BoundaryConditionData &boundaryConditions, const BoundaryPatches::ENUMDATA boundaryPatch, 
                                const Fields::ENUMDATA field, const Fields::ENUMDATA fieldToCheck, Fields::ENUMDATA orthogonalField1, const Fields::ENUMDATA orthogonalField2)
{
    // Check if continuity equation implies a zero gradient boundary condition. This occurs if both orthogonal fields have a uniform BC
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
    using BP = BoundaryPatches::ENUMDATA;
    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;
    
    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;
    indexVector3 endIndexVector = { mesh.nCells(X) - 1, mesh.nCells(Y) - 1 , mesh.nCells(Z) - 1};
    intType iEnd;

    TransportCoefficients::ENUMDATA east, west;     // These are just names, they can be north, south etc.
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    Axis::ENUMDATA axis;
    BoundaryConditions::ENUMDATA positivePatchBC, negativePatchBC;  // Store these since the continuity equation can override a BC to be zeroGradient


    // Diffusion in each axis is calculated in the same way
    for (int axisNum = 0; axisNum != Axis::count; axisNum++) {

        axis = static_cast<Axis::ENUMDATA>(axisNum);
        iEnd = endIndexVector(axis);

        // Set variables unique to each axis
        if         (axis == X) {
            positivePatch = BP::xPositive;
            negativePatch = BP::xNegative;
            east = e;
            west = w;    
            positivePatchBC = GetDiffusionBC(boundaryConditions, positivePatch, field, F::U, F::V, F::W);
            negativePatchBC = GetDiffusionBC(boundaryConditions, negativePatch, field, F::U, F::V, F::W); 
              
        } else if (axis == Y) {
            positivePatch = BP::yPositive;
            negativePatch = BP::yNegative;
            east = n;
            west = s;
            positivePatchBC = GetDiffusionBC(boundaryConditions, positivePatch, field, F::V, F::W, F::U);
            negativePatchBC = GetDiffusionBC(boundaryConditions, negativePatch, field, F::V, F::W, F::U); 

        } else if (axis == Z) {
            positivePatch = BP::zPositive;
            negativePatch = BP::zNegative;
            east = t;
            west = b;
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
                boundaryConstants[positivePatch] += 2*mesh.cellLengthsInv[axis](iEnd) * boundaryConditions[field][positivePatch].value; // This is on the RHS
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
                boundaryConstants[negativePatch] += 2*mesh.cellLengthsInv[axis](0) * boundaryConditions[field][negativePatch].value; // This is on the RHS
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


void SetPicardCoefficients(ArrayAllocator<CFD::TransportCoefficients, CFD::array3D> &coeffs, const ArrayAllocator<Fields> &faceVelocities,
                           const Mesh &mesh, const InputData &inputData, const Fields::ENUMDATA field)
{
    using BC = BoundaryConditions::ENUMDATA;
    using BP = BoundaryPatches::ENUMDATA;
    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;

    // Coefficients are updated by looping through cell faces. This makes boundary conditions simpler and halves the number of upwind checks.
    // ALl directions can contribute to the p coefficient, so initialise it to zero
    coeffs[p].setConstant(0);
    
    // X normal faces
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

    
    // Y normal faces
    floatType coeff_n, coeff_s;
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


    // Z normal faces
    floatType coeff_t, coeff_b;
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


    // +x boundary
    intType iEnd = mesh.nCells[X] - 1;
    switch ( boundaryConditions[field][BP::xPositive].type ) {
        
        case BC::zeroGradient:
            coeffs[p].chip(iEnd, X) += faceVelocities[F::U].chip(iEnd+1, X) * faceVelocities[F::U].constant( mesh.cellLengthsInv[X](iEnd) );
            break;

        case BC::uniform:
            break;

        case BC::extrapolated:
            break;

        default:
            break;
    }


    // -x boundary
    switch ( boundaryConditions[field][BP::xNegative].type ) {
        
        case BC::zeroGradient:
            coeffs[p].chip(0, X) += faceVelocities[F::U].chip(0, X) * faceVelocities[F::U].constant( mesh.cellLengthsInv[X](0) );
            break;

        case BC::uniform:
            break;

        case BC::extrapolated:
            break;

        default:
            break;
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
    using BP = BoundaryPatches::ENUMDATA;
    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;

    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    TransportCoefficients::ENUMDATA east, west;     // These are just names, they can be north, south etc.
    if         (axis == X) {
        positivePatch = BP::xPositive;
        negativePatch = BP::xNegative;
        east = e;
        west = w;
    } else if (axis == Y) {
        positivePatch = BP::yPositive;
        negativePatch = BP::yNegative;
        east = n;
        west = s;
    } else if (axis == Z) {
        positivePatch = BP::zPositive;
        negativePatch = BP::zNegative;
        east = t;
        west = b;
    }
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
    SetPicardCoefficients(fvCoeffs.Umom.AU, faceVelocities, mesh, inputData, Fields::U);


    // Momentum pressure terms
    SetFaceInterpolatedCoefficients(fvCoeffs.Umom.AP, fvCoeffs.Umom.boundaryP, mesh, inputData, Fields::P, Axis::X);
    SetFaceInterpolatedCoefficients(fvCoeffs.Vmom.AP, fvCoeffs.Vmom.boundaryP, mesh, inputData, Fields::P, Axis::Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.Wmom.AP, fvCoeffs.Wmom.boundaryP, mesh, inputData, Fields::P, Axis::Z);

    // Continuity velocity terms
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AU, fvCoeffs.Cont.boundaryVel, mesh, inputData, Fields::U, Axis::X);
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AV, fvCoeffs.Cont.boundaryVel, mesh, inputData, Fields::V, Axis::Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AW, fvCoeffs.Cont.boundaryVel, mesh, inputData, Fields::W, Axis::Z);

    // Continuity pressure terms (Rhie-Chow interpolation)

    // Source terms
    /* NULL */

    // Boundary constants in source terms

}


void UpdateFVCoefficients(FVCoefficients &fvCoeffs, const Mesh &mesh, const ArrayAllocator<Fields> &faceVelocities)
{

}


}   // end namespace CFD



