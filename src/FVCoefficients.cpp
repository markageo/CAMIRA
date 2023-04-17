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
BoundaryConditions::ENUMDATA GetDiffusionBC( const InputData::BoundaryConditionData &boundaryConditions, 
                                             const BoundaryPatches::ENUMDATA boundaryPatch, 
                                             const Fields::ENUMDATA field, 
                                             const Fields::ENUMDATA fieldToCheck, 
                                             const Fields::ENUMDATA orthogonalField1, 
                                             const Fields::ENUMDATA orthogonalField2)
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
void SetDiffusionCoeffients(std::vector< ArrayAllocator<TransportCoefficients, array1D> > &diff, 
                            std::vector<floatType> &boundaryConstants, 
                            const Mesh &mesh, 
                            const InputData &inputData, 
                            const Fields::ENUMDATA field)
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
void UpwindXnormal( ArrayAllocator<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
                    const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                    const Mesh &mesh)
{

    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType uf, coeff_e, coeff_w;
    for (iterType k = 0; k != faceVelocities[F::U].dimension(Z); k++) {
        for (iterType j = 0; j != faceVelocities[F::U].dimension(Y); j++) {
            coeffs[p](0, j, k) = 0;     // This one doesn't get reset in the loop
            for (iterType i = 1; i != faceVelocities[F::U].dimension(X)-1; i++) {
                
                uf = faceVelocities[F::U](i, j, k);
                coeff_w = uf * mesh.cellLengthsInv[X](i-1);
                coeff_e = uf * mesh.cellLengthsInv[X](i);

                // Cell on west side
                coeffs[e](i-1, j, k) = std::min( coeff_w, 0.0 );
                coeffs[p](i-1, j, k) += coeff_w - coeffs[e](i, j, k);

                // Cell on east side
                coeffs[w](i, j, k)  = std::max( coeff_e, 0.0 );
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
    for (iterType k = 0; k != faceVelocities[F::U].dimension(Z); k++) {
        for (iterType j = 1; j != faceVelocities[F::U].dimension(Y)-1; j++) {
            for (iterType i = 0; i != faceVelocities[F::U].dimension(X); i++) {
                
                uf = faceVelocities[F::V](i, j, k);
                coeff_s = uf * mesh.cellLengthsInv[Y](j-1);
                coeff_n = uf * mesh.cellLengthsInv[Y](j);

                // Cell on south side
                coeffs[n](i, j-1, k) = std::min( coeff_s, 0.0 );
                coeffs[p](i, j-1, k) += coeff_s - coeffs[n](i, j, k);

                // Cell on north side
                coeffs[s](i, j, k)  = std::max( coeff_n, 0.0 );
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
    for (iterType k = 1; k != faceVelocities[F::U].dimension(Z)-1; k++) {
        for (iterType j = 0; j != faceVelocities[F::U].dimension(Y); j++) {
            for (iterType i = 0; i != faceVelocities[F::U].dimension(X); i++) {
                
                uf = faceVelocities[F::W](i, j, k);
                coeff_b = uf * mesh.cellLengthsInv[Z](k-1);
                coeff_t = uf * mesh.cellLengthsInv[Z](k);

                // Cell on bottom side 
                coeffs[t](i, j, k-1) = std::min( coeff_b, 0.0 );
                coeffs[p](i, j, k-1) += coeff_b - coeffs[t](i, j, k); 

                // Cell on top side
                coeffs[b](i, j, k)  = std::max( coeff_t, 0.0 );
                coeffs[p](i, j, k)  += coeff_t - coeffs[b](i, j, k);

            }
        }
    }

}


void AdvectionBoundaryConditions( ArrayAllocator<TransportCoefficients, array3D> &coeffs, 
                                  std::vector<array2D> &boundaryConstants,
                                  const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                                  const Mesh &mesh,  
                                  const std::vector<InputData::BoundaryConditionStruct> &boundaryConditionStructs, 
                                  const BoundaryPatches::ENUMDATA boundaryPatch,
                                  const intType iCellBound,
                                  const intType iFaceBound)
{

    using BC = BoundaryConditions::ENUMDATA;
    using F  = Fields::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    static const std::array<Fields::ENUMDATA, 3> faceVelocityFields = {F::U, F::V, F::W}; // Used to get corresponding velocity field from axis
    static const std::array<TransportCoefficients::ENUMDATA, 6> interiorCoeffs = {w,  // xPositive
                                                                                  e,  // xNegative
                                                                                  s,  // yPositive
                                                                                  n,  // yNegative
                                                                                  b,  // zPositive
                                                                                  t}; // zNegative

    const Axis::ENUMDATA axis = BoundaryPatchAxis[boundaryPatch];
    const TransportCoefficients::ENUMDATA interiorCoeff = interiorCoeffs[boundaryPatch];
    const Fields::ENUMDATA axisVel = faceVelocityFields[axis];

    // Axis positive boundary
    switch ( boundaryConditionStructs[boundaryPatch].type ) {
        
        case BC::zeroGradient:
            coeffs[p].chip(iCellBound, axis) += faceVelocities[axisVel].chip(iFaceBound, axis) * faceVelocities[axisVel].chip(iFaceBound, axis).constant( mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::uniform:
            boundaryConstants[boundaryPatch] += faceVelocities[axisVel].chip(iFaceBound, axis)
                                                * boundaryConstants[boundaryPatch].constant( boundaryConditionStructs[boundaryPatch].value * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        case BC::extrapolated:
            coeffs[p   ].chip(iCellBound, axis) += faceVelocities[axisVel].chip(iFaceBound, axis) 
                                                 * faceVelocities[axisVel].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].p * mesh.cellLengthsInv[axis](iCellBound) );
            coeffs[interiorCoeff].chip(iCellBound, axis) += faceVelocities[axisVel].chip(iFaceBound, axis) 
                                                          * faceVelocities[axisVel].chip(iFaceBound, axis).constant( mesh.extrapFactors[boundaryPatch].a * mesh.cellLengthsInv[axis](iCellBound) );
            break;

        default:
            break;
    }

}


void SetAdvectionCoefficients( ArrayAllocator<TransportCoefficients, array3D> &coeffs, 
                               std::vector<array2D> &boundaryConstants,
                               const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                               const Mesh &mesh, 
                               const InputData &inputData, 
                               const Fields::ENUMDATA field)
{
    using enum TransportCoefficients::ENUMDATA;

    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;
    
    // Coefficients are updated by looping through cell faces. This makes boundary conditions simpler and halves the number of upwind checks.
    UpwindXnormal(coeffs, faceVelocities, mesh);
    UpwindYnormal(coeffs, faceVelocities, mesh);
    UpwindZnormal(coeffs, faceVelocities, mesh);


    // Boundary conditions by axis
    intType iCellBound, iFaceBound;    // Index of cell and face at boundary
    for (int axis = 0; axis != Axis::count; axis++) {

        // Axis positive
        iCellBound = mesh.nCells(axis) - 1;
        iFaceBound = iCellBound + 1;
        AdvectionBoundaryConditions(coeffs, boundaryConstants, faceVelocities, mesh, boundaryConditions[field], positivePatches[axis], iCellBound, iFaceBound);

        // Axis negative
        iCellBound = 0;
        iFaceBound = iCellBound;
        AdvectionBoundaryConditions(coeffs, boundaryConstants, faceVelocities, mesh, boundaryConditions[field], negativePatches[axis], iCellBound, iFaceBound);

    }

}



/*---------------------------------------------------------------------------------------------------------------*\
                                           Add Diffusion Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

void AddDiffusion( ArrayAllocator<TransportCoefficients, array3D> &velCoeffs, 
                   std::vector< array2D > &boundaryVel,
                   const std::vector< ArrayAllocator<TransportCoefficients, array1D> > &diffCoeffs, 
                   const std::vector< floatType > &boundaryDiff,
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
void SetFaceInterpolatedCoefficients( ArrayAllocator<CFD::TransportCoefficients, CFD::array1D> &coeffs, 
                                      std::vector<floatType> &boundaryConstants, 
                                      const Mesh &mesh, 
                                      const InputData &inputData, 
                                      const Fields::ENUMDATA field, 
                                      const Axis::ENUMDATA axis)
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
    const intType i = idx(X);
    const intType j = idx(Y);
    const intType k = idx(Z);
    const TransportCoefficients::ENUMDATA east = eastCoefficients[axis];
    const TransportCoefficients::ENUMDATA west = westCoefficients[axis];
    const floatType d = 1.0 / ( AUU[p](i, j, k) + AUU[p](i-1, j, k) );

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
    fvCoeffs.Cont.AP[w ].chip(0, X) = fvCoeffs.Cont.AP[w ].chip(0, X).constant( 0 );
    fvCoeffs.Cont.AP[p ].chip(0, X) = fvCoeffs.Cont.AP[p ].chip(0, X).constant( 0 );
    fvCoeffs.Cont.AP[e ].chip(0, X) = fvCoeffs.Cont.AP[e ].chip(0, X).constant( 0 );

    for (iterType k = 0; k != mesh.nCells(Z); k++) {
        for (iterType j = 0; j != mesh.nCells(Y); j++) {
            for (iterType i = 1; i != mesh.nCells(X)-1; i++) {

                // Coefficients vector, in order of westmost to east most
                coeffs = MWICoeffs({i, j, k}, fvCoeffs.Umom.AU, fvCoeffs.Umom.AP, mesh, rho, X); 

                // Cell on west side 
                fvCoeffs.Cont.AP[w ](i-1, j, k) += coeffs[0];
                fvCoeffs.Cont.AP[p ](i-1, j, k) += coeffs[1];
                fvCoeffs.Cont.AP[e ](i-1, j, k) += coeffs[2];
                fvCoeffs.Cont.AP[ee](i-1, j, k)  = coeffs[3];

                // Cell on east side
                fvCoeffs.Cont.AP[ww](i, j, k) = coeffs[0];
                fvCoeffs.Cont.AP[w ](i, j, k) = coeffs[1];
                fvCoeffs.Cont.AP[p ](i, j, k) = coeffs[2];  // Shouldn't be += since this is the first time it is set
                fvCoeffs.Cont.AP[e ](i, j, k) = coeffs[3];

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

    std::vector<floatType> coeffs;

    // These will get written to multiple times and need to be reset to zero
    fvCoeffs.Cont.AP[n ].chip(0, Y) = fvCoeffs.Cont.AP[n ].chip(0, Y).constant( 0 );
    fvCoeffs.Cont.AP[s ].chip(0, Y) = fvCoeffs.Cont.AP[s ].chip(0, Y).constant( 0 );

    for (iterType k = 0; k != mesh.nCells(Z); k++) {
        for (iterType j = 1; j != mesh.nCells(Y)-1; j++) {
            for (iterType i = 0; i != mesh.nCells(X); i++) {
                
                // Coefficients vector, in order of westmost to east most
                coeffs = MWICoeffs({i, j, k}, fvCoeffs.Vmom.AV, fvCoeffs.Vmom.AP, mesh, rho, Y); 

                // Cell on south
                fvCoeffs.Cont.AP[s ](i-1, j, k) += coeffs[0];
                fvCoeffs.Cont.AP[p ](i-1, j, k) += coeffs[1];
                fvCoeffs.Cont.AP[n ](i-1, j, k) += coeffs[2];
                fvCoeffs.Cont.AP[nn](i-1, j, k)  = coeffs[3];

                // Cell on north
                fvCoeffs.Cont.AP[ss](i, j, k)  = coeffs[0];
                fvCoeffs.Cont.AP[s ](i, j, k)  = coeffs[1];
                fvCoeffs.Cont.AP[p ](i, j, k) += coeffs[2];
                fvCoeffs.Cont.AP[n ](i, j, k)  = coeffs[3];

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
    fvCoeffs.Cont.AP[t ].chip(0, Z) = fvCoeffs.Cont.AP[t ].chip(0, Z).constant( 0 );
    fvCoeffs.Cont.AP[b ].chip(0, Z) = fvCoeffs.Cont.AP[b ].chip(0, Z).constant( 0 );

    for (iterType k = 1; k != mesh.nCells(Z)-1; k++) {
        for (iterType j = 0; j != mesh.nCells(Y); j++) {
            for (iterType i = 0; i != mesh.nCells(X); i++) {

                // Coefficients vector, in order of westmost to east most
                coeffs = MWICoeffs({i, j, k}, fvCoeffs.Wmom.AW, fvCoeffs.Wmom.AP, mesh, rho, Z); 

                // Cell on west side 
                fvCoeffs.Cont.AP[b ](i-1, j, k) += coeffs[0];
                fvCoeffs.Cont.AP[p ](i-1, j, k) += coeffs[1];
                fvCoeffs.Cont.AP[t ](i-1, j, k) += coeffs[2];
                fvCoeffs.Cont.AP[tt](i-1, j, k)  = coeffs[3];

                // Cell on east side
                fvCoeffs.Cont.AP[bb](i, j, k)  = coeffs[0];
                fvCoeffs.Cont.AP[b ](i, j, k)  = coeffs[1];
                fvCoeffs.Cont.AP[p ](i, j, k) += coeffs[2];
                fvCoeffs.Cont.AP[t ](i, j, k)  = coeffs[3];

            }
        }
    }

}


void SetMomentumInterpolationCoefficients( FVCoefficients &fvCoeffs, 
                                           const Mesh &mesh, 
                                           const InputData &inputData)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const floatType rho = inputData.rho;
    std::vector<floatType> coeffs;

    // Internal cells
    MWInterpolationXnormal(fvCoeffs, mesh, rho);
    MWInterpolationYnormal(fvCoeffs, mesh, rho);
    MWInterpolationZnormal(fvCoeffs, mesh, rho);

    // Boundary conditions


}
 

}   // end anonymous namespace



namespace CFD 
{


void InitialiseFVCoefficients(FVCoefficients &fvCoeffs, 
                              const Mesh &mesh, 
                              const ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
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
    // SetMomentumInterpolationCoefficients(fvCoeffs, mesh, inputData);
    

    // Source terms
    /* NULL */

    // Add boundary constants to source terms

}


void UpdateFVCoefficients(FVCoefficients &fvCoeffs, 
                          const Mesh &mesh, 
                          const ArrayAllocator<Fields, CFD::array3D> &faceVelocities)
{

}


}   // end namespace CFD



