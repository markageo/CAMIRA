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

    const Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[boundaryPatch];

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
void DiffusionPositiveBoundary( EnumVector< Axis,  EnumVector<TransportCoefficients, array1D> > &diff, 
                                EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs,
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
void DiffusionNegativeBoundary( EnumVector< Axis, EnumVector<TransportCoefficients, array1D> > &diff, 
                                EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs,
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
void SetDiffusionCoeffients(EnumVector< Axis, EnumVector<TransportCoefficients, array1D> > &diff, 
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
        boundaryConstants[ LUT::PositivePatch[axis] ] *= mesh.cellLengthsInv[axis]( mesh.nCells(axis)-1 );
        boundaryConstants[ LUT::NegativePatch[axis] ] *= mesh.cellLengthsInv[axis]( 0 );

        // Multiply by viscosity
        for (intType i = 0; i != mesh.nCells(axis); i++) {
            diff[axis][p   ](i) *= inputData.nu;
            diff[axis][east](i) *= inputData.nu;
            diff[axis][west](i) *= inputData.nu;
        }
        boundaryConstants[ LUT::PositivePatch[axis] ] *= inputData.nu;
        boundaryConstants[ LUT::NegativePatch[axis] ] *= inputData.nu;

    } );
}





/*---------------------------------------------------------------------------------------------------------------*\
                                         Momentum Advection Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


// Upwind coefficients
void UpwindInterior( EnumVector<CFD::TransportCoefficients, CFD::array3D> &coeffs, 
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

    TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis], 
                                    west = LUT::LoCoeff[axis];

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


void AdvectionPositiveBoundary( EnumVector<TransportCoefficients, array3D> &coeffs, 
                                EnumVector<BoundaryPatches, array2D> &boundaryConstants,
                                const EnumVector<Axis, array3D> &faceFluxes, 
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::PositivePatch[axis];
    const TransportCoefficients::ENUMDATA west = LUT::LoCoeff[axis];
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


void AdvectionNegativeBoundary( EnumVector<TransportCoefficients, array3D> &coeffs, 
                                EnumVector<BoundaryPatches, array2D> &boundaryConstants,
                                const EnumVector<Axis, CFD::array3D> &faceFluxes, 
                                const Mesh &mesh,  
                                const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs,
                                const Axis::ENUMDATA axis)
{
    using BC = BoundaryConditions::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    const BoundaryPatches::ENUMDATA boundaryPatch = LUT::NegativePatch[axis];
    const TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis];
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


void SetAdvectionCoefficients( EnumVector<TransportCoefficients, array3D> &coeffs, 
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
        UpwindInterior(coeffs, faceFluxes, mesh, axis);

        // Boundary faces
        AdvectionPositiveBoundary(coeffs, boundaryConstants, faceFluxes, mesh, boundaryConditions.U[velocityComponent], axis);
        AdvectionNegativeBoundary(coeffs, boundaryConstants, faceFluxes, mesh, boundaryConditions.U[velocityComponent], axis);

    } );

}





/*---------------------------------------------------------------------------------------------------------------*\
                                           Add Diffusion Coefficients
\*---------------------------------------------------------------------------------------------------------------*/

void AddDiffusion( EnumVector< TransportCoefficients, array3D > &velCoeffs, 
                   EnumVector< BoundaryPatches, array2D > &boundaryVel,
                   const EnumVector< Axis, EnumVector<TransportCoefficients, array1D> > &diffCoeffs, 
                   const EnumVector< BoundaryPatches, floatType > &boundaryDiff,
                   const Mesh &mesh)
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    for (intType k = 0; k != mesh.nCells(Z); k++) {

        floatType zpk = diffCoeffs[Z][p](k),
                  ztk = diffCoeffs[Z][t](k),
                  zbk = diffCoeffs[Z][b](k);

        for (intType j = 0; j != mesh.nCells(Y); j++) {

            floatType ypj = diffCoeffs[Y][p](j),
                      ynj = diffCoeffs[Y][n](j),
                      ysj = diffCoeffs[Y][s](j);

            #pragma GCC ivdep
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

        boundaryVel[patch] += boundaryVel[patch].constant( boundaryDiff[patch] );
        
    } );

}





/*---------------------------------------------------------------------------------------------------------------*\
                                        Linear Interpolated Coefficients
\*---------------------------------------------------------------------------------------------------------------*/


void InterpolationPositiveBoundary( EnumVector< TransportCoefficients, array1D > &coeffs, 
                                    EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                    const Mesh &mesh,  
                                    const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs,
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


void InterpolationNegativeBoundary( EnumVector< TransportCoefficients, array1D > &coeffs, 
                                    EnumVector< BoundaryPatches, floatType > &boundaryConstants,
                                    const Mesh &mesh,  
                                    const EnumVector< BoundaryPatches, InputData::BoundaryConditionData > &boundaryConditionStructs, 
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
void SetFaceInterpolatedCoefficients( EnumVector<CFD::TransportCoefficients, CFD::array1D> &coeffs, 
                                      EnumVector< BoundaryPatches, floatType > &boundaryConstants, 
                                      const Mesh &mesh, 
                                      const EnumVector< BoundaryPatches, InputData::BoundaryConditionData> &boundaryConditionStructs, 
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
    InterpolationPositiveBoundary(coeffs, boundaryConstants, mesh, boundaryConditionStructs, axis);
    InterpolationNegativeBoundary(coeffs, boundaryConstants, mesh, boundaryConditionStructs, axis);
    
    // Multiply by inverse of cell length
    for (intType i = 0; i != mesh.nCells(axis); i++) {
        coeffs[p   ](i) *= mesh.cellLengthsInv[axis](i);
        coeffs[east](i) *= mesh.cellLengthsInv[axis](i);
        coeffs[west](i) *= mesh.cellLengthsInv[axis](i); 
    }
    boundaryConstants[ LUT::PositivePatch[axis] ] *= mesh.cellLengthsInv[axis]( mesh.nCells(axis)-1 );
    boundaryConstants[ LUT::NegativePatch[axis] ] *= mesh.cellLengthsInv[axis]( 0 );
}



// Divide pressure terms in momentum equaqtions by density
void DivideMomentumPressureByDensity( EnumVector<CFD::TransportCoefficients, CFD::array1D> &coeffs, 
                                      EnumVector< BoundaryPatches, floatType > &boundaryConstants, 
                                      const floatType rho, 
                                      const Axis::ENUMDATA axis)
{ 
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    TransportCoefficients::ENUMDATA east = LUT::HiCoeff[axis],    // These are just names, they can be north, south etc.
                                    west = LUT::LoCoeff[axis];  

    
    coeffs[p   ] /= coeffs[p   ].constant( rho );
    coeffs[east] /= coeffs[west].constant( rho );
    coeffs[west] /= coeffs[west].constant( rho );

    boundaryConstants[ LUT::PositivePatch[axis] ] /= rho;
    boundaryConstants[ LUT::NegativePatch[axis] ] /= rho;
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
                              +  mesh.interpFactors[axis](i)      * momentumPressureCoeffs[west](i  ) );

        mwiSparseCoeffs[2](i) = ( (1 - mesh.interpFactors[axis](i)) * momentumPressureCoeffs[east](i-1)
                              + mesh.interpFactors[axis](i)       * momentumPressureCoeffs[p   ](i  ) );

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

        mwiCompactCoeffs[0](i) = mesh.cellCenterDiffInv[axis](i)  * rhoInv;

        mwiCompactCoeffs[1](i) = - mesh.cellCenterDiffInv[axis](i)   * rhoInv;
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
template< MomentumInterpolation MI >
void MWInterpolationInteriorImplicit( ContinuityEquation<MI> &continuityEquation,
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

    // Starting index and number of faces to iterate over
    iVector3 startIndex, nFaces;
    EnumFor<Axis>( [&] ( Axis::ENUMDATA a) {
        startIndex[a] = 0;
        nFaces[a] = mesh.nCells[a];
    } );
    startIndex[axis] += 1;


    // Cell indexing
    TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                    eeast = LUT::HiHiCoeff[axis],
                                    west  = LUT::LoCoeff[axis],
                                    wwest = LUT::LoLoCoeff[axis];

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D HiIndex = { i, j, k },
                             LoIndex = { i, j, k };
                LoIndex[axis] -= 1;

                floatType d = MWIWeightingCoeff( {i, j, k}, momentumDiagCoeffInv, axis );

                // Coefficients for westmost to eastmost cell
                intType idx = HiIndex[axis];
                floatType coeff0 = d * mwiSparseCoeffs[0](idx),
                          coeff1 = d * mwiSparseCoeffs[1](idx) + mwiCompactCoeffs[0](idx),
                          coeff2 = d * mwiSparseCoeffs[2](idx) + mwiCompactCoeffs[1](idx),
                          coeff3 = d * mwiSparseCoeffs[3](idx);

                // Cell on west side 
                floatType LoCellLengthInv = mesh.cellLengthsInv[axis]( LoIndex[axis] );
                continuityPressureCoeffs[west ](LoIndex) += coeff0 * LoCellLengthInv;
                continuityPressureCoeffs[p    ](LoIndex) += coeff1 * LoCellLengthInv;
                continuityPressureCoeffs[east ](LoIndex) += coeff2 * LoCellLengthInv;
                continuityPressureCoeffs[eeast](LoIndex) += coeff3 * LoCellLengthInv;

                // Cell on east side
                floatType HiCellLengthInv = mesh.cellLengthsInv[axis]( HiIndex[axis] );
                continuityPressureCoeffs[wwest](HiIndex) -= coeff0 * HiCellLengthInv;
                continuityPressureCoeffs[west ](HiIndex) -= coeff1 * HiCellLengthInv;
                continuityPressureCoeffs[p    ](HiIndex) -= coeff2 * HiCellLengthInv;
                continuityPressureCoeffs[east ](HiIndex) -= coeff3 * HiCellLengthInv;

            }
        }
    }

}



// Semi explicit momentum interpolation coefficient for internal faces
template< MomentumInterpolation MI > 
void MWInterpolationInteriorSemiExplicit( ContinuityEquation<MI> &continuityEquation, 
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

    // Starting index and number of faces to iterate over
    iVector3 startIndex, nFaces;
    EnumFor<Axis>( [&] ( Axis::ENUMDATA a) {
        startIndex[a] = 0;
        nFaces[a] = mesh.nCells[a];
    } );
    startIndex[axis] += 1;

    // For getting the index of a neighbouring cell
    auto NeighbourIndex = [] ( arrayIndex3D index, intType shift, Axis::ENUMDATA shiftAxis ) { 
        index[shiftAxis] += shift; 
        return index; 
    };


    // Cell indexing
    TransportCoefficients::ENUMDATA east  = LUT::HiCoeff[axis], 
                                    west  = LUT::LoCoeff[axis];

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D HiIndex = { i, j, k },
                             LoIndex = { i, j, k };
                LoIndex[axis] -= 1;
                intType idx = HiIndex[axis];

                floatType d = MWIWeightingCoeff( {i, j, k}, momentumDiagCoeffInv, axis );

                floatType LoCellLengthInv = mesh.cellLengthsInv[axis]( LoIndex[axis] ),
                          HiCellLengthInv = mesh.cellLengthsInv[axis]( HiIndex[axis] );

                // Implicit compact difference --------------------------------------------------------------------------
                
                // Coefficients for westmost to eastmost cell
                floatType coeffCompact0 = d * mwiCompactCoeffs[0](idx),
                          coeffCompact1 = d * mwiCompactCoeffs[1](idx);

                // Cell on west side 
                continuityPressureCoeffs[p    ](LoIndex) += coeffCompact0 * LoCellLengthInv;
                continuityPressureCoeffs[east ](LoIndex) += coeffCompact1 * LoCellLengthInv;

                // Cell on east side
                continuityPressureCoeffs[west ](HiIndex) -= coeffCompact0 * HiCellLengthInv;
                continuityPressureCoeffs[p    ](HiIndex) -= coeffCompact1 * HiCellLengthInv;

                // ------------------------------------------------------------------------------------------------------


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
                continuitySourceTerm(LoIndex) += ( coeffSparse0 * P( G(LoWest)  )
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



// Boundary constants that come from MWI
void MWInterpolationBoundary( EnumVector<BoundaryPatches, array2D> &continuityBoundaryPressure,
                              const EnumVector<BoundaryPatches, floatType> &momentumBoundaryPressure,
                              const array3D &momentumDiagCoeffInv, 
                              const Mesh &mesh,
                              const Axis::ENUMDATA axis )
{
    // Boundary condition contribution comes from cell face one off the boundary, since MWI correction is taken as zero
    // at the boundary.

    using enum Axis::ENUMDATA;

    BoundaryPatches::ENUMDATA positivePatch = LUT::PositivePatch[ axis ],
                              negativePatch = LUT::NegativePatch[ axis ];

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



// Set momentum interpolation coefficients
template< MomentumInterpolation MI >
void SetMomentumInterpolationCoefficients( FVCoefficients<MI> &fvCoeffs,
                                           const Mesh &mesh,
                          [[maybe_unused]] const array3D &P )
{
    // Assumes that coefficients are set to zero
    // Correction is zero at the boundary faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {


        if        constexpr ( MI == MomentumInterpolation::SemiExplicit ) {

            MWInterpolationInteriorImplicit(fvCoeffs.Cont, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);

        } else if constexpr ( MI == MomentumInterpolation::Implicit ) {

            MWInterpolationInteriorSemiExplicit(fvCoeffs.Cont, P, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);

        }

        // Boundary constants
        MWInterpolationBoundary(fvCoeffs.Cont.boundaryP, fvCoeffs.Mom[axis].boundaryP, fvCoeffs.Mom[axis].diagCoeffInv, mesh, axis);

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
        momCoeffs.B.chip( 0   , axis ) -= momCoeffs.boundaryVel[negativePatch]
                                        + momCoeffs.B.chip( 0   , axis ).constant( momCoeffs.boundaryP[negativePatch] );

        // Positive side boundary
        momCoeffs.B.chip( iEnd, axis ) -= momCoeffs.boundaryVel[positivePatch]
                                        + momCoeffs.B.chip( iEnd, axis ).constant( momCoeffs.boundaryP[positivePatch] );
    }
}


template< MomentumInterpolation MI >
void AddContinuityBoundaryConstants( ContinuityEquation<MI> &contCoeffs )
{
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;
    intType iEnd;

    // Each axis
    for (intType axis = 0; axis != Axis::count; axis++) {

        positivePatch = LUT::PositivePatch[ static_cast<size_t>( axis ) ];
        negativePatch = LUT::NegativePatch[ static_cast<size_t>( axis ) ];
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
template< MomentumInterpolation MI >
FVCoefficients<MI> InitialiseFVCoefficients( const Mesh &mesh,
                                             const FieldData<array3D> &fields,
                                             const EnumVector<Axis, array3D> &faceFluxes, 
                                             const InputData &inputData)
{
    using TC = TransportCoefficients::ENUMDATA;

    // Default construct the coefficients class
    FVCoefficients<MI> fvCoeffs(mesh.nCells);

    // Momentum equations
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

         // Diffusion coefficients
        SetDiffusionCoeffients(fvCoeffs.Mom[axis].diff, fvCoeffs.Mom[axis].boundaryDiff, mesh, inputData, axis);

        // Advection terms
        SetAdvectionCoefficients(fvCoeffs.Mom[axis].AU[axis], fvCoeffs.Mom[axis].boundaryVel, faceFluxes, mesh, inputData, axis);

        // Add diffusion to velocity terms
        AddDiffusion(fvCoeffs.Mom[axis].AU[axis], fvCoeffs.Mom[axis].boundaryVel, fvCoeffs.Mom[axis].diff, fvCoeffs.Mom[axis].boundaryDiff, mesh);

        // Inverse of AP coefficient
        fvCoeffs.Mom[axis].diagCoeffInv = fvCoeffs.Mom[axis].AU[axis][TC::p].inverse();

        // Momentum pressure terms
        SetFaceInterpolatedCoefficients(fvCoeffs.Mom[axis].AP, fvCoeffs.Mom[axis].boundaryP, mesh, inputData.boundaryConditions.P, axis);
        DivideMomentumPressureByDensity(fvCoeffs.Mom[axis].AP, fvCoeffs.Mom[axis].boundaryP, inputData.rho, axis);

        // Add boundary constants to source terms
        AddMomentumBoundaryConstants(fvCoeffs.Mom[axis]);

        // Relaxation factor
        fvCoeffs.Mom[axis].relaxation = inputData.schemes.implicitRelaxation.U[axis];

    } );


    // Continuity equation
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        // Velocity terms
        SetFaceInterpolatedCoefficients(fvCoeffs.Cont.AU[axis], fvCoeffs.Cont.boundaryVel, mesh, inputData.boundaryConditions.U[axis], axis);

        // Momentum weighted interpolation constants in the coefficients
        SetMomentumInterpolationSparseConstants(fvCoeffs.Cont.mwiSparseCoeffs[axis], fvCoeffs.Mom[axis].AP, mesh, axis);
        SetMomentumInterpolationCompactConstants(fvCoeffs.Cont.mwiCompactCoeffs[axis], inputData.rho, mesh, axis);

    } );

    // Momentum Weighted interpolation
    SetMomentumInterpolationCoefficients(fvCoeffs, mesh, fields.P);
    
    // Add boundary constants to source terms
    AddContinuityBoundaryConstants(fvCoeffs.Cont);

    // Relaxation factor
    fvCoeffs.Cont.relaxation = inputData.schemes.implicitRelaxation.P;


    return fvCoeffs;
}
template  FVCoefficients<MomentumInterpolation::Implicit> InitialiseFVCoefficients( const Mesh &, const FieldData<array3D> &, const EnumVector<Axis, array3D> &, const InputData &);
template  FVCoefficients<MomentumInterpolation::SemiExplicit> InitialiseFVCoefficients( const Mesh &, const FieldData<array3D> &, const EnumVector<Axis, array3D> &, const InputData &);





// Update linearisation in momenum and continuity equations
template< MomentumInterpolation MI >
void UpdateFVCoefficients(FVCoefficients<MI> &fvCoeffs, 
                          const Mesh &mesh,
                          const FieldData<array3D> &fields,
                          const EnumVector<Axis, array3D> &faceFluxes,
                          const InputData &inputData)
{
    using TC = TransportCoefficients::ENUMDATA;

    // Zero momentum equations
    EnumFor<Axis> ( [&] (Axis::ENUMDATA axis) {

        EnumFor<TransportCoefficients>( [&] (TransportCoefficients::ENUMDATA tc) {
                fvCoeffs.Mom[axis].AU[axis][tc].setZero();
        } );

        EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
            fvCoeffs.Mom[axis].boundaryVel[bp].setZero();
        } );

        fvCoeffs.Mom[axis].B.setZero();

    } );

    // Zero continuity equation
    EnumFor<TransportCoefficients>( [&] (TransportCoefficients::ENUMDATA tc) {
            fvCoeffs.Cont.AP[tc].setZero();
    } );

    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        fvCoeffs.Cont.boundaryP[bp].setZero();
    } );
    fvCoeffs.Cont.B.setZero();



    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        // Advection terms
        SetAdvectionCoefficients(fvCoeffs.Mom[axis].AU[axis], fvCoeffs.Mom[axis].boundaryVel, faceFluxes, mesh, inputData, axis);

        // Add diffusion to velocity terms
        AddDiffusion(fvCoeffs.Mom[axis].AU[axis], fvCoeffs.Mom[axis].boundaryVel, fvCoeffs.Mom[axis].diff, fvCoeffs.Mom[axis].boundaryDiff, mesh);

        // Inverse of AP coefficient
        fvCoeffs.Mom[axis].diagCoeffInv = fvCoeffs.Mom[axis].AU[axis][TC::p].inverse();

        // Add boundary constants to source terms
        AddMomentumBoundaryConstants(fvCoeffs.Mom[axis]);

    } );


    // Set the momentum interpolation coefficients
    SetMomentumInterpolationCoefficients(fvCoeffs, mesh, fields.P);

    // Add boundary constants to source terms
    AddContinuityBoundaryConstants(fvCoeffs.Cont);

}
template void UpdateFVCoefficients(FVCoefficients<MomentumInterpolation::Implicit> &, const Mesh &, const FieldData<array3D> &, const EnumVector<Axis, array3D> &, const InputData &);
template void UpdateFVCoefficients(FVCoefficients<MomentumInterpolation::SemiExplicit> &, const Mesh &, const FieldData<array3D> &, const EnumVector<Axis, array3D> &, const InputData &);


}   // end namespace CFD
