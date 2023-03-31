#include "FiniteVolumeFunctions.h"

// Implementation file for finite volume coefficient structure and update functions

 
namespace
{

using namespace CFD;


bool CheckDiffusionZeroGradient(const InputData::BoundaryConditionData &boundaryConditions)
{

}


// Set diffusion coefficients for a given momentum equation
void SetDiffusionCoeffients(std::vector< ArrayAllocator<TransportCoefficients, array1D> > &diff, const Mesh &mesh, const InputData &inputData, 
                            const Fields::ENUMDATA field)
{

    // !!!!!!!!!!!! NEED TO FLIP THE SIGN OF THE COEFFICIENTS !!!!!!!!!!!!

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
    bool positivePatchZeroGradient, negativePatchZeroGradient;   // Override if boundary conditions in continuity equation imply zeroGradient BC


    // Diffusion in each axis is calculated in the same way
    for (int axisNum = 0; axisNum != Axis::count; axisNum++) {

        axis = static_cast<Axis::ENUMDATA>(axisNum);
        iEnd = endIndexVector(axis);
        positivePatchZeroGradient = false;
        negativePatchZeroGradient = false;

        if         (axis == X) {
            positivePatch = BP::xPositive;
            negativePatch = BP::xNegative;
            east = e;
            west = w;

            // Check if continuity equation implies a zero gradient boundary condition
            if (field == F::U) {

                // Check if V and W boundary conditions are uniform
                if (boundaryConditions[F::V][positivePatch].type == BC::uniform && 
                    boundaryConditions[F::W][positivePatch].type == BC::uniform) {
                    positivePatchZeroGradient = true;
                }

                if (boundaryConditions[F::V][negativePatch].type == BC::uniform && 
                    boundaryConditions[F::W][negativePatch].type == BC::uniform) {
                    negativePatchZeroGradient = true;
                }
            }

        } else if (axis == Y) {
            positivePatch = BP::yPositive;
            negativePatch = BP::yNegative;
            east = n;
            west = s;

            if (field == F::V) {

                if (boundaryConditions[F::U][positivePatch].type == BC::uniform && 
                    boundaryConditions[F::W][positivePatch].type == BC::uniform) {
                    positivePatchZeroGradient = true;
                }

                if (boundaryConditions[F::U][negativePatch].type == BC::uniform && 
                    boundaryConditions[F::W][negativePatch].type == BC::uniform) {
                    negativePatchZeroGradient = true;
                }
            }

        } else if (axis == Z) {
            positivePatch = BP::zPositive;
            negativePatch = BP::zNegative;
            east = t;
            west = b;

            if (field == F::W) {

                if (boundaryConditions[F::U][positivePatch].type == BC::uniform && 
                    boundaryConditions[F::V][positivePatch].type == BC::uniform) {
                    positivePatchZeroGradient = true;
                }

                if (boundaryConditions[F::U][negativePatch].type == BC::uniform && 
                    boundaryConditions[F::V][negativePatch].type == BC::uniform) {
                    negativePatchZeroGradient = true;
                }
            }

        }


        // Internal cells
        for (iterType i = 1; i != iEnd; i++) {
            diff[axis][p   ](i) = - ( mesh.cellCenterDiffInv[axis](i+1) + mesh.cellCenterDiffInv[axis](i) );
            diff[axis][east](i) = mesh.cellCenterDiffInv[axis](i+1);
            diff[axis][west](i) = mesh.cellCenterDiffInv[axis](i);
        }


        // Axis positive boundary
        switch ( boundaryConditions[field][positivePatch].type ) {
            
            case BC::zeroGradient: 
                diff[axis][p   ](iEnd) = - mesh.cellCenterDiffInv[axis](iEnd);
                diff[axis][east](iEnd) = 0;
                diff[axis][west](iEnd) = mesh.cellCenterDiffInv[axis](iEnd);
                break;

            case BC::uniform:

                break;

            case BC::extrapolated:

                break;

            default:
                break;
        }


        // Axis negative boundary
        switch ( boundaryConditions[field][negativePatch].type ) {
            
            case BC::zeroGradient: 
                diff[axis][p   ](0) = - mesh.cellCenterDiffInv[axis](1);
                diff[axis][east](0) = mesh.cellCenterDiffInv[axis](1);
                diff[axis][west](0) = 0;
                break;

            case BC::uniform:

                break;

            case BC::extrapolated:

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


// Set coefficients for quantities that are intrpolated linearly onto faces.
void SetFaceInterpolatedCoefficients(ArrayAllocator<CFD::TransportCoefficients, CFD::array1D> &coeffs, array3D &sourceTerm, const Mesh &mesh, 
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
            sourceTerm.chip( iEnd, axis ) += sourceTerm.constant( -boundaryConditions[field][positivePatch].value );
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
            sourceTerm.chip( 0, axis ) += sourceTerm.constant( boundaryConditions[field][negativePatch].value );
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
    SetDiffusionCoeffients(fvCoeffs.diffu, mesh, inputData, Fields::U);
    SetDiffusionCoeffients(fvCoeffs.diffv, mesh, inputData, Fields::V);
    SetDiffusionCoeffients(fvCoeffs.diffw, mesh, inputData, Fields::W);

    // Momentum velocity terms


    // Momentum pressure terms
    SetFaceInterpolatedCoefficients(fvCoeffs.aup, fvCoeffs.bu, mesh, inputData, Fields::P, Axis::X);
    SetFaceInterpolatedCoefficients(fvCoeffs.avp, fvCoeffs.bv, mesh, inputData, Fields::P, Axis::Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.avp, fvCoeffs.bv, mesh, inputData, Fields::P, Axis::Z);

    // Continuity velocity terms
    SetFaceInterpolatedCoefficients(fvCoeffs.acu, fvCoeffs.bc, mesh, inputData, Fields::U, Axis::X);
    SetFaceInterpolatedCoefficients(fvCoeffs.acv, fvCoeffs.bc, mesh, inputData, Fields::V, Axis::Y);
    SetFaceInterpolatedCoefficients(fvCoeffs.acw, fvCoeffs.bc, mesh, inputData, Fields::W, Axis::Z);

    // Continuity pressure terms (Rhie-Chow interpolation)

}


void UpdateFVCoefficients(FVCoefficients &fvCoeffs, const Mesh &mesh, const ArrayAllocator<Fields> &faceVelocities)
{

}


}   // end namespace CFD



