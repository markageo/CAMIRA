#include "FiniteVolumeFunctions.h"

// Implementation file for finite volume coefficient structure and update functions

 
namespace
{

using namespace CFD;

void SetPressureMomentum(FVCoefficients &fvCoeffs, const Mesh &mesh, const InputData::BoundaryConditionData &boundaryConditions)
{

    using F = Fields::ENUMDATA;
    using BC = BoundaryConditions::ENUMDATA;
    using BP = BoundaryPatches::ENUMDATA;
    using AX = Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    iterType iEnd = mesh.nCells(0) - 1,
             jEnd = mesh.nCells(1) - 1,
             kEnd = mesh.nCells(2) - 1;

    // U momentum 
    for (iterType i = 0; i != iEnd; i++) {
        fvCoeffs.aup[p](i) = 1 - mesh.interpFactors[AX::X](i) - mesh.interpFactors[AX::X](i-1);
        fvCoeffs.aup[e](i) = mesh.interpFactors[AX::X](i);
        fvCoeffs.aup[w](i) = - ( 1 - mesh.interpFactors[AX::X](i-1) ); 
    }

    // V momentum 
    for (iterType j = 0; j != jEnd; j++) {
        fvCoeffs.avp[p](j) = 1 - mesh.interpFactors[AX::Y](j) - mesh.interpFactors[AX::Y](j-1);
        fvCoeffs.avp[n](j) = mesh.interpFactors[AX::Y](j);
        fvCoeffs.avp[s](j) = - ( 1 - mesh.interpFactors[AX::Y](j-1) ); 
    }

    // W momentum 
    for (iterType k = 0; k != kEnd; k++) {
        fvCoeffs.awp[p](k) = 1 - mesh.interpFactors[AX::Z](k) - mesh.interpFactors[AX::Z](k-1);
        fvCoeffs.awp[t](k) = mesh.interpFactors[AX::Z](k);
        fvCoeffs.awp[b](k) = - ( 1 - mesh.interpFactors[AX::Z](k-1) ); 
    }


    // +x boundary
    switch ( boundaryConditions[F::P][BP::xPositive].type ) {
        
        case BC::zeroGradient:
            fvCoeffs.aup[p]( iEnd ) += fvCoeffs.aup[e]( iEnd );
            fvCoeffs.aup[e]( iEnd ) = 0;
            break;

        case BC::uniform:
            fvCoeffs.aup[p]( iEnd ) += - ( 1 - mesh.interpFactors[AX::X]( iEnd ) );
            fvCoeffs.aup[e]( iEnd ) = 0;
            fvCoeffs.bu( iEnd ) += boundaryConditions[F::P][BP::xPositive].value; 
            break;

        case BC::extrapolated:
            
            break;

        default:
            break;
    }


    // -x boundary
    switch ( boundaryConditions[F::P][BP::xNegative].type ) {
        
        case BC::zeroGradient:

            break;

        case BC::uniform:

            break;

        case BC::extrapolated:
            
            break;

        default:
            break;
    }

}


}   // end anonymous namespace



namespace CFD 
{


void InitialiseFVCoefficients(FVCoefficients &fvCoeffs, const Mesh &mesh, const ArrayAllocator<Fields> &faceVelocities, 
    const InputData::BoundaryConditionData &boundaryConditions)
{

    // Momentum velocity terms

    // Momentum pressure terms
    SetPressureMomentum(fvCoeffs, mesh, boundaryConditions);

    // Continuity velocity terms

    // Continuity pressure terms

}


void UpdateFVCoefficients(FVCoefficients &fvCoeffs, const Mesh &mesh, const ArrayAllocator<Fields> &faceVelocities)
{

}


}   // end namespace CFD



