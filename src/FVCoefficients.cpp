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
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    iterType iEnd = mesh.nCells(0) - 1,
             jEnd = mesh.nCells(1) - 1,
             kEnd = mesh.nCells(2) - 1;


    // ---------------------------------------------------- U momentum ---------------------------------------------------- //

    // Internal cells
    for (iterType i = 1; i != iEnd-1; i++) {
        fvCoeffs.aup[p](i) = 1 - mesh.interpFactors[X](i+1) - mesh.interpFactors[X](i);
        fvCoeffs.aup[e](i) = mesh.interpFactors[X](i+1);
        fvCoeffs.aup[w](i) = - ( 1 - mesh.interpFactors[X](i) ); 
    }

    // +x boundary
    switch ( boundaryConditions[F::P][BP::xPositive].type ) {
        
        case BC::zeroGradient:
            fvCoeffs.aup[p]( iEnd ) =  ( 1 - mesh.interpFactors[X](iEnd) );
            fvCoeffs.aup[e]( iEnd ) = 0;
            fvCoeffs.aup[w]( iEnd ) = -( 1 - mesh.interpFactors[X](iEnd) );
            break;

        case BC::uniform:
            fvCoeffs.aup[p]( iEnd ) = -mesh.interpFactors[X](iEnd);
            fvCoeffs.aup[e]( iEnd ) = 0;
            fvCoeffs.aup[w]( iEnd ) = -( 1 - mesh.interpFactors[X](iEnd) );
            fvCoeffs.bu.chip( iEnd, X ) += fvCoeffs.bu.constant( -boundaryConditions[F::P][BP::xPositive].value );
            break;

        case BC::extrapolated:
            fvCoeffs.aup[p]( iEnd ) = mesh.extrapFactors[BP::xPositive].p + mesh.interpFactors[X](iEnd);
            fvCoeffs.aup[e]( iEnd ) = 0;
            fvCoeffs.aup[w]( iEnd ) = mesh.extrapFactors[BP::xPositive].a - ( 1 - mesh.interpFactors[X](iEnd) );
            break;

        default:
            break;
    }


    // -x boundary
    switch ( boundaryConditions[F::P][BP::xNegative].type ) {
        
        case BC::zeroGradient:
            fvCoeffs.aup[p]( 0 ) = -mesh.interpFactors[X](1);
            fvCoeffs.aup[e]( 0 ) =  mesh.interpFactors[X](1);
            fvCoeffs.aup[w]( 0 ) = 0;
            break;

        case BC::uniform:
            fvCoeffs.aup[p]( 0 ) = 1 - mesh.interpFactors[X](1);
            fvCoeffs.aup[e]( 0 ) = mesh.interpFactors[X](1);
            fvCoeffs.aup[w]( 0 ) = 0;
            fvCoeffs.bu.chip( 0, X ) += fvCoeffs.bu.constant( boundaryConditions[F::P][BP::xNegative].value );
            break;

        case BC::extrapolated:
            fvCoeffs.aup[p]( 0 ) = 1 - mesh.interpFactors[X](1) - mesh.extrapFactors[BP::xNegative].p;
            fvCoeffs.aup[e]( 0 ) = mesh.interpFactors[X](1) - mesh.extrapFactors[BP::xNegative].a;
            fvCoeffs.aup[w]( 0 ) = 0;
            break;

        default:
            break;
    } 


    // ---------------------------------------------------- V momentum ---------------------------------------------------- //

    // Internal cells
    for (iterType j = 1; j != jEnd-1; j++) {
        fvCoeffs.avp[p](j) = 1 - mesh.interpFactors[Y](j+1) - mesh.interpFactors[Y](j);
        fvCoeffs.avp[n](j) = mesh.interpFactors[Y](j+1);
        fvCoeffs.avp[s](j) = - ( 1 - mesh.interpFactors[Y](j) ); 
    }

    // +y boundary
    switch ( boundaryConditions[F::P][BP::yPositive].type ) {
        
        case BC::zeroGradient:
            fvCoeffs.avp[p]( jEnd ) =  ( 1 - mesh.interpFactors[Y](jEnd) );
            fvCoeffs.avp[n]( jEnd ) = 0;
            fvCoeffs.avp[s]( jEnd ) = -( 1 - mesh.interpFactors[Y](jEnd) );
            break;

        case BC::uniform:
            fvCoeffs.avp[p]( jEnd ) = -mesh.interpFactors[Y](jEnd);
            fvCoeffs.avp[n]( jEnd ) = 0;
            fvCoeffs.avp[s]( jEnd ) = -( 1 - mesh.interpFactors[Y](jEnd) );
            fvCoeffs.bv.chip( jEnd, Y ) += fvCoeffs.bv.constant( -boundaryConditions[F::P][BP::yPositive].value );
            break;

        case BC::extrapolated:
            fvCoeffs.avp[p]( jEnd ) = mesh.extrapFactors[BP::yPositive].p + mesh.interpFactors[Y](jEnd);
            fvCoeffs.avp[n]( jEnd ) = 0;
            fvCoeffs.avp[s]( jEnd ) = mesh.extrapFactors[BP::yPositive].a - ( 1 - mesh.interpFactors[Y](jEnd) );
            break;

        default:
            break;
    }


    // -y boundary
    switch ( boundaryConditions[F::P][BP::yNegative].type ) {
        
        case BC::zeroGradient:
            fvCoeffs.avp[p]( 0 ) = -mesh.interpFactors[Y](1);
            fvCoeffs.avp[n]( 0 ) =  mesh.interpFactors[Y](1);
            fvCoeffs.avp[s]( 0 ) = 0;
            break;

        case BC::uniform:
            fvCoeffs.avp[p]( 0 ) = 1 - mesh.interpFactors[Y](1);
            fvCoeffs.avp[n]( 0 ) = mesh.interpFactors[Y](1);
            fvCoeffs.avp[s]( 0 ) = 0;
            fvCoeffs.bv.chip( 0, Y ) += fvCoeffs.bv.constant( boundaryConditions[F::P][BP::yNegative].value );
            break;

        case BC::extrapolated:
            fvCoeffs.avp[p]( 0 ) = 1 - mesh.interpFactors[Y](1) - mesh.extrapFactors[BP::yNegative].p;
            fvCoeffs.avp[n]( 0 ) = mesh.interpFactors[Y](1) - mesh.extrapFactors[BP::yNegative].a;
            fvCoeffs.avp[s]( 0 ) = 0;
            break;

        default:
            break;
    } 


    // ---------------------------------------------------- W momentum ---------------------------------------------------- //

    // Internal cells
    for (iterType k = 1; k != kEnd-1; k++) {
        fvCoeffs.awp[p](k) = 1 - mesh.interpFactors[Z](k+1) - mesh.interpFactors[Z](k);
        fvCoeffs.awp[t](k) = mesh.interpFactors[Z](k+1);
        fvCoeffs.awp[b](k) = - ( 1 - mesh.interpFactors[Z](k) ); 
    }

    // +z boundary
    switch ( boundaryConditions[F::P][BP::zPositive].type ) {
        
        case BC::zeroGradient:
            fvCoeffs.awp[p]( kEnd ) =  ( 1 - mesh.interpFactors[Z](kEnd) );
            fvCoeffs.awp[t]( kEnd ) = 0;
            fvCoeffs.awp[b]( kEnd ) = -( 1 - mesh.interpFactors[Z](kEnd) );
            break;

        case BC::uniform:
            fvCoeffs.awp[p]( kEnd ) = -mesh.interpFactors[Z](kEnd);
            fvCoeffs.awp[t]( kEnd ) = 0;
            fvCoeffs.awp[b]( kEnd ) = -( 1 - mesh.interpFactors[Z](kEnd) );
            fvCoeffs.bw.chip( kEnd, Z ) += fvCoeffs.bw.constant( -boundaryConditions[F::P][BP::zPositive].value );
            break;

        case BC::extrapolated:
            fvCoeffs.awp[p]( kEnd ) = mesh.extrapFactors[BP::zPositive].p + mesh.interpFactors[Z](kEnd);
            fvCoeffs.awp[t]( kEnd ) = 0;
            fvCoeffs.awp[b]( kEnd ) = mesh.extrapFactors[BP::zPositive].a - ( 1 - mesh.interpFactors[Z](kEnd) );
            break;

        default:
            break;
    }


    // -z boundary
    switch ( boundaryConditions[F::P][BP::zNegative].type ) {
        
        case BC::zeroGradient:
            fvCoeffs.awp[p]( 0 ) = -mesh.interpFactors[Z](1);
            fvCoeffs.awp[t]( 0 ) =  mesh.interpFactors[Z](1);
            fvCoeffs.awp[b]( 0 ) = 0;
            break;

        case BC::uniform:
            fvCoeffs.awp[p]( 0 ) = 1 - mesh.interpFactors[Z](1);
            fvCoeffs.awp[t]( 0 ) = mesh.interpFactors[Z](1);
            fvCoeffs.awp[b]( 0 ) = 0;
            fvCoeffs.bw.chip( 0, Z ) += fvCoeffs.bw.constant( boundaryConditions[F::P][BP::zNegative].value );
            break;

        case BC::extrapolated:
            fvCoeffs.awp[p]( 0 ) = 1 - mesh.interpFactors[Z](1) - mesh.extrapFactors[BP::zNegative].p;
            fvCoeffs.awp[t]( 0 ) = mesh.interpFactors[Z](1) - mesh.extrapFactors[BP::zNegative].a;
            fvCoeffs.awp[b]( 0 ) = 0;
            break;

        default:
            break;
    } 
    

}



void SetVelocityContinuity(FVCoefficients &fvCoeffs, const Mesh &mesh, const InputData::BoundaryConditionData &boundaryConditions)
{
    using F = Fields::ENUMDATA;
    using BC = BoundaryConditions::ENUMDATA;
    using BP = BoundaryPatches::ENUMDATA;
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    iterType iEnd = mesh.nCells(0) - 1,
             jEnd = mesh.nCells(1) - 1,
             kEnd = mesh.nCells(2) - 1;

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
    SetVelocityContinuity(fvCoeffs, mesh, boundaryConditions);

    // Continuity pressure terms

}


void UpdateFVCoefficients(FVCoefficients &fvCoeffs, const Mesh &mesh, const ArrayAllocator<Fields> &faceVelocities)
{

}


}   // end namespace CFD



