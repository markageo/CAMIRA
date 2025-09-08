#include "ImmersedBoundary.h"

#include "../Core/FVTools.h"
#include "../Core/FVLookups.h"

#include <cmath>

namespace CAMIRA
{

namespace 
{

// ------------------------------------------- Functions for IB without wall functions -------------------------------------------
FieldData<floatType> GetIBFieldValues( const TensorIndex3D &cellIndex,
                                       const IBCell::SourceTermData &sourceTermData, 
                                       const FieldData<Tensor3D> &fields )
{
    using CAMIRA::FVT::G;

    // The value of velocity on the boundary, just hard code this to zero to be a solid wall
    // The pressure at the immersed boundary surface is not used
    FieldData<floatType> ibValues( 0.0f );

    // Extrapolate pressure onto the immersed boundary
    ibValues.P = sourceTermData.directionalIBDataPtr->ibExtrapCoeff_p * fields.P( G(cellIndex) )
               + sourceTermData.directionalIBDataPtr->ibExtrapCoeff_a * fields.P( G(sourceTermData.cellIndex_a) );

    return ibValues;
}



FieldData<floatType> ReconstructFaceValues( const TensorIndex3D &cellIndex,
                                            const IBCell::SourceTermData &sourceTermData, 
                                            const FieldData<Tensor3D> &fields )
{
    using CAMIRA::FVT::G;

    FieldData<floatType> faceValues( 0.0f );

    // Velocity reconstruction using immersed boundary
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        faceValues.U[axis] = sourceTermData.directionalIBDataPtr->faceReconstructionCoeff_p  * fields.U[axis]( G(cellIndex) )
                           + sourceTermData.directionalIBDataPtr->faceReconstructionCoeff_a  * fields.U[axis]( G(sourceTermData.cellIndex_a) )
                           + sourceTermData.directionalIBDataPtr->faceReconstructionCoeff_ib * sourceTermData.ibValues.U[axis];
    } );

    // Pressure extrapolated to face
    faceValues.P = sourceTermData.directionalIBDataPtr->faceExtrapCoeff_p  *  fields.P( G(cellIndex) )
                 + sourceTermData.directionalIBDataPtr->faceExtrapCoeff_a  *  fields.P( G(sourceTermData.cellIndex_a) );

    return faceValues;
}

// --------------------------------------------------------------------------------------------------------------------------------



// --------------------------------------------- Functions for IB with wall functions ---------------------------------------------


void WallFunctionFaceVelocitiesAndWallShear( IBCell::SourceTermData &sourceTermData, 
                                             const FieldData<Tensor3D> &fields,
                                             const floatType nu, 
                                             const floatType rho )
{
    using enum Axis::ENUMDATA;

    const fVector3 &normalVector = sourceTermData.wallFunctionDataPtr->normalVector;
    floatType &yPlusImagePoint   = sourceTermData.wallFunctionDataPtr->yPlusImagePoint;
    floatType &yImagePoint       = sourceTermData.wallFunctionDataPtr->imagePointDistance;
    floatType &yFace             = sourceTermData.wallFunctionDataPtr->faceCenterDistance;

    // Log law constants
    const floatType kappa = 0.41f;
    const floatType CPlus = 5.0f;

    // Calculate tangential velocity at the image point
    fVector3 imagePointVelocity = { sourceTermData.wallFunctionDataPtr->fieldProbePtr->GetFieldValue( fields.U[X] ),
                                    sourceTermData.wallFunctionDataPtr->fieldProbePtr->GetFieldValue( fields.U[Y] ),
                                    sourceTermData.wallFunctionDataPtr->fieldProbePtr->GetFieldValue( fields.U[Z] ) };

    
    fVector3 imagePointWallNormalVelocity = imagePointVelocity.dot( normalVector ) * normalVector;

    fVector3 imagePointWallTangentialVelocity = imagePointVelocity - imagePointWallNormalVelocity;

    floatType imagePointWallTangentialVelocityMagnitude = imagePointWallTangentialVelocity.norm();

    floatType Rey = imagePointWallTangentialVelocityMagnitude * yImagePoint / nu;

    // Solve the log law for y+ using Newtons method
    intType maxIters = 10;
    floatType tol = 1e-5;
    for ( intType k = 0; k != maxIters; k++ ) {

        // Set the old guess
        floatType yPlusImagePointOld = yPlusImagePoint;

        // Update guess
        yPlusImagePoint = ( yPlusImagePointOld + kappa * Rey )
                        / ( std::log( yPlusImagePointOld ) + 1.0f + kappa * CPlus );

        // Check tolerence
        floatType eps = abs( yPlusImagePoint - yPlusImagePointOld );
        if ( eps < tol ) {
            break;
        }
        // std::cout << eps << std::endl;
            

    }

    // std::cout << yPlusImagePoint << ", ";

    if ( yPlusImagePoint > 400 )
        yPlusImagePoint = 400;

    if ( yPlusImagePoint < 20 )
        yPlusImagePoint = 20;

    if ( !std::isfinite( yPlusImagePoint ) )   
        yPlusImagePoint = 50;

    // yPlusImagePoint = 50;

    // std::cout << yPlusImagePoint << std::endl;

    // Calculate friction velocity
    floatType frictionVelocity = yPlusImagePoint * nu / yImagePoint;
    
    // Tangential velocity at cell face from log law
    floatType yPlusFace = yFace * frictionVelocity / nu;

    floatType uPlusFace = ( 1.0f / kappa ) * std::log( yPlusFace ) + CPlus;
    fVector3 faceWallTangentialVelocity = uPlusFace * frictionVelocity * imagePointWallTangentialVelocity.normalized();

    // Normal velocity at cell face
    fVector3 faceWallNormalVelocity = ( yFace / yImagePoint ) * imagePointWallNormalVelocity;

    // Cartesian velocity components at face
    fVector3 faceVelocity = faceWallNormalVelocity + faceWallTangentialVelocity;

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        sourceTermData.faceValues.U[axis] = faceVelocity[axis];
    } );

    // // DEBUGGING
    // EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
    //     sourceTermData.faceValues.U[axis] = ( yFace / yImagePoint ) * imagePointVelocity[axis];
    // } );
    

    // Calculate the wall shear stress magntiude
    floatType wallShearMagnitude = rho * frictionVelocity * frictionVelocity;

    // Components of tangential velocity that are tangent to the cell face
    fVector3 onFaceTangentVector = faceWallTangentialVelocity;
    onFaceTangentVector[sourceTermData.direction] = 0.0f;
    onFaceTangentVector.normalize();
    
    sourceTermData.wallShearStress = wallShearMagnitude * onFaceTangentVector;

}



void WallFunctionFacePressure( IBCell::SourceTermData &sourceTermData, 
                               const FieldData<Tensor3D> &fields )
{

    floatType imagePointPressure = sourceTermData.wallFunctionDataPtr->fieldProbePtr->GetFieldValue( fields.P );

    // Neumann condition
    sourceTermData.faceValues.P = imagePointPressure;

}



// --------------------------------------------------------------------------------------------------------------------------------




// ---------------------------------------------- Shared functions for all IB methods ---------------------------------------------

floatType CalculateVelocityFluxError( std::vector<IBCell> &ibCellsComponent )
{
    const floatType velocityFluxIB = 0.0f;
    floatType velocityFluxError = 0.0f;
    floatType velocityFlux = 0.0f;
    
    for ( auto &ibCell : ibCellsComponent ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            Axis::ENUMDATA axis = sourceTermData.direction;
            velocityFlux += sourceTermData.faceValues.U[axis] * sourceTermData.faceAreaComponent;

        }
    }

    velocityFluxError = velocityFluxIB - velocityFlux;

    return velocityFluxError;
}



void CorrectIBFaceVelocities( std::vector<IBCell> &ibCellsComponent )
{
    // Calculate the velocity flux error over the entire immersed boundary
    floatType velocityFluxError = CalculateVelocityFluxError( ibCellsComponent );

    // Add the corrections to each face velocity
    #pragma omp parallel for
    for ( auto &ibCell : ibCellsComponent ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            floatType correction = sourceTermData.velocityFluxCorrectionCoeff * velocityFluxError;
            sourceTermData.faceValues.U[ sourceTermData.direction ] += correction;

        }
    }
}



FieldData<floatType> ExtrapolateFaceToGhostCells( const TensorIndex3D &cellIndex,
                                                  const IBCell::SourceTermData &sourceTermData, 
                                                  const FieldData<Tensor3D> &fields )
{
    using CAMIRA::FVT::G;

    FieldData<floatType> ghostCellValues( 0.0f );

    ForAllFieldData( [&] (intType f) {
        ghostCellValues[f] = sourceTermData.ghostExtrapCoeff_p  * fields[f]( G(cellIndex) )
                           + sourceTermData.ghostExtrapCoeff_f  * sourceTermData.faceValues[f];         
    } );

    return ghostCellValues;
}



floatType GetFarPressureGhostCellValue( const TensorIndex3D &cellIndex,
                                        const IBCell::SourceTermData &sourceTermData,
                                        const FieldData<Tensor3D> &fields )
{
    using CAMIRA::FVT::G;

    return sourceTermData.farPressureCoeff_p * fields.P( G(cellIndex) )
         + sourceTermData.farPressureCoeff_a * fields.P( G(sourceTermData.cellIndex_a) )
         + sourceTermData.farPressureCoeff_g * sourceTermData.ghostCellValues.P;    
}

// --------------------------------------------------------------------------------------------------------------------------------


}   // end anonymous namespace




void UpdateIBData( IBData &ibData, 
                   const FieldData<Tensor3D> &fields )
{
    // Go through each component, we do this so each one gets its own correction
    for ( auto &ibCellsComponent : ibData.ibCells ) {

        if ( ibData.useWallFunctions ) {

            #pragma omp parallel for 
            for ( auto &ibCell : ibCellsComponent ) { 
                for ( auto &sourceTermData : ibCell.sourceTermsData ) {

                    // Use wall function to update velocities at cell face and wall shear stress
                    WallFunctionFaceVelocitiesAndWallShear( sourceTermData, fields, ibData.nu, ibData.rho );

                    // Set the face pressure
                    WallFunctionFacePressure( sourceTermData, fields );


                }
            }

        } else {

            #pragma omp parallel for 
            for ( auto &ibCell : ibCellsComponent ) { 
                for ( auto &sourceTermData : ibCell.sourceTermsData ) {

                    // Update values on the immersed boundary
                    sourceTermData.ibValues = GetIBFieldValues( ibCell.cellIndex, sourceTermData, fields );

                    // Use new immersed boundary values to update the face values
                    sourceTermData.faceValues = ReconstructFaceValues( ibCell.cellIndex, sourceTermData, fields );

                }
            }

        }


        // Correct them to globally conserve mass
        CorrectIBFaceVelocities( ibCellsComponent );

        #pragma omp parallel for
        for ( auto &ibCell : ibCellsComponent ) { 
            for ( auto &sourceTermData : ibCell.sourceTermsData ) {

                // Extrapolate them to ghost cells
                sourceTermData.ghostCellValues           = ExtrapolateFaceToGhostCells( ibCell.cellIndex, sourceTermData, fields );
                sourceTermData.farPressureGhostCellValue = GetFarPressureGhostCellValue( ibCell.cellIndex, sourceTermData, fields );

            }
        }

    }
}



void MaskFields( FieldData<Tensor3D> &fields, 
                 const Tensor3D &mask )
{
    ForAllFieldData( [&] (intType f) {
        fields[f] *= mask;
    } );
}


}   // end namespace CAMIRA