#include "ImmersedBoundary.h"

#include "../../Core/FVTools.h"
#include "../../Core/FVLookups.h"

namespace CAMIRA
{

namespace 
{


FieldData<floatType> GetIBFieldValues( const TensorIndex3D &cellIndex,
                                       const IBCell::SourceTermData &sourceTermData, 
                                       const FieldData<Tensor3D> &fields )
{
    using CAMIRA::FVT::G;

    // The value of velocity on the boundary, just hard code this to zero to be a solid wall
    // The pressure at the immersed boundary surface is not used
    FieldData<floatType> ibValues( 0.0f );

    // Extrapolate pressure onto the immersed boundary
    ibValues.P = sourceTermData.ibExtrapCoeff_p * fields.P( G(cellIndex) )
               + sourceTermData.ibExtrapCoeff_a * fields.P( G(sourceTermData.cellIndex_a) );

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
        faceValues.U[axis] = sourceTermData.faceReconstructionCoeff_p  * fields.U[axis]( G(cellIndex) )
                           + sourceTermData.faceReconstructionCoeff_a  * fields.U[axis]( G(sourceTermData.cellIndex_a) )
                           + sourceTermData.faceReconstructionCoeff_ib * sourceTermData.ibValues.U[axis];
    } );

    // Pressure extrapolated to face
    faceValues.P = sourceTermData.faceExtrapCoeff_p  *  fields.P( G(cellIndex) )
                 + sourceTermData.faceExtrapCoeff_a  *  fields.P( G(sourceTermData.cellIndex_a) );

    return faceValues;
}



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


}   // end anonymous namespace




void UpdateIBData( IBData &ibData, 
                   const FieldData<Tensor3D> &fields )
{
    // Go through each component, we do this so each one gets its own correction
    for ( auto &ibCellsComponent : ibData.ibCells ) {

        #pragma omp parallel for 
        for ( auto &ibCell : ibCellsComponent ) { 
            for ( auto &sourceTermData : ibCell.sourceTermsData ) {

                // Update values on the immersed boundary
                sourceTermData.ibValues = GetIBFieldValues( ibCell.cellIndex, sourceTermData, fields );

                // Use new immersed boundary values to update the face values
                sourceTermData.faceValues = ReconstructFaceValues( ibCell.cellIndex, sourceTermData, fields );

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