#include "ImmersedBoundary.h"

#include "../Tools/FVTools.h"
#include "../Tools/FVLookups.h"

namespace CFD
{

namespace 
{


FieldData<floatType> GetIBFieldValues( const TensorIndex3D &cellIndex,
                                       const IBCell::SourceTermData &sourceTermData, 
                                       const FieldData<Tensor3D> &fields )
{
    using CFD::FVT::G;

    // The value of velocity on the boundary, just hard code this to zero to be a solid wall
    FieldData<floatType> ibValues( 0.0f );

    // Extrapolate pressure onto the immersed boundary
    ibValues.P = sourceTermData.ibExtrapFactor_p * fields.P( G(cellIndex) )
               + sourceTermData.ibExtrapFactor_a * fields.P( G(sourceTermData.cellIndex_a) );

    return ibValues;
}



FieldData<floatType> ExtrapolateFaceValues( const TensorIndex3D &cellIndex,
                                            const IBCell::SourceTermData &sourceTermData, 
                                            const FieldData<Tensor3D> &fields )
{
    using CFD::FVT::G;

    FieldData<floatType> faceValues( 0.0f );

    ForAllFieldData( [&] (intType f) {
        faceValues[f] = sourceTermData.faceExtrapCoeff_p  * fields[f]( G(cellIndex) )
                      + sourceTermData.faceExtrapCoeff_a  * fields[f]( G(sourceTermData.cellIndex_a) )
                      + sourceTermData.faceExtrapCoeff_ib * sourceTermData.ibValues[f];
    } );

    return faceValues;
}



floatType CalculateVelocityFluxError( IBData &ibData )
{
    const floatType velocityFluxIB = 0.0f;
    floatType velocityFluxError = 0.0f;
    floatType velocityFlux = 0.0f;
    
    for ( auto &ibCell : ibData.ibCells ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            Axis::ENUMDATA axis = sourceTermData.direction;
            velocityFlux += sourceTermData.faceValues.U[axis] * sourceTermData.faceAreaComponent;

        }
    }

    velocityFluxError = velocityFluxIB - velocityFlux;

    return velocityFluxError;
}



void CorrectIBFaceVelocities( IBData &ibData )
{

    // Calculate the velocity flux error over the entire immersed boundary
    floatType velocityFluxError = CalculateVelocityFluxError( ibData );

    // Add the corrections to each face velocity
    for ( auto &ibCell : ibData.ibCells ) { 
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
    using CFD::FVT::G;

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
    using CFD::FVT::G;

    return sourceTermData.farPressureCoeff_p * fields.P( G(cellIndex) )
         + sourceTermData.farPressureCoeff_a * fields.P( G(sourceTermData.cellIndex_a) )
         + sourceTermData.farPressureCoeff_g * sourceTermData.ghostCellValues.P;       
}


}   // end anonymous namespace




void UpdateIBData( IBData &ibData, 
                   const FieldData<Tensor3D> &fields )
{
    for ( auto &ibCell : ibData.ibCells ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            // Update values on the immersed boundary
            sourceTermData.ibValues = GetIBFieldValues( ibCell.cellIndex, sourceTermData, fields );

            // Use new immersed boundary values to update the face values
            sourceTermData.faceValues = ExtrapolateFaceValues( ibCell.cellIndex, sourceTermData, fields );

        }
    }

    // Correct them to globally conserve mass
    CorrectIBFaceVelocities( ibData );

    for ( auto &ibCell : ibData.ibCells ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            // Extrapolate them to ghost cells
            sourceTermData.ghostCellValues           = ExtrapolateFaceToGhostCells( ibCell.cellIndex, sourceTermData, fields );
            sourceTermData.farPressureGhostCellValue = GetFarPressureGhostCellValue( ibCell.cellIndex, sourceTermData, fields );

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


}   // end namespace CFD