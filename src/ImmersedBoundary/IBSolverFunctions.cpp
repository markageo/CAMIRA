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
    FieldData<floatType> ibFieldValues( 0.0f );

    // Extrapolate pressure onto the immersed boundary
    ibFieldValues.P = sourceTermData.ibExtrapFactor_p * fields.P( G(cellIndex) )
                    + sourceTermData.ibExtrapFactor_a * fields.P( G(sourceTermData.cellIndex_a) );

    return ibFieldValues;
}



FieldData<floatType> GetGhostCellValues( const TensorIndex3D &cellIndex,
                                         const IBCell::SourceTermData &sourceTermData, 
                                         const FieldData<Tensor3D> &fields )
{
    using CFD::FVT::G;

    FieldData<floatType> ghostCellValues( 0.0f );

    ForAllFieldData( [&] (intType f) {
        ghostCellValues[f] = sourceTermData.ghostExtrapCoeff_p  * fields[f]( G(cellIndex) )
                           + sourceTermData.ghostExtrapCoeff_a  * fields[f]( G(sourceTermData.cellIndex_a) )
                           + sourceTermData.ghostExtrapCoeff_ib * sourceTermData.ibFieldValues[f];
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



floatType MomentumIBSource( const Axis::ENUMDATA momentumAxis,
                            const IBCell::SourceTermData &sourceTermData, 
                            const TensorIndex3D &cellIndex,
                            const FVCoefficients &fvCoeffs ) 
{
    Axis::ENUMDATA faceNormal = sourceTermData.direction;
    TransportCoefficients::ENUMDATA coeff = ( sourceTermData.directionIndex == +1 ) ?  LUT::HiCoeff[faceNormal] : LUT::LoCoeff[faceNormal];

    // Velocity term
    floatType ibSource = - fvCoeffs.Mom[momentumAxis].AU[momentumAxis][coeff](cellIndex) * sourceTermData.ghostCellValues.U[momentumAxis];

    // Pressure term
    if ( momentumAxis == faceNormal ) {
        ibSource += - fvCoeffs.Mom[momentumAxis].AP[coeff](cellIndex[faceNormal]) * sourceTermData.ghostCellValues.P;
    }

    return ibSource;
}



floatType ContinuityIBSource( const IBCell::SourceTermData &sourceTermData, 
                              const TensorIndex3D &cellIndex,
                              const FVCoefficients &fvCoeffs ) 
{
    Axis::ENUMDATA faceNormal = sourceTermData.direction;
    TransportCoefficients::ENUMDATA coeff  = ( sourceTermData.directionIndex == +1 ) ? LUT::HiCoeff[faceNormal]   : LUT::LoCoeff[faceNormal];
    TransportCoefficients::ENUMDATA ccoeff = ( sourceTermData.directionIndex == +1 ) ? LUT::HiHiCoeff[faceNormal] : LUT::LoLoCoeff[faceNormal];

    // Divergence term
    floatType ibSource = - fvCoeffs.Cont.AU[faceNormal][coeff](cellIndex[faceNormal]) * sourceTermData.ghostCellValues.U[faceNormal];

    // Pressure terms
    ibSource += - fvCoeffs.Cont.AP[coeff ](cellIndex) * sourceTermData.ghostCellValues.P
                - fvCoeffs.Cont.AP[ccoeff](cellIndex) * sourceTermData.farPressureGhostCellValue;

    return ibSource;
}


}   // end anonymous namespace



void AddIBSourceTerms( FVCoefficients &fvCoeffs,
                       const IBData &ibData )
{

    // Iterate through each forced cell
    for ( auto &ibCell : ibData.ibCells ) { 

        TensorIndex3D cellIndex = ibCell.cellIndex;

        // A source term is added for each forced face
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            // Momentum equations
            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                fvCoeffs.Mom[axis].B( cellIndex ) += MomentumIBSource( axis, sourceTermData, cellIndex, fvCoeffs );
            } );

            // Continuity equation
            fvCoeffs.Cont.B( cellIndex ) += ContinuityIBSource( sourceTermData, cellIndex, fvCoeffs );

        }

    }

}



void SetIBFaceFluxes( EnumVector<Axis, Tensor3D> &faceFluxes,
                      const IBData &ibData,
                      const FieldData<Tensor3D> &fields ) 
{
    using CFD::FVT::G;

    for ( auto &ibCell : ibData.ibCells ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            Axis::ENUMDATA axis = sourceTermData.direction;
            TensorIndex3D faceIndex = ibCell.cellIndex;    
            faceIndex[axis] += sourceTermData.faceDirectionIndex;

            faceFluxes[axis](faceIndex) = sourceTermData.faceInterpCoeff_p * fields.U[axis]( G(ibCell.cellIndex) )
                                        + sourceTermData.faceInterpCoeff_g * sourceTermData.ghostCellValues.U[axis];
        }
    }
}



void UpdateIBData( IBData &ibData, 
                   const FieldData<Tensor3D> &fields )
{
    for ( auto &ibCell : ibData.ibCells ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            // Update values on the immersed boundary
            sourceTermData.ibFieldValues = GetIBFieldValues( ibCell.cellIndex, sourceTermData, fields );

            // Use the new immersed boundary values to update the ghost cells
            sourceTermData.ghostCellValues           = GetGhostCellValues( ibCell.cellIndex, sourceTermData, fields );
            sourceTermData.farPressureGhostCellValue = GetFarPressureGhostCellValue( ibCell.cellIndex, sourceTermData, fields );

        }
    }
}



void MaskFields( FieldData<Tensor3D> &fields, 
                 const Tensor3D &mask )
{
    Eigen::array<Eigen::Index, 3> offsets = { nGhost, nGhost, nGhost };

    ForAllFieldData( [&] (intType f) {

        Eigen::array<Eigen::Index, 3> extents = { fields[f].dimension(0) - 2*nGhost,
                                                  fields[f].dimension(1) - 2*nGhost,
                                                  fields[f].dimension(2) - 2*nGhost };

        fields[f].slice( offsets, extents ) *= mask;

    } );
    
}


}   // end namespace CFD