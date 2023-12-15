#include "ImmersedBoundary.h"

#include "../Tools/FVTools.h"
#include "../Tools/FVLookups.h"

namespace CFD
{

namespace 
{


FieldData<floatType> GetForcedFaceValues( const TensorIndex3D &cellIndex,
                                          const IBCell::FaceData &faceData, 
                                          const FieldData<Tensor3D> &fields )
{
    // Values of the field on the forced face
    FieldData<floatType> forcedFaceFieldValues;

    // The value of velocity on the boundary, just hard code this to zero to be a solid wall
    FieldData<floatType> ibFieldValues( 0.0f );

    // Interpolate velocities onto the forced face
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                
        forcedFaceFieldValues.U[axis] = faceData.interpCoeffCell * fields.U[axis](cellIndex)
                                      + faceData.interpCoeffIB   * ibFieldValues.U[axis];

    } );

    // Extrapolate pressure onto the immersed boundary
    ibFieldValues.P = faceData.extrapFactor_p * fields.P(cellIndex)
                    + faceData.extrapFactor_a * fields.P(faceData.adjacentCellIndex);

    // Interpolate the pressure onto the forced face
    forcedFaceFieldValues.P = faceData.interpCoeffCell * fields.P(cellIndex)
                            + faceData.interpCoeffIB   * ibFieldValues.P;

    return forcedFaceFieldValues;
}



floatType LoSideMomentumIBSource( const Axis::ENUMDATA momentumAxis,
                                  const IBCell::FaceData &faceData, 
                                  const TensorIndex3D &forcedCellIndex,
                                  const EnumVector<Axis, Tensor3D> &faceFluxes,
                                  const FieldData<Tensor3D> &fields,
                                  const Mesh &mesh,
                                  const InputData &inputData ) 
{
    floatType ibSource = 0.0f;

    intType faceIndexOffset   = faceData.faceIndexOffset;
    Axis::ENUMDATA faceNormal = faceData.faceNormal;
    TensorIndex3D forcedFaceIndex = forcedCellIndex;
    forcedFaceIndex[faceNormal] += faceIndexOffset;

    // Advective term
    ibSource += mesh.cellLengthsInv[faceNormal](forcedCellIndex[faceNormal]) 
                * faceFluxes[faceNormal](forcedFaceIndex) 
                * faceData.fieldValues.U[momentumAxis]

                - mesh.cellLengthsInv[faceNormal](forcedCellIndex[faceNormal]) 
                * std::min( faceFluxes[faceNormal](forcedFaceIndex), static_cast<floatType>(0.0f) ) 
                * fields.U[momentumAxis](forcedCellIndex);

    // Diffusive term
    ibSource += inputData.nu * (
                2 * mesh.cellLengthsInv[faceNormal](forcedCellIndex[faceNormal]) 
                    * ( faceData.fieldValues.U[momentumAxis] - fields.U[momentumAxis](forcedCellIndex) )

                    - mesh.cellCenterDiffInv[faceNormal](forcedFaceIndex[faceNormal]) * fields.U[momentumAxis](forcedCellIndex)
                );

    // Pressure gradient term
    ibSource += ( 
                mesh.cellLengthsInv[faceNormal](forcedCellIndex[faceNormal]) * faceData.fieldValues.P
    
                - mesh.interpFactors[faceNormal](forcedFaceIndex[faceNormal]) * fields.P(forcedCellIndex)
                ) / inputData.rho;

    return ibSource;
}



floatType HiSideMomentumIBSource( const Axis::ENUMDATA momentumAxis,
                                  const IBCell::FaceData &faceData, 
                                  const TensorIndex3D &forcedCellIndex,
                                  const EnumVector<Axis, Tensor3D> &faceFluxes,
                                  const FieldData<Tensor3D> &fields,
                                  const Mesh &mesh,
                                  const InputData &inputData ) 
{
    floatType ibSource = 0.0f;

    intType faceIndexOffset   = faceData.faceIndexOffset;
    Axis::ENUMDATA faceNormal = faceData.faceNormal;
    TensorIndex3D forcedFaceIndex = forcedCellIndex;
    forcedFaceIndex[faceNormal] += faceIndexOffset;

    // Advective term
    ibSource += - mesh.cellLengthsInv[faceNormal](forcedCellIndex[faceNormal]) 
                    * faceFluxes[faceNormal](forcedFaceIndex) 
                    * faceData.fieldValues.U[momentumAxis]

                + mesh.cellLengthsInv[faceNormal](forcedCellIndex[faceNormal]) 
                    * std::max( faceFluxes[faceNormal](forcedFaceIndex), static_cast<floatType>(0.0f) ) 
                    * fields.U[momentumAxis](forcedCellIndex);

    // Diffusive term
    ibSource += inputData.nu * (
                2 * mesh.cellLengthsInv[faceNormal](forcedCellIndex[faceNormal]) 
                    * ( faceData.fieldValues.U[momentumAxis] - fields.U[momentumAxis](forcedCellIndex) )

                    - mesh.cellCenterDiffInv[faceNormal](forcedFaceIndex[faceNormal]) * fields.U[momentumAxis](forcedCellIndex)
                );

    // Pressure gradient term
    ibSource += ( 
                - mesh.cellLengthsInv[faceNormal](forcedCellIndex[faceNormal]) * faceData.fieldValues.P
    
                + ( 1 - mesh.interpFactors[faceNormal](forcedFaceIndex[faceNormal]) ) * fields.P(forcedCellIndex)
                ) / inputData.rho;

    return ibSource;
}



floatType LoSideContinuityIBSource( const IBCell::FaceData &faceData, 
                                    const TensorIndex3D &forcedCellIndex,
                                    const EnumVector<Axis, Tensor3D> &faceFluxes,
                                    const FieldData<Tensor3D> &fields,
                                    const Mesh &mesh,
                                    const InputData &inputData ) 
{

    floatType ibSource = 0.0f;


    intType faceIndexOffset   = faceData.faceIndexOffset;
    Axis::ENUMDATA faceNormal = faceData.faceNormal;
    TensorIndex3D forcedFaceIndex = forcedCellIndex;
    forcedFaceIndex[faceNormal] += faceIndexOffset;

    // Divergence term
    ibSource += mesh.cellLengthsInv[faceNormal](forcedCellIndex[faceNormal]) * faceData.fieldValues.U[faceNormal]
    
                - mesh.interpFactors[faceNormal](forcedFaceIndex[faceNormal]) * fields.U[faceNormal](forcedCellIndex);

    // MWI term
    ibSource += 0.0f;

    return ibSource;

}



floatType HiSideContinuityIBSource( const IBCell::FaceData &faceData, 
                                    const TensorIndex3D &forcedCellIndex,
                                    const EnumVector<Axis, Tensor3D> &faceFluxes,
                                    const FieldData<Tensor3D> &fields,
                                    const Mesh &mesh,
                                    const InputData &inputData ) 
{

    floatType ibSource = 0.0f;


    intType faceIndexOffset   = faceData.faceIndexOffset;
    Axis::ENUMDATA faceNormal = faceData.faceNormal;
    TensorIndex3D forcedFaceIndex = forcedCellIndex;
    forcedFaceIndex[faceNormal] += faceIndexOffset;

    // Divergence term
    ibSource += - mesh.cellLengthsInv[faceNormal](forcedCellIndex[faceNormal]) * faceData.fieldValues.U[faceNormal]
    
                + ( 1 - mesh.interpFactors[faceNormal](forcedFaceIndex[faceNormal]) ) * fields.U[faceNormal](forcedCellIndex);

    // MWI term
    ibSource += 0.0f;

    return ibSource;

}



}   // end anonymous namespace





void AddIBSourceTerms( FVCoefficients &fvCoeffs,
                       const IBData &ibData, 
                       const EnumVector<Axis, Tensor3D> &faceFluxes,
                       const FieldData<Tensor3D> &fields,
                       const Mesh &mesh,
                       const InputData &inputData )
{

    // Iterate through each forced cell
    for ( auto &ibCell : ibData.IBCells ) { 

        TensorIndex3D forcedCellIndex = ibCell.cellIndex;

        // A source term is added for each forced face
        for ( auto &faceData : ibCell.facesData ) {

            intType faceIndexOffset   = faceData.faceIndexOffset;
            Axis::ENUMDATA faceNormal = faceData.faceNormal;
            TensorIndex3D forcedFaceIndex = forcedCellIndex;
            forcedFaceIndex[faceNormal] += faceIndexOffset;


            if ( faceData.faceIndexOffset == 0 ) {          // Forced face on low side

                // Momentum equations
                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                    fvCoeffs.Mom[axis].B( forcedCellIndex ) += LoSideMomentumIBSource( axis, faceData, forcedCellIndex, faceFluxes, fields, mesh, inputData );
                } );

                // Continuity equation
                fvCoeffs.Cont.B( forcedCellIndex ) += LoSideContinuityIBSource( faceData, forcedCellIndex, faceFluxes, fields, mesh, inputData );

            } else if ( faceData.faceIndexOffset == 1 ) {   // Forced face on high side

                 // Momentum equations
                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                    fvCoeffs.Mom[axis].B( forcedCellIndex ) += HiSideMomentumIBSource( axis, faceData, forcedCellIndex, faceFluxes, fields, mesh, inputData );
                } );

                // Continuity equation
                fvCoeffs.Cont.B( forcedCellIndex ) += HiSideContinuityIBSource( faceData, forcedCellIndex, faceFluxes, fields, mesh, inputData );

            }

        }

    }

}



void SetIBFaceVelocities( EnumVector<Axis, Tensor3D> &faceFluxes,
                          const IBData &ibData ) 
{
    for ( auto &ibCell : ibData.IBCells ) { 
        for ( auto &faceData : ibCell.facesData ) {

            TensorIndex3D faceIndex = ibCell.cellIndex;    
            faceIndex[faceData.faceNormal] += faceData.faceIndexOffset;

            faceFluxes[faceData.faceNormal](faceIndex) = faceData.fieldValues.U[faceData.faceNormal];

        }
    }
}



void UpdateForcedFaceFieldValues( IBData &ibData, 
                                  const FieldData<Tensor3D> &fields )
{
    for ( auto &ibCell : ibData.IBCells ) { 
        for ( auto &faceData : ibCell.facesData ) {

            faceData.fieldValues = GetForcedFaceValues( ibCell.cellIndex, faceData, fields );

        }
    }
}



void MaskFields( FieldData<Tensor3D> &fields, 
                 const Tensor3D &mask )
{
    Eigen::array<Eigen::Index, 3> offsets = { nGhost, nGhost, nGhost };

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        Eigen::array<Eigen::Index, 3> extents = { fields.U[axis].dimension(0) - 2*nGhost,
                                                  fields.U[axis].dimension(1) - 2*nGhost,
                                                  fields.U[axis].dimension(2) - 2*nGhost };

        fields.U[axis].slice( offsets, extents ) *= mask;
    } );
}


}   // end namespace CFD