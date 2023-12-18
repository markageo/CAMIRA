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
                                  const FVCoefficients &fvCoeffs,
                                  const FieldData<Tensor3D> &fields,
                                  const Mesh &mesh ) 
{
    Axis::ENUMDATA faceNormal = faceData.faceNormal;
    intType fidx = forcedCellIndex[faceNormal] + faceData.faceIndexOffset;
    TransportCoefficients::ENUMDATA loCoeff = LUT::LoCoeff[faceNormal];
    floatType interpFactor = mesh.interpFactors[faceNormal](fidx);

    // Velocity terms
    floatType ibSource = - fvCoeffs.Mom[momentumAxis].AU[momentumAxis][loCoeff](forcedCellIndex) 
                         * ( 1.0f / interpFactor )
                         * faceData.fieldValues.U[momentumAxis]
               
                       +    fvCoeffs.Mom[momentumAxis].AU[momentumAxis][loCoeff](forcedCellIndex) 
                         * ( (1.0f - interpFactor) / interpFactor )
                         * fields.U[momentumAxis](forcedCellIndex);

    // Pressure terms
    if ( momentumAxis == faceNormal ) {

        ibSource += - fvCoeffs.Mom[momentumAxis].AP[loCoeff](forcedCellIndex[faceNormal]) 
                    * ( 1.0f / interpFactor )
                    * faceData.fieldValues.P
               
                  +    fvCoeffs.Mom[momentumAxis].AP[loCoeff](forcedCellIndex[faceNormal]) 
                    * ( (1.0f - interpFactor) / interpFactor )
                    * fields.P(forcedCellIndex);

    }

    return ibSource;
}



floatType HiSideMomentumIBSource( const Axis::ENUMDATA momentumAxis,
                                  const IBCell::FaceData &faceData, 
                                  const TensorIndex3D &forcedCellIndex,
                                  const FVCoefficients &fvCoeffs,
                                  const FieldData<Tensor3D> &fields,
                                  const Mesh &mesh ) 
{
    Axis::ENUMDATA faceNormal = faceData.faceNormal;
    intType fidx = forcedCellIndex[faceNormal] + faceData.faceIndexOffset;
    TransportCoefficients::ENUMDATA loCoeff = LUT::HiCoeff[faceNormal];
    floatType interpFactor = mesh.interpFactors[faceNormal](fidx);

    // Velocity terms
    floatType ibSource = - fvCoeffs.Mom[momentumAxis].AU[momentumAxis][loCoeff](forcedCellIndex) 
                         * ( 1.0f / ( 1.0f - interpFactor ) )
                         * faceData.fieldValues.U[momentumAxis]
               
                       +    fvCoeffs.Mom[momentumAxis].AU[momentumAxis][loCoeff](forcedCellIndex) 
                         * ( interpFactor / (1.0f - interpFactor) )
                         * fields.U[momentumAxis](forcedCellIndex);

    // Pressure terms
    if ( momentumAxis == faceNormal ) {

        ibSource += - fvCoeffs.Mom[momentumAxis].AP[loCoeff](forcedCellIndex[faceNormal]) 
                    * ( 1.0f / ( 1.0f - interpFactor ) )
                    * faceData.fieldValues.P
               
                  +    fvCoeffs.Mom[momentumAxis].AP[loCoeff](forcedCellIndex[faceNormal]) 
                    * ( interpFactor / (1.0f - interpFactor) )
                    * fields.P(forcedCellIndex);

    }

    return ibSource;
}



floatType LoSideContinuityIBSource( const IBCell::FaceData &faceData, 
                                    const TensorIndex3D &forcedCellIndex,
                                    const FVCoefficients &fvCoeffs,
                                    const FieldData<Tensor3D> &fields,
                                    const Mesh &mesh ) 
{
    Axis::ENUMDATA faceNormal = faceData.faceNormal;
    intType fidx = forcedCellIndex[faceNormal] + faceData.faceIndexOffset;
    TransportCoefficients::ENUMDATA loCoeff = LUT::LoCoeff[faceNormal];
    floatType interpFactor = mesh.interpFactors[faceNormal](fidx);

    floatType ibSource = - fvCoeffs.Cont.AU[faceNormal][loCoeff](forcedCellIndex[faceNormal]) 
                         * ( 1.0f / interpFactor )
                         * faceData.fieldValues.U[faceNormal]
               
                       +    fvCoeffs.Cont.AU[faceNormal][loCoeff](forcedCellIndex[faceNormal]) 
                         * ( (1.0f - interpFactor) / interpFactor )
                         * fields.U[faceNormal](forcedCellIndex);

    // Pressure terms
    intType ffidx = fidx - 1;
    TransportCoefficients::ENUMDATA loloCoeff = LUT::LoLoCoeff[faceNormal];
    floatType lw  = 1.0f / mesh.cellCenterDiffInv[faceNormal](fidx),
              lww = 1.0f / mesh.cellCenterDiffInv[faceNormal](ffidx);
    floatType wideExtrapCoeffCell = 1 - ( (lw + lww) / lw ) * ( 1 + interpFactor / (1-interpFactor) ),
              wideExtrapCoeffFace =     ( (lw + lww) / lw ) * ( 1 / (1 - interpFactor) );

    ibSource += - ( fvCoeffs.Cont.AP[loCoeff](forcedCellIndex) * ( 1.0f / (1 - interpFactor) )
                  + fvCoeffs.Cont.AP[loloCoeff](forcedCellIndex) * wideExtrapCoeffFace 
                  ) * faceData.fieldValues.P

                + ( fvCoeffs.Cont.AP[loCoeff](forcedCellIndex) * ( interpFactor / (1 - interpFactor) )
                  - fvCoeffs.Cont.AP[loloCoeff](forcedCellIndex) * wideExtrapCoeffCell 
                  ) * fields.P(forcedCellIndex);

    return ibSource;
}



floatType HiSideContinuityIBSource( const IBCell::FaceData &faceData, 
                                    const TensorIndex3D &forcedCellIndex,
                                    const FVCoefficients &fvCoeffs,
                                    const FieldData<Tensor3D> &fields,
                                    const Mesh &mesh ) 
{
    Axis::ENUMDATA faceNormal = faceData.faceNormal;
    intType fidx = forcedCellIndex[faceNormal] + faceData.faceIndexOffset;
    TransportCoefficients::ENUMDATA hiCoeff = LUT::HiCoeff[faceNormal];
    floatType interpFactor = mesh.interpFactors[faceNormal](fidx);

    // Velocity terms
    floatType ibSource = - fvCoeffs.Cont.AU[faceNormal][hiCoeff](forcedCellIndex[faceNormal]) 
                         * ( 1.0f / ( 1.0f - interpFactor ) )
                         * faceData.fieldValues.U[faceNormal]
               
                       +    fvCoeffs.Cont.AU[faceNormal][hiCoeff](forcedCellIndex[faceNormal]) 
                         * ( interpFactor / (1.0f - interpFactor) )
                         * fields.U[faceNormal](forcedCellIndex);

    // Pressure terms
    intType ffidx = fidx + 1;
    TransportCoefficients::ENUMDATA hihiCoeff = LUT::HiHiCoeff[faceNormal];
    floatType le  = 1.0f / mesh.cellCenterDiffInv[faceNormal](fidx),
              lee = 1.0f / mesh.cellCenterDiffInv[faceNormal](ffidx);
    floatType wideExtrapCoeffCell = 1 - ( (le + lee) / le ) * ( 1 + (1-interpFactor) / interpFactor ),
              wideExtrapCoeffFace =     ( (le + lee) / le ) * ( 1 / interpFactor );

    ibSource += - ( fvCoeffs.Cont.AP[hiCoeff](forcedCellIndex) * ( 1.0f / interpFactor )
                  + fvCoeffs.Cont.AP[hihiCoeff](forcedCellIndex) * wideExtrapCoeffFace 
                  ) * faceData.fieldValues.P

                + ( fvCoeffs.Cont.AP[hiCoeff](forcedCellIndex) * ( (1-interpFactor) / interpFactor )
                  - fvCoeffs.Cont.AP[hihiCoeff](forcedCellIndex) * wideExtrapCoeffCell 
                  ) * fields.P(forcedCellIndex);

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

            if ( faceData.faceIndexOffset == 0 ) {          // Forced face on low side

                // Momentum equations
                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                    fvCoeffs.Mom[axis].B( forcedCellIndex ) += LoSideMomentumIBSource( axis, faceData, forcedCellIndex, fvCoeffs, fields, mesh );
                } );

                // Continuity equation
                fvCoeffs.Cont.B( forcedCellIndex ) += LoSideContinuityIBSource( faceData, forcedCellIndex, fvCoeffs, fields, mesh);

            } else if ( faceData.faceIndexOffset == 1 ) {   // Forced face on high side

                 // Momentum equations
                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                    fvCoeffs.Mom[axis].B( forcedCellIndex ) += HiSideMomentumIBSource( axis, faceData, forcedCellIndex, fvCoeffs, fields, mesh );
                } );

                // Continuity equation
                fvCoeffs.Cont.B( forcedCellIndex ) += HiSideContinuityIBSource( faceData, forcedCellIndex, fvCoeffs, fields, mesh );

            }

        }

    }

}



void SetIBFaceFluxes( EnumVector<Axis, Tensor3D> &faceFluxes,
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

    ForAllFieldData( [&] (intType f) {

        Eigen::array<Eigen::Index, 3> extents = { fields[f].dimension(0) - 2*nGhost,
                                                  fields[f].dimension(1) - 2*nGhost,
                                                  fields[f].dimension(2) - 2*nGhost };

        fields[f].slice( offsets, extents ) *= mask;

    } );
    
}


}   // end namespace CFD