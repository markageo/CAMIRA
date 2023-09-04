#include "ImmersedBoundary.h"

#include "../Tools/FVTools.h"

namespace CFD
{


void SetGhostCellValues( FieldData<Tensor3D> &fields, 
                         const IBData &ibData )
{
    using FVT::G;
    for ( const IBGhostCell &ibGhostCell : ibData.ghostCells ) {

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

            floatType imagePointValue = ibGhostCell.fieldProbe.GetFieldValue( fields[axis] );

            fVector3 imagePointGradientVector = ibGhostCell.fieldProbe.GetFieldGradient( fields[axis] );

            floatType imagePointNormalGradient = imagePointGradientVector.dot( ibGhostCell.normalUnitVector );

            fields[axis]( G( ibGhostCell.ghostCellIndex ) ) = ibGhostCell.extrapImageVelocityCoeff * imagePointValue
                                                            + ibGhostCell.extrapImageGradientCoeff * imagePointNormalGradient;

        } ); 
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