#include "ImmersedBoundary.h"

namespace CFD
{


void SetGhostCellValues( FieldData<Tensor3D> &fields, 
                         const IBData &ibData )
{
    using enum Axis::ENUMDATA;
    using ipVector = Eigen::Matrix< floatType, IBGhostCell::numInterpPoints, 1 >;

    for ( const IBGhostCell &ibGhostCell : ibData.ghostCells ) {

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

            // Vector of field values
            ipVector fieldValues = { 0.0f,  // Boundary value
                                     fields.U[axis]( ibGhostCell.fluidCellIndices[0] ),
                                     fields.U[axis]( ibGhostCell.fluidCellIndices[1] ),
                                     fields.U[axis]( ibGhostCell.fluidCellIndices[2] ) };

            // Polynomial weighting coefficients
            ipVector weightingCoeffs = ibGhostCell.pointsMatrixInv * fieldValues;

            // Ghost cell value
            floatType ghostCellValue = weightingCoeffs(0)
                                     + weightingCoeffs(1) * ibGhostCell.imagePointCoordinates(0)
                                     + weightingCoeffs(2) * ibGhostCell.imagePointCoordinates(1)
                                     + weightingCoeffs(3) * ibGhostCell.imagePointCoordinates(2);

            // Set the ghost cell
            fields.U[axis]( ibGhostCell.ghostCellIndex ) = ghostCellValue;

        } ); 
    }

}



}   // end namespace CFD