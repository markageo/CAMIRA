#include "ImmersedBoundary.h"


namespace CFD
{


namespace 
{



}   // end anonymous namespace





IBData CreateImmersedBoundaryData( const CellIDTensor3D &cellID, 
                                   const Mesh &mesh )
{

    using enum Axis::ENUMDATA;

    IBData ibData;
    ibData.mask = Tensor3D( mesh.nCells[X], mesh.nCells[Y], mesh.nCells[Z] ).setZero();

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                if ( cellID( i, j, k ) != CellType::Ghost ) {

                    if ( cellID(i, j, k) == CellType::Fluid ) {
                        ibData.mask(i, j, k) = 1.0f;
                    }
                    continue;
                }

                // Find coordinates of nearest immersed boundary point


                // Get the normal vector

                // Determine coordinates of image point

                // Find nearest fluid cells for interpolation

                // Form the interpolation coefficient matrix


            }
        }
    }

    return ibData;

}


}