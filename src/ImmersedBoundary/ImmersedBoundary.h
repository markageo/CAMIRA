#ifndef CFD_IMMERSED_BOUNDARY
#define CFD_IMMERSED_BOUNDARY

#include "../Types.h"
#include "../Geometry/Geometry.h"

namespace CFD {

enum class CellType {
    Fluid, Ghost, Solid
};
using CellIDTensor3D = Eigen::Tensor< CellType, 3 >;


// ----------------------------------------- Definition in IBData.cpp ---------------------------------------- //

struct IBGhostCell {
    static constexpr intType numInterpolationPoints = 4;
    fVector3 boundaryPointCoordinates;
    fVector3 imagePointCoordiantes;
    std::array< TensorIndex3D, numInterpolationPoints > fluidCellIndices;
    TensorIndex3D ghostCellIndex;
    Eigen::Matrix< floatType, numInterpolationPoints, numInterpolationPoints > interpolationMatrixInv;
};

struct IBData {
    std::vector< IBGhostCell > ghostCells;
    Tensor3D mask;
};


IBData CreateImmersedBoundaryData( const CellIDTensor3D &, const Mesh & );




// -------------------------------------- Definition in IBGeometry.cpp -------------------------------------- //

// Create and ID array to identify different cells in immersed boundary
CellIDTensor3D TagCells( const Mesh &, const Polyhedron & );


}   // end namespace CFD



#endif  // CFD_IMMERSED_BOUNDARY