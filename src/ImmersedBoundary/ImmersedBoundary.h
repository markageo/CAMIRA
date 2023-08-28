#ifndef CFD_IMMERSED_BOUNDARY
#define CFD_IMMERSED_BOUNDARY

#include "../Types.h"
#include "../Geometry/Geometry.h"

namespace CFD {

enum class CellType {
    Fluid = 0, 
    Ghost = 1, 
    Solid = 2
};
using CellIDTensor3D = Eigen::Tensor< CellType, 3 >;


// --------------------------------------- Definition in IBData.cpp -------------------------------------- //

struct IBGhostCell {
    static constexpr intType numInterpPoints = 4;
    static constexpr intType numFluidInterpPoints = numInterpPoints - 1;
    using InterpMatrix = Eigen::Matrix< floatType, numInterpPoints, numInterpPoints >;

    std::array< TensorIndex3D, numFluidInterpPoints > fluidCellIndices;
    TensorIndex3D ghostCellIndex;
    fVector3 imagePointCoordinates;
    InterpMatrix pointsMatrixInv;
};

struct IBData {
    std::vector< IBGhostCell > ghostCells;
    Tensor3D mask;
};


IBData CreateImmersedBoundaryData( const InputData &, const Mesh & );





// --------------------------------- Definition in IBGeometry.cpp -------------------------------------- //

// Create and ID array to identify different cells in immersed boundary
CellIDTensor3D TagCells( const Mesh &, const Polyhedron & );





// ------------------------------- Definition in IBSolverFunctions.cpp --------------------------------- //

// Set ghost cell values in velocity fields
void SetGhostCellValues( FieldData<Tensor3D> &, const IBData & );




}   // end namespace CFD



#endif  // CFD_IMMERSED_BOUNDARY