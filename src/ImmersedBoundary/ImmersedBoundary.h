#ifndef CFD_IMMERSED_BOUNDARY
#define CFD_IMMERSED_BOUNDARY

#include "../Types.h"
#include "../Geometry/Geometry.h"
#include "../Tools/FieldProbe.h"

namespace CFD {

enum class CellType {
    Fluid = 0, 
    Ghost = 1, 
    Solid = 2
};
using CellIDTensor3D = Eigen::Tensor< CellType, 3 >;


// --------------------------------------- Definition in IBData.cpp -------------------------------------- //

struct IBGhostCell {
    FieldProbe fieldProbe;
    TensorIndex3D ghostCellIndex;
    floatType extrapImageVelocityCoeff, extrapImageGradientCoeff;
    fVector3 normalUnitVector;
};

struct IBData {
    std::vector< IBGhostCell > ghostCells;
    Tensor3D mask;
    CellIDTensor3D cellID;
};


IBData CreateImmersedBoundaryData( const InputData &, const Mesh & );





// --------------------------------- Definition in IBGeometry.cpp -------------------------------------- //

// Create and ID array to identify different cells in immersed boundary
CellIDTensor3D TagCells( const Mesh &, const Polyhedron & );





// ------------------------------- Definition in IBSolverFunctions.cpp --------------------------------- //

// Set ghost cell values in velocity fields
void SetGhostCellValues( FieldData<Tensor3D> &, const IBData & );


// Set all solid cells in the mask to zero 
void MaskFields( FieldData<Tensor3D> &, const Tensor3D & );

}   // end namespace CFD



#endif  // CFD_IMMERSED_BOUNDARY