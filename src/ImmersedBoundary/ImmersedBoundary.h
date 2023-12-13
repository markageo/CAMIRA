#ifndef CFD_IMMERSED_BOUNDARY
#define CFD_IMMERSED_BOUNDARY

#include "../Types.h"
#include "../Geometry/Geometry.h"
#include "../Tools/FieldProbe.h"

#include <vector>

namespace CFD {

enum class CellType {
    Fluid    = 0, 
    Solid    = 1
};
using CellIDTensor3D = Eigen::Tensor< CellType, 3 >;


// --------------------------------------- Definition in IBData.cpp -------------------------------------- //

struct IBCell {
    TensorIndex3D cellIndex;
    
    struct FaceData {
        TransportCoefficients::ENUMDATA face;
        floatType     interpCoeffCell,    // Interpolation coeff that multiplies with cell value
                      interpCoeffIB;      // Interpolation coefficient that multiplies with the IB value
        floatType     extrapFactor_p,     // Boundary cell extrapolation coefficient 
                      extrapFactor_a;     // One from boundary cell extrapolation coefficient    
        TensorIndex3D adjacentCellIndex;  // One from boundary cell index, for extrapolation 
    };
    std::vector< FaceData > facesData;
       
};

struct IBData {
    std::vector< IBCell > IBCells;
    Tensor3D mask;
};


IBData CreateImmersedBoundaryData( const InputData &, const Mesh & );



// ------------------------------- Definition in IBSolverFunctions.cpp --------------------------------- //


// Set all solid cells in the mask to zero 
void MaskFields( FieldData<Tensor3D> &, const Tensor3D & );

}   // end namespace CFD



#endif  // CFD_IMMERSED_BOUNDARY