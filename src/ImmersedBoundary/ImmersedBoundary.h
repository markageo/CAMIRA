#ifndef CFD_IMMERSED_BOUNDARY
#define CFD_IMMERSED_BOUNDARY

#include "../Types.h"
#include "../Geometry/Geometry.h"
#include "../Tools/FieldProbe.h"
#include "../FiniteVolume/FiniteVolume.h"

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
        Axis::ENUMDATA faceNormal;
        intType faceIndexOffset;          // Converts cell index to face index. 0: face on low side, 1: face on high side
        floatType     interpCoeffCell,    // Interpolation coeff that multiplies with cell value
                      interpCoeffIB;      // Interpolation coefficient that multiplies with the IB value
        floatType     extrapFactor_p,     // Boundary cell extrapolation coefficient 
                      extrapFactor_a;     // One from boundary cell extrapolation coefficient    
        TensorIndex3D adjacentCellIndex;  // One from boundary cell index, for extrapolation 
        FieldData<floatType> fieldValues;
    };
    std::vector< FaceData > facesData;
       
};

struct IBData {
    std::vector< IBCell > IBCells;
    Tensor3D mask;
};


IBData CreateImmersedBoundaryData( const InputData &, const Mesh & );



// ------------------------------- Definition in IBSolverFunctions.cpp --------------------------------- //

// Add source terms to finite volume equations which include the effect of the immersed boundary
void AddIBSourceTerms( FVCoefficients &, const IBData &, const FieldData<Tensor3D> &, const Mesh & );


// Set the face velocities to their interpolated values according to the immersed boundary
void SetIBFaceFluxes( EnumVector<Axis, Tensor3D> &, const IBData & );


// Re-extrapolate and interpolate new field values onto the forced faces and store them in IBData
void UpdateForcedFaceFieldValues( IBData &, const FieldData<Tensor3D> & );


// Set all solid cells in the mask to zero 
void MaskFields( FieldData<Tensor3D> &, const Tensor3D & );

}   // end namespace CFD



#endif  // CFD_IMMERSED_BOUNDARY