#ifndef CFD_IMMERSED_BOUNDARY
#define CFD_IMMERSED_BOUNDARY

#include "../Types.h"
#include "../Geometry/Geometry.h"
#include "../Tools/FieldProbe.h"
#include "../FiniteVolume/FiniteVolume.h"

#include <vector>

namespace CFD {

enum CellType {
    Solid    = 0,
    Fluid    = 1 
};


// --------------------------------------- Definition in IBData.cpp -------------------------------------- //


// For new method
struct IBCell {
    TensorIndex3D cellIndex;
    
    struct SourceTermData {
        Axis::ENUMDATA direction;
        intType        directionIndex,        // Cell index offset, either +1 or -1.
                       faceDirectionIndex;    // Face index offset, either 0 for lo side, or 1 for hi side
        TensorIndex3D cellIndex_a;            // One from boundary cell index, for extrapolation 

        // Coefficients for extrapolating onto the ghost cell
        floatType ghostExtrapCoeff_p,         // Multiplies with immediate cell value
                  ghostExtrapCoeff_a,         // Multiplies with first interior cell
                  ghostExtrapCoeff_ib;        // Multiplies with IB value

        // Coefficients for interpolating onto the face between the ghost cell
        floatType faceInterpCoeff_p,          // Multiplies with cell value
                  faceInterpCoeff_g;          // Multiples with ghost cell value

        // Coefficients for extrapolating onto the immersed boundary surface
        floatType ibExtrapFactor_p,           // Boundary cell
                  ibExtrapFactor_a;           // One interior of boundary cell  
        
        // Coefficients for the far pressure ghost cell to allow for correct MWI at the immersed boundary
        floatType farPressureCoeff_p,        // Boundary cell
                  farPressureCoeff_a,        // One interior of boundary cell
                  farPressureCoeff_g;        // Ghost cell       

        FieldData<floatType> ibFieldValues; // Value of fields on the immersed boundary surface at the intersection point

        FieldData<floatType> ghostCellValues;
        floatType farPressureGhostCellValue;
    };
    std::vector< SourceTermData > sourceTermsData;
       
};

struct IBData {
    std::vector< IBCell > ibCells;
    Tensor3D mask;
};


IBData CreateImmersedBoundaryData( const InputData &, const Mesh & );



// ------------------------------- Definition in IBSolverFunctions.cpp --------------------------------- //


// Add source terms to finite volume equations which include the effect of the immersed boundary
void AddIBSourceTerms( FVCoefficients &, const IBData & );


// Set the face velocities to their interpolated values according to the immersed boundary
void SetIBFaceFluxes( EnumVector<Axis, Tensor3D> &, const IBData &, const FieldData<Tensor3D> & );


// Update values of ghost cells by re-interpolating from the immersed boundary
void UpdateIBData( IBData &, const FieldData<Tensor3D> & );


// Set all solid cells in the mask to zero 
void MaskFields( FieldData<Tensor3D> &, const Tensor3D & );



}   // end namespace CFD



#endif  // CFD_IMMERSED_BOUNDARY