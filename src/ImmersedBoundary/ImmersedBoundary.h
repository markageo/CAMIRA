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

// Contains data for a single forced cell
struct IBCell {
    TensorIndex3D cellIndex;
    
    struct SourceTermData {
        Axis::ENUMDATA direction;
        intType        directionIndex,        // Cell index offset, either +1 or -1.
                       faceDirectionIndex;    // Face index offset, either 0 for lo side, or 1 for hi side
        TensorIndex3D cellIndex_a;            // One from boundary cell index, for extrapolation 

        // Coefficients for extrapolating onto face
        floatType faceExtrapCoeff_p,         // Multiplies with immediate cell value
                  faceExtrapCoeff_a,         // Multiplies with first interior cell
                  faceExtrapCoeff_ib;        // Multiplies with IB value

        // Coefficients for extrapolating from face to ghost cell
        floatType ghostExtrapCoeff_p,       // Multiplies with immediate cell vale
                  ghostExtrapCoeff_f;       // Multiplies with the face value


        // Coefficients for extrapolating onto the immersed boundary surface
        floatType ibExtrapFactor_p,           // Boundary cell
                  ibExtrapFactor_a;           // One interior of boundary cell  
        
        // Coefficients for the far pressure ghost cell to allow for correct MWI at the immersed boundary
        floatType farPressureCoeff_p,        // Boundary cell
                  farPressureCoeff_a,        // One interior of boundary cell
                  farPressureCoeff_g;        // Ghost cell       

        // Coordinate distance squared to the immersed boundary surface
        floatType ibDistance2;

        // Cell face normal area vector component. This has a sign.
        floatType faceAreaComponent;

        // Coefficient that multiplies with total velocity flux error
        floatType velocityFluxCorrectionCoeff;

        FieldData<floatType> ibValues,          // Field variables on immersed boundary surface at intersection point 
                             faceValues,        // Field variables on face between ghost cell and interior cell  
                             ghostCellValues;   // Field variables at ghost cell
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
void SetIBFaceFluxes( EnumVector<Axis, Tensor3D> &, const IBData & );


// Update values of ghost cells by re-interpolating from the immersed boundary
void UpdateIBData( IBData &, const FieldData<Tensor3D> & );


// Set all solid cells in the mask to zero 
void MaskFields( FieldData<Tensor3D> &, const Tensor3D & );



}   // end namespace CFD



#endif  // CFD_IMMERSED_BOUNDARY