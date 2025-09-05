#ifndef CAMIRA_IMMERSED_BOUNDARY
#define CAMIRA_IMMERSED_BOUNDARY

#include "../Core/Types.h"
#include "../Geometry/Geometry.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"
#include "../FiniteVolume/Mesh.h"
#include "../DerivedQuantities/FieldProbe.h"

#include <vector>
#include <memory>

namespace CAMIRA {

enum CellType {
    Solid    = 0,
    Fluid    = 1 
};


// --------------------------------------- Definition in IBData.cpp -------------------------------------- //


// Specific coefficients and data for directional immersed boundary method
struct DirectionalIBData
{
    // Coefficients for reconstructing velocity field onto face
    floatType faceReconstructionCoeff_p,         // Multiplies with immediate cell value
              faceReconstructionCoeff_a,         // Multiplies with first interior cell
              faceReconstructionCoeff_ib;        // Multiplies with IB value

    // For extrapolating directly onto AB face from fluid cells (does not account for IB), i.e. not a reconstruction
    floatType faceExtrapCoeff_p,
              faceExtrapCoeff_a;

    // Coefficients for extrapolating onto the immersed boundary surface from fluid cells
    floatType ibExtrapCoeff_p,           // Boundary cell
              ibExtrapCoeff_a;           // One interior of boundary cell  

    // Distance of center to IB surface along coordinate direction
    floatType ibDistance;

};


// Specific coefficients and data for wall functions
struct WallFunctionData
{
    FieldProbe fieldProbe;
    fVector3 normalVector;

    // Distance of image point to IB surface along normal direction
    floatType imagePointDistance;

    // Distance of of face center to IB surface along normal direction
    floatType faceCenterDistance;
};


// Contains data for a single forced cell
struct IBCell {
    TensorIndex3D cellIndex;
    
    struct SourceTermData {
        Axis::ENUMDATA direction;
        intType        directionIndex,        // Cell index offset, either +1 or -1.
                       faceDirectionIndex;    // Face index offset, either 0 for lo side, or 1 for hi side
        TensorIndex3D cellIndex_g,            // Ghost cell
                      cellIndex_a;            // One from boundary cell index, for extrapolation 

        
        // Coefficients for extrapolating from face to ghost cell
        floatType ghostExtrapCoeff_p,       // Multiplies with immediate cell value
                  ghostExtrapCoeff_f;       // Multiplies with the face value


        // Coefficients for the far pressure ghost cell to allow for correct MWI at the immersed boundary
        floatType farPressureCoeff_p,        // Boundary cell
                  farPressureCoeff_a,        // One interior of boundary cell
                  farPressureCoeff_g;        // Ghost cell       

        // Store pointers to this data because only one of them can be used at a time
        std::unique_ptr<DirectionalIBData> directionalIBDataPtr;
        std::unique_ptr<WallFunctionData> wallFunctionDataPtr;

        // Cell face normal area vector component. This has a sign, positive in the positive coordinate direction.
        floatType faceAreaComponent;

        // Coefficient that multiplies with total velocity flux error
        floatType velocityFluxCorrectionCoeff;

        FieldData<floatType> ibValues,          // Field variables on immersed boundary surface at intersection point 
                             faceValues,        // Field variables on face between ghost cell and interior cell  
                             ghostCellValues;   // Field variables at ghost cell
        floatType farPressureGhostCellValue;
        floatType wallShear;                    // Wall shear value at the face, only needed in wall functions
    };
    std::vector< SourceTermData > sourceTermsData;
       
};

struct IBData {
    std::vector< std::vector< IBCell > > ibCells;   // Outer vector is for each geometry component, inner vector is the cells
    Tensor3D mask;
};


void SetImmersedBoundaryData( IBData &, const InputData &, const AxisTransformationMap &, const Mesh &);

void WriteGeometryToFile( const InputData &, const AxisTransformationMap & );


// ------------------------------- Definition in IBSolverFunctions.cpp --------------------------------- //


// Update values of ghost cells by re-interpolating from the immersed boundary
void UpdateIBData( IBData &, const FieldData<Tensor3D> & );


// Set all solid cells in the mask to zero 
void MaskFields( FieldData<Tensor3D> &, const Tensor3D & );



}   // end namespace CAMIRA



#endif  // CAMIRA_IMMERSED_BOUNDARY