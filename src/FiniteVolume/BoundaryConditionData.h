#ifndef CAMIRA_BOUNDARY_CONDITION_DATA   
#define CAMIRA_BOUNDARY_CONDITION_DATA

#include "../Core/Types.h"

#include <memory>

namespace CAMIRA
{

// -------------------------------------- Definition in FiniteVolumeStructures.cpp -------------------------------------- //

// Structure to store all domain boundary condition information
struct BoundaryConditionData {
    struct Patch {
        BoundaryConditions::ENUMDATA type;
        Tensor2D value;
    };
    using Patches = EnumVector< BoundaryPatches, Patch >;
    FieldData< Patches > fields; 
    bool pressureFieldIsFloating;
};



} // end namespace CAMIRA

#endif // CAMIRA_BOUNDARY_CONDITION_DATA