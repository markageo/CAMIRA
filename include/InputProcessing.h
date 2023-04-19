#ifndef INPUT_PROCESSING
#define INPUT_PROCESSING

#include "Types.h"
#include "boost/property_tree/ptree.hpp"
#include <vector>
#include <utility>
#include <map>

namespace CFD
{


struct InputData
{
    
    // Constructor
    InputData();

    // Model
    floatType rho, nu;


    // Domain size
    floatVector3 domainSize;


    // Mesh
    struct MeshSegment {
        floatType lowerBound;
        floatType upperBound;
        intType nCells;
        floatType biasFactor;
    };
    EnumVector< Axis, std::vector< MeshSegment > > meshSegments;


    // Boundary conditions
    struct BoundaryConditionStruct {
        BoundaryConditions::ENUMDATA type;
        CFD::floatType value;
    };
    using BoundaryConditionData = EnumVector< Fields, EnumVector< BoundaryPatches, BoundaryConditionStruct > >;
    BoundaryConditionData boundaryConditions;
    std::map< BoundaryPatches::ENUMDATA, BoundaryPatches::ENUMDATA > axisTransformation; // Code patch -> user patch

    
};

InputData ReadInputData(const std::string &);

InputData InputDataFromCommandLine(int, char const **);

} // end namespace CFD

#endif  // INPUT_PROCESSING