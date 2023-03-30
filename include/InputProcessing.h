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
    std::vector< std::vector<MeshSegment> > meshSegments;   // meshSegments[Axis][segment]



    // Boundary conditions
    struct BoundaryConditionStruct {
        BoundaryConditions::ENUMDATA type;
        CFD::floatType value;
    };
    using BoundaryConditionData = std::vector< std::vector< BoundaryConditionStruct > >;
    BoundaryConditionData boundaryConditions; // boundaryConditions[Field][Patch]
    std::map< BoundaryPatches::ENUMDATA, BoundaryPatches::ENUMDATA > axisTransformation; // Code patch -> user patch

    
};

InputData ReadInputData(const std::string &);

} // end namespace CFD

#endif  // INPUT_PROCESSING