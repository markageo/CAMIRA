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

    // Domain size
    CFD::floatType domainSize_x, domainSize_y, domainSize_z;


    // Mesh
    struct MeshSegment {
        CFD::floatType lowerBound;
        CFD::floatType upperBound;
        CFD::intType nCells;
        CFD::floatType biasFactor;
    };
    std::vector< std::vector<MeshSegment> > meshSegments;   // meshSegments[Axis][segment]



    // Boundary conditions
    struct BoundaryConditionStruct {
        BoundaryConditions::ENUMDATA type;
        CFD::floatType value;
    };
    using BoundaryConditionData = std::vector< std::vector< BoundaryConditionStruct > >;
    BoundaryConditionData boundaryConditions; // boundaryConditions[Field][Patch]
    std::map< CFD::BoundaryPatches::ENUMDATA, CFD::BoundaryPatches::ENUMDATA > axisTransformation; // Code patch -> user patch

    
};

InputData ReadInputData(const std::string &);

} // end namespace CFD

#endif  // INPUT_PROCESSING