#ifndef INPUT_PROCESSING
#define INPUT_PROCESSING

#include "Types.h"
#include "boost/property_tree/ptree.hpp"
#include <vector>
#include <utility>

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
    std::vector<MeshSegment> meshSegments_x, meshSegments_y, meshSegments_z;


    // Boundary conditions
    struct BoundaryConditionStruct {
        BoundaryConditions::ENUMDATA type;
        CFD::floatType value;
    };
    std::vector< std::vector< BoundaryConditionStruct > > boundaryConditions; 
    
};

InputData ReadInputData(const std::string &);

} // end namespace CFD

#endif  // INPUT_PROCESSING