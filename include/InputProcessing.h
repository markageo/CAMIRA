#ifndef INPUT_PROCESSING
#define INPUT_PROCESSING

#include "SimulationParameters.h"
#include "boost/property_tree/ptree.hpp"
#include <vector>
#include <utility>
#include <optional>


struct InputData
{

    SIM::floatType domainSize_x, domainSize_y, domainSize_z;

    struct MeshSegment {
        SIM::floatType lowerBound;
        SIM::floatType upperBound;
        SIM::intType nCells;
        SIM::floatType biasFactor;
    };
    std::vector<MeshSegment> meshSegments_x, meshSegments_y, meshSegments_z;

};

std::optional<InputData> ReadInputData(const std::string &);

#endif  // INPUT_PROCESSING