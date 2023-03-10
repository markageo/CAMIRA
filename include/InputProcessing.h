#ifndef INPUT_PROCESSING
#define INPUT_PROCESSING

#include "SimulationParameters.h"
#include "boost/property_tree/ptree.hpp"
#include <vector>
#include <utility>
#include <optional>

namespace CFD
{

struct InputData
{

    CFD::floatType domainSize_x, domainSize_y, domainSize_z;

    struct MeshSegment {
        CFD::floatType lowerBound;
        CFD::floatType upperBound;
        CFD::intType nCells;
        CFD::floatType biasFactor;
    };
    std::vector<MeshSegment> meshSegments_x, meshSegments_y, meshSegments_z;

};

std::optional<InputData> ReadInputData(const std::string &);

} // end namespace CFD

#endif  // INPUT_PROCESSING