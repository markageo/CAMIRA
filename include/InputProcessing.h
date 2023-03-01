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

    // mesh
    struct Mesh{
        std::vector<std::pair<SIM::floatType, SIM::floatType>> segmentBounds_x, segmentBounds_y, segmentBounds_z;
        std::vector<SIM::floatType> biasFactors_x, biasFactors_y, biasFactors_z;
        std::vector<SIM::intType> nCells_x, nCells_y, nCells_z;
    } mesh;

};

std::optional<InputData> ReadInputData(const std::string &);

#endif  // INPUT_PROCESSING