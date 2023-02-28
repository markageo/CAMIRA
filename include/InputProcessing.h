#ifndef INPUT_PROCESSING
#define INPUT_PROCESSING

#include "SimulationParameters.h"
#include "boost/property_tree/ptree.hpp"
#include <vector>
#include <utility>
#include <optional>


struct InputData
{

    SIM::floatType domainSizeX, domainSizeY, domainSizeZ;

    // mesh
    struct {
        std::vector<std::pair<SIM::floatType, SIM::floatType>> xSegmentBounds, ySegmentBounds, zSegmentBounds;
        std::vector<SIM::floatType> xBiasFactors, yBiasFactors, zBiasFactors;
        std::vector<SIM::floatType> xCells, yCells, zCells;
    } mesh;

};

std::optional<InputData> ReadInputData(const std::string &);

#endif  // INPUT_PROCESSING