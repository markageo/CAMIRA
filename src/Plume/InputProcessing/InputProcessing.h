#ifndef CAMIRA_PLUME_INPUT_PROCESSING
#define CAMIRA_PLUME_INPUT_PROCESSING

#include "Core/Types.h"
#include "Core/Mesh/Mesh.h"
#include "Core/Geometry/Geometry.h"
#include "Plume/ConfigEnums.h"
#include "boost/property_tree/ptree.hpp"
#include <vector>
#include <utility>
#include <map>
#include <tuple>
#include <memory>

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{


// -------------------------------------- Definition in InputProcessing.cpp -------------------------------------- //

// Holds all use input data from input file
struct InputData
{

    // File directory of input file
    std::string inputFileDirectory;

    // File directory of velocity field
    std::string velocityFieldFilename;

    // File directory of geometry
    std::string stlGeometryFilename;

    // Model
    floatType diffusionCoeff;
    intType numberOfTimeSteps;
    floatType initialParticlesPerUnitMass;
    floatType timeStepSize;

    
    // Sources
    struct ContinuousReleasePointData {
        fArray3 location;
        floatType massFlowRate;
    };
    std::vector<ContinuousReleasePointData> continuousReleasePoints;

    struct InstantaneousReleasePointData {
        fArray3 location;
        floatType totalMass;
    };
    std::vector<InstantaneousReleasePointData> instantaneousReleasePoints;


    // Boundary conditions
    struct BoundaryConditionInputData {
        BoundaryConditions type;
    };
    EnumVector< BoundaryPatches, BoundaryConditionInputData > boundaryConditions;

    struct ParallelSettings {
        intType numberOfThreads;
    };
    ParallelSettings parallelSettings;

    // Ouput
    enum class OutputFormatType {
        BINARY, ASCII
    };
    OutputFormatType outputFormatType;
    std::string fieldOutputFilename;
    std::string profilingFilename;
    intType fieldWriteInterval;

    struct ProbeData {
        std::string filename;
        fArray3 location;
    };
    std::vector< ProbeData > probes;
};


InputData ReadInputData(const std::string &);

InputData InputDataFromCommandLine(int, char const **);


}   // end namespace PLUME

}   // end namespace CAMIRA

#endif  // CAMIRA_PLUME_INPUT_PROCESSING