#ifndef INPUT_PROCESSING
#define INPUT_PROCESSING

#include "Core/Types.h"
#include "Core/Mesh/Mesh.h"
#include "Core/Geometry/Geometry.h"
#include "boost/property_tree/ptree.hpp"
#include <vector>
#include <utility>
#include <map>
#include <tuple>
#include <memory>

namespace CAMIRA
{


// -------------------------------------- Definition in InputProcessing.cpp -------------------------------------- //

// Holds all use input data from input file
struct InputData
{
    // Constructor
    InputData();

    // File directory of input file
    std::string inputFileDirectory;

    // Model
    floatType rho, nu;
    bool transient;

    // Domain bounds
    fArray3 domainLowerBounds,
            domainUpperBounds;

    TurbulenceModels turbulenceModel;
    floatType eddyViscosityRelaxation;
    struct ZEQ1ModelData {
        floatType reynoldsNumberBuildingHeight,                     // Inflow reynolds number at average building height
                  inflowTurbulenceIntensityBuildingHeight;          // Inflow turbulence intensity at average building height
    } zeq1ModelData;

    struct ZEQ2ModelData {
        Axis::ENUMDATA heightAxis;                                  // The axis which corresponds to height from ground level
        floatType averageBuildingHeight,                            // Average building height
                  inflowVelocityBuildingHeight,                     // Inflow velocity magnitude at average building height
                  inflowTKEBuildingHeight,                          // Inflow turbulence kinetic energy at average building height
                  inflowIntergralTimeScaleBuildingHeight,           // Inflow integral timescale at average building height
                  roughnessLength;   
    } zeq2ModelData;

    struct ZEQ3ModelData {
        Axis::ENUMDATA heightAxis;                                  // The axis which corresponds to height from ground level
        floatType averageBuildingHeight,                            // Average building height
                  inflowVelocityBuildingHeight,                     // Inflow velocity magnitude at average building height
                  roughnessLength;   
    } zeq3ModelData;

    struct ZEQ4ModelData {
        Axis::ENUMDATA heightAxis;                                  // The axis which corresponds to height from ground level
        floatType averageBuildingHeight,                            // Average building height
                  averageBuildingWidth,                             // Average building width
                  referenceHeight;                                  // Reference height for city
    } zeq4ModelData;

    // Mesh
    EnumVector< Axis, std::vector< MeshSegment > > meshSegments;

    // Solid Geometry
    GeometryBoundaryTreatment geoemtryBoundaryTreatment;
    GeometryData geometryData;
    bool hasIBGeometry;


    // Boundary conditions
    struct Profile1D {
        Axis::ENUMDATA axis;
        Tensor1D coordinates, values;
    };

    struct BoundaryConditionInputData {
        BoundaryConditions::ENUMDATA type;
        floatType uniformValue;    
        Profile1D profile1D;
        bool hasUniformValue = false,
             hasProfile1D    = false;
    };
    using FieldBoundaryConditions = EnumVector< BoundaryPatches, BoundaryConditionInputData >;
    FieldData< FieldBoundaryConditions > boundaryConditions;

    // Initial conditions
    enum class InitialConditionTypes { 
        uniform, vtkFile
    };
    InitialConditionTypes initialConditionType;
    std::string initialConditionsFieldFilename;
    FieldData<floatType> constantInitialConditions;

    // Solver
    struct Schemes {
        TimeSchemes timeScheme;
        floatType timeStep;
        intType numberOfTimesteps;
        MomentumInterpolation momentumInterpolation;
        AdvectionSchemes advectionScheme;
        floatType advectionBlendingFactor;

        intType maxOuterIterations;
        FieldData<floatType> maxOuterResiduals;
    };
    Schemes schemes;

    struct SmootherSettings {
        Smoothers type;
        intType maxIterations;
        FieldData<floatType> maxResiduals;
        FieldData<floatType> relaxation;

        BoundaryPatches::ENUMDATA lineSweepDirection;
        BoundaryPatches::ENUMDATA planeSweepDirection;
    }; 
    SmootherSettings smootherSettings;

    struct MultigridSettings {
        MultigridCycleType cycle;
        size_t maxCoarseLevels; // Make size_t because this is used with std::vector
        intType preSmoothingIterations,
                postSmoothingIterations,
                maxCoarseGridIterations,
                fineGridIterations;
        FieldData<floatType> maxCoarseGridResiduals;
    };
    MultigridSettings multigridSettings;

    struct ParallelSettings {
        intType numberOfThreads;
    };
    ParallelSettings parallelSettings;

    // Ouput
    enum class OutputFormatType {
        BINARY, ASCII
    };
    OutputFormatType outputFormatType;
    std::string residualHistoryFilename;
    std::string fieldOutputFilename;
    std::string profilingFilename;
    std::string geometryOutputFilename;
    bool outputGeometry;
    intType fieldWriteInterval;

    struct ProbeData {
        std::string filename;
        fArray3 location;
    };
    std::vector< ProbeData > probes;

    bool calculateForces;
    std::string forceCalculatorFilename;

    bool calculateYPlus;
    std::string yPlusCalculatorFilename;

};


InputData ReadInputData(const std::string &);

InputData InputDataFromCommandLine(int, char const **);

std::tuple< BoundaryPatches::ENUMDATA, BoundaryPatches::ENUMDATA > ReadSweepDirections( const std::string & );


} // end namespace CAMIRA

#endif  // INPUT_PROCESSING