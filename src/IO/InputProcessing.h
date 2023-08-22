#ifndef INPUT_PROCESSING
#define INPUT_PROCESSING

#include "../Types.h"
#include "boost/property_tree/ptree.hpp"
#include <vector>
#include <utility>
#include <map>
#include <tuple>

namespace CFD
{


// -------------------------------------- Definition in InputProcessing.cpp -------------------------------------- //


struct InputData
{
    // Constructor
    InputData();

    // File directory of input file
    std::string inputFileDirectory;

    // Model
    floatType rho, nu;


    // Domain size
    fArray3 domainSize;


    // Mesh
    struct MeshSegment {
        floatType lowerBound;
        floatType upperBound;
        intType nCells;
        floatType biasFactor;
    };
    EnumVector< Axis, std::vector< MeshSegment > > meshSegments;

    // Solid Geometry
    struct SolidBlockData {
        fArray3 centerPosition, 
                 dimensions, 
                 rotation;
    };
    std::vector< SolidBlockData > solidBlocks;

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
    FieldData<floatType> initialConditions;

    // Solver
    struct Schemes {
        Linearisation linearisation;
        MomentumInterpolation momentumInterpolation;
        AdvectionSchemes advectionScheme;
        FaceInterpolationSchemes faceInterpolationScheme;

        FieldData<floatType> implicitRelaxation;

        intType maxOuterIterations;
        FieldData<floatType> maxOuterResiduals;
    };
    Schemes schemes;

    struct LinearSolverSettings {
        LinearSolvers type;
        intType maxIterations;
        FieldData<floatType> maxResiduals;
        FieldData<floatType> relaxation;

        BoundaryPatches::ENUMDATA lineSweepDirection;
        BoundaryPatches::ENUMDATA planeSweepDirection;
    }; 
    LinearSolverSettings linearSolverSettings;


    // Ouput
    std::string residualHistoryFilename;
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

std::tuple< BoundaryPatches::ENUMDATA, BoundaryPatches::ENUMDATA > ReadSweepDirections( const std::string & );


} // end namespace CFD

#endif  // INPUT_PROCESSING