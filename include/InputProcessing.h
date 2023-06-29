#ifndef INPUT_PROCESSING
#define INPUT_PROCESSING

#include "Types.h"
#include "boost/property_tree/ptree.hpp"
#include <vector>
#include <utility>
#include <map>

namespace CFD
{


// -------------------------------------- Definition in InputProcessing.cpp -------------------------------------- //


struct InputData
{
    // Constructor
    InputData();

    // Model
    floatType rho, nu;


    // Domain size
    fVector3 domainSize;


    // Mesh
    struct MeshSegment {
        floatType lowerBound;
        floatType upperBound;
        intType nCells;
        floatType biasFactor;
    };
    EnumVector< Axis, std::vector< MeshSegment > > meshSegments;

    // Boundary conditions
    struct BoundaryConditionData {
        BoundaryConditions::ENUMDATA type;
        floatType value;    
    };
    using FieldBoundaryConditions =  FieldData< EnumVector< BoundaryPatches, BoundaryConditionData  > >;
    FieldBoundaryConditions boundaryConditions;


    // Initial conditions
    FieldData<floatType> initialConditions;

    // Solver
    struct Schemes {
        Linearisation linearisation;
        FieldData<floatType> implicitRelaxation;

        AdvectionSchemes advectionScheme;
        FaceInterpolationSchemes faceInterpolationScheme;

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

};


InputData ReadInputData(const std::string &);

InputData InputDataFromCommandLine(int, char const **);


} // end namespace CFD

#endif  // INPUT_PROCESSING