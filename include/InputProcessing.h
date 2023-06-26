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


    // A general struct for holding values corresponding to different fields
    template < typename dataType >
    struct FieldData {
        EnumVector<Axis, dataType> U;
        dataType P;
    };

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

    struct LineSolverSettings {
        LineSolvers type;
        intType maxIterations;
        FieldData<floatType> maxResiduals;
        FieldData<floatType> relaxation;
        BoundaryPatches::ENUMDATA sweepDirection;
    };

    struct PlaneSolverSettings {
        PlaneSolvers type;
        intType maxIterations;
        FieldData<floatType> maxResiduals;
        FieldData<floatType> relaxation;
        BoundaryPatches::ENUMDATA sweepDirection;

        LineSolverSettings lineSolverSettings;
    };

    
    struct LinearSolverSettings {
        LinearSolvers type;
        intType maxIterations;
        FieldData<floatType> maxResiduals;
        FieldData<floatType> relaxation;

        PlaneSolverSettings planeSolverSettings;
    }; 
    LinearSolverSettings linearSolverSettings;

};


InputData ReadInputData(const std::string &);

InputData InputDataFromCommandLine(int, char const **);


} // end namespace CFD

#endif  // INPUT_PROCESSING