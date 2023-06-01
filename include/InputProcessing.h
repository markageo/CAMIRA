#ifndef INPUT_PROCESSING
#define INPUT_PROCESSING

#include "Types.h"
#include "boost/property_tree/ptree.hpp"
#include <vector>
#include <utility>
#include <map>

namespace CFD
{

struct InputData
{
    // Constructor
    Inputcd ..Data();

    // Model
    floatType rho, nu;


    // Domain size
    floatVector3 domainSize;


    // Mesh
    struct MeshSegment {
        floatType lowerBound;
        floatType upperBound;
        intType nCells;
        floatType biasFactor;
    };
    EnumVector< Axis, std::vector< MeshSegment > > meshSegments;


    // Boundary conditions
    struct BoundaryConditionStruct {
        BoundaryConditions::ENUMDATA type;
        CFD::floatType value;
    };
    using BoundaryConditionData = EnumVector< Fields, EnumVector< BoundaryPatches, BoundaryConditionStruct > >;
    BoundaryConditionData boundaryConditions;

    // Initial conditions
    EnumVector<Fields, floatType> initialConditions;

    // Solver
    struct Schemes {
        Linearisation linearisation;
        EnumVector<Fields, floatType> implicitRelaxation;

        intType maxOuterIterations;
        EnumVector<Fields, floatType> maxOuterResiduals;
    } schemes;

    struct LinearSolverSettings {
        LinearSolvers type;
        intType maxIterations;
        EnumVector<Fields, floatType> maxResiduals;
        EnumVector<Fields, floatType> relaxation;
    } linearSolverSettings;

    struct PlaneSolverSettings {
        PlaneSolvers type;
        intType maxIterations;
        EnumVector<Fields, floatType> maxResiduals;
        EnumVector<Fields, floatType> relaxation;
        BoundaryPatches::ENUMDATA sweepDirection;
    } planeSolverSettings;

    struct LineSolverSettings {
        LineSolvers type;
        intType maxIterations;
        EnumVector<Fields, floatType> maxResiduals;
        EnumVector<Fields, floatType> relaxation;
        BoundaryPatches::ENUMDATA sweepDirection;
    } lineSolverSettings;


    // Structure for storing axis transformation, is just a one-to-one map
    class AxisTransformationMap
    {
        using BP = BoundaryPatches::ENUMDATA;
        public:
            AxisTransformationMap();

            // Setting values
            void Set(const BP codePatch, const BP userPatch)
            {
                m_codeMap[codePatch] = userPatch;
                m_userMap[userPatch] = codePatch; 
            }

            // Code patch from user patch
            const BP &CodePatch(const BP userPatch) const 
            { return m_userMap.at( userPatch ); }

            // User patch from the code patch
            const BP &UserPatch(const BP codePatch) const
            { return m_codeMap.at( codePatch ); }

        private:
            std::map< BP, BP> m_codeMap;    // Code patch -> user patch
            std::map< BP, BP> m_userMap;    // User patch -> code patch
    };
    AxisTransformationMap axisTransformation;

};



InputData ReadInputData(const std::string &);

InputData InputDataFromCommandLine(int, char const **);


} // end namespace CFD

#endif  // INPUT_PROCESSING