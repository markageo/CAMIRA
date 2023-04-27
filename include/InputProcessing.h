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
    InputData();

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


    // Structure for storing axis transformation, is just a one-to-one map
    class AxisTransformationMap
    {
        using BP = BoundaryPatches::ENUMDATA;
        public:
            AxisTransformationMap();

            // Code patch from user patch
            BP &CodePatch(const BP userPatch)
            { return m_userMap[ userPatch ]; }

            const BP &CodePatch(const BP userPatch) const 
            { return m_userMap.at( userPatch ); }

            // User patch from the code patch
            BP &UserPatch(const BP codePatch)
            { return m_codeMap[ codePatch ]; }

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