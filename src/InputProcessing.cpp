#include "InputProcessing.h"

#define NDEBUG
#include "InputParser.h"
#include "FiniteVolumeStructures.h"
#include "Types.h"
#include "boost/property_tree/ptree.hpp"
#include <utility>
#include <optional>
#include <iostream>
#include <algorithm>
#include <map>

#define VECTOR_START_CHAR       '('
#define VECTOR_END_CHAR         ')'
#define VECTOR_DELIMITER_CHAR   ','
#define MULTI_DELIMITER_CHAR    ','

namespace pt = boost::property_tree;

// InputData constructor
CFD::InputData::InputData() :
    boundaryConditions(CFD::Fields::ENUMDATA::count, std::vector< CFD::InputData::BoundaryConditionStruct >(CFD::BoundaryPatches::ENUMDATA::count) )
    {};


namespace
{

    using namespace CFD;

    /*-------------------------------------------------------------------------------------*\
                                         Helper Functions
    \*-------------------------------------------------------------------------------------*/

    // Convert string to given numeric type T.
    template <typename T> T 
    String2Type(const std::string &str)
    {
        // NOTE: This does not work for ints in scientific notation.
        std::istringstream strstream(str);
        T num;
        strstream >> num;
        return num;
    }

    // Parse vector string into an std::vector
    std::vector<CFD::floatType> ParseVectorString(const std::string &vecString, const int &dim)
    {
        std::vector<CFD::floatType> vec;
        std::string::const_iterator stringIterator = vecString.begin();
        std::string valueString;
        int dimCount = 0;

        if (*stringIterator != VECTOR_START_CHAR) {
            // throw ERROR - expecting vector
        }
        ++stringIterator;
        
        while(stringIterator != vecString.end()) {
            if (*stringIterator == VECTOR_END_CHAR) {
                vec.push_back( String2Type<CFD::floatType>(valueString) );
                break;
            }

            if (*stringIterator == VECTOR_DELIMITER_CHAR) {
                vec.push_back( String2Type<CFD::floatType>(valueString) );
                valueString.clear();
                dimCount++;
            } else {
                valueString += *stringIterator;
            }
            stringIterator++;
        }

        if (*stringIterator != VECTOR_END_CHAR) {
            // throw ERROR - vector not closed
        }

        if (dimCount != dim) {
            // throw ERROR - expecting vector of size dim
        }

        return vec;

    }

    /*-------------------------------------------------------------------------------------*\
                                                Mesh
    \*-------------------------------------------------------------------------------------*/

    void ValidateSegmentBounds(const std::vector<InputData::MeshSegment> &meshSegment)
    {
        using v_size_type = std::vector<InputData::MeshSegment>::size_type;

        // The first one must start at zero
        if (meshSegment[0].lowerBound != 0.0) {
            // throw ERROR - must have mesh segment starting at coordinate 0
        }

        // Must be no gaps or overlaps in the segments
        for (v_size_type i = 1; i != meshSegment.size(); i++) {
            if ( meshSegment[i].lowerBound != meshSegment[i-1].upperBound ) {
                // throw ERROR - segments must not overlap or have gaps
            }
        }

    }


    void ReadGrid(const pt::ptree &meshTree, std::vector<InputData::MeshSegment> &meshSegments, const std::string &gridString)
    {
        const pt::ptree &gridTree = meshTree.get_child(gridString);
        std::string boundsString, nCellsString, biasFactorString;
        std::vector<CFD::floatType> tempBoundsVector;
        InputData::MeshSegment tempMeshSegment;
        for (auto segment : gridTree) {
            if (segment.first != "Segment") {
                // throw ERROR invalid child name
            }
            nCellsString     = segment.second.get<std::string>("nCells");
            boundsString     = segment.second.get<std::string>("bounds");
            biasFactorString = segment.second.get<std::string>("biasFactor");
            tempBoundsVector = ParseVectorString(boundsString, 2);

            tempMeshSegment.nCells = String2Type<CFD::floatType>(nCellsString);
            tempMeshSegment.biasFactor = String2Type<CFD::floatType>(biasFactorString);
            tempMeshSegment.lowerBound = tempBoundsVector[0];
            tempMeshSegment.upperBound = tempBoundsVector[1];

            meshSegments.push_back(tempMeshSegment);
        }
    }


    void ReadMesh(InputData &inputData, const pt::ptree &tree)
    {
        const pt::ptree &meshTree = tree.get_child("Mesh");

        // Domain
        const std::string &domainSizeString = meshTree.get<std::string>("domain");
        std::vector<CFD::floatType> domainSize = ParseVectorString(domainSizeString, 3);
        inputData.domainSize_x = domainSize[0];
        inputData.domainSize_y = domainSize[1];
        inputData.domainSize_z = domainSize[2];

        // Grids
        ReadGrid(meshTree, inputData.meshSegments_x, "GridX");
        ReadGrid(meshTree, inputData.meshSegments_y, "GridY");
        ReadGrid(meshTree, inputData.meshSegments_z, "GridZ");

        // Sort in increasing order of lower bound
        auto sortComparison = [](const auto& i, const auto& j) { return i.lowerBound < j.lowerBound; };
        std::sort( inputData.meshSegments_x.begin(), inputData.meshSegments_x.end(), sortComparison);
        std::sort( inputData.meshSegments_y.begin(), inputData.meshSegments_y.end(), sortComparison );
        std::sort( inputData.meshSegments_z.begin(), inputData.meshSegments_z.end(), sortComparison );

        // Check that the bounds are valid
        ValidateSegmentBounds(inputData.meshSegments_x);
        ValidateSegmentBounds(inputData.meshSegments_y);
        ValidateSegmentBounds(inputData.meshSegments_z);

    }

    /*-------------------------------------------------------------------------------------*\
                                       Boundary Conditions
    \*-------------------------------------------------------------------------------------*/

    InputData::BoundaryConditionStruct readBoundaryValueString(const std::string &boundaryString)
    {
        using BC = BoundaryConditions::ENUMDATA;

        InputData::BoundaryConditionStruct bcStruct;
        std::string bcTypeString;
        std::string bcValueString;
        std::string::const_iterator stringIterator = boundaryString.begin();

        // Read the boundary condition type
        while (stringIterator != boundaryString.end()) {

            if (*stringIterator == MULTI_DELIMITER_CHAR) 
                break;
            
            bcTypeString += *stringIterator;
            stringIterator++;

        }


        // Store the type, some don't need a value
        if        (bcTypeString == "uniform") {

            bcStruct.type = BC::uniform;

        } else if (bcTypeString == "zeroGradient") {

            bcStruct.type = BC::zeroGradient;
            bcStruct.value = 0.0;
            return bcStruct;

        } else if (bcTypeString == "extrapolated") {
            bcStruct.type = BC::extrapolated;
            bcStruct.value = 0.0;
            return bcStruct;
        } else {
            // throw ERROR -  invalid BC type
        }


        // Check for expected value
        if (*stringIterator != MULTI_DELIMITER_CHAR) {
            // throw ERROR - expected value for boudnary condition
        }
        stringIterator++;

        // Read the value to end of line
        while (stringIterator != boundaryString.end()) {            
            bcValueString += *stringIterator;
            stringIterator++;
        }

        // Convert and store value
        if (bcStruct.type == BC::uniform) {
            bcStruct.value = std::stod(bcValueString);
        }

        return bcStruct;
    }


    void ReadBoundaryConditions(InputData &inputData, const pt::ptree &tree)
    {
        using BP = BoundaryPatches::ENUMDATA;
        using F  = Fields::ENUMDATA;

        const pt::ptree &boundaryConditionsTree = tree.get_child("BoundaryConditions");

        // Boundary name strings
        std::map<BP, std::string> boundaryPatchMap{ {BP::xPositive, "+x"}, {BP::xNegative, "-x"},
                                                    {BP::yPositive, "+y"}, {BP::yNegative, "-y"},
                                                    {BP::zPositive, "+z"}, {BP::zNegative, "-z"} };

        // Field name strings
        std::map<F, std::string> fieldMap { {F::U, "u"},
                                            {F::V, "v"},
                                            {F::W, "w"},
                                            {F::P, "p"} };

        // Iterate boundary patches
        std::string valueString;
        const pt::ptree *boundaryPatchTreePointer = nullptr;
        for (auto [patchEnum, patchString] : boundaryPatchMap) {
            boundaryPatchTreePointer = &( boundaryConditionsTree.get_child(patchString) );

            // Iterate fields
            for (auto [fieldEnum, fieldString] : fieldMap) {

                valueString = boundaryPatchTreePointer->get<std::string>(fieldString);
                inputData.boundaryConditions[fieldEnum][patchEnum] = readBoundaryValueString(valueString);

            }
        }

    }


}   // end anonymous namepsace


CFD::InputData CFD::ReadInputData(const std::string &inputFileName) 
{
    pt::ptree tree = ParseFile(inputFileName);
    
    InputData inputData;
    ReadMesh(inputData, tree);
    ReadBoundaryConditions(inputData, tree);
    
    return inputData;

}




