#include "InputProcessing.h"

#define NDEBUG
#include "InputParser.h"
#include "SimulationParameters.h"
#include "boost/property_tree/ptree.hpp"
#include <utility>
#include <optional>
#include <iostream>

#define VECTOR_START_CHAR       '('
#define VECTOR_END_CHAR         ')'
#define VECTOR_DELIMITER_CHAR   ','

namespace pt = boost::property_tree;

namespace
{

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
    std::vector<SIM::floatType> ParseVectorString(const std::string &vecString, const int &dim)
    {
        std::vector<SIM::floatType> vec;
        std::string::const_iterator stringIterator = vecString.begin();
        std::string valueString;
        int dimCount = 0;

        if (*stringIterator != VECTOR_START_CHAR) {
            // throw ERROR - expecting vector
        }
        ++stringIterator;
        
        while(stringIterator != vecString.end()) {
            if (*stringIterator == VECTOR_END_CHAR) {
                vec.push_back( String2Type<SIM::floatType>(valueString) );
                break;
            }

            if (*stringIterator == VECTOR_DELIMITER_CHAR) {
                vec.push_back( String2Type<SIM::floatType>(valueString) );
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

    void ReadGrid(const pt::ptree &meshTree, std::vector<SIM::floatType> &nCells, std::vector<SIM::floatType> &biasFactors, 
        std::vector<std::pair<SIM::floatType, SIM::floatType>> &segmentBounds, const std::string &gridString)
    {
        const pt::ptree &gridTree = meshTree.get_child(gridString);
        std::string boundsString, nCellsString, biasFactorString;
        std::pair<SIM::floatType, SIM::floatType> boundsPair;
        std::vector<SIM::floatType> tempBoundsVector;
        for (auto segment : gridTree) {
            if (segment.first != "Segment") {
                // throw ERROR invalid child name
            }
            nCellsString     = segment.second.get<std::string>("nCells");
            boundsString     = segment.second.get<std::string>("bounds");
            biasFactorString = segment.second.get<std::string>("biasFactor");

            nCells.push_back( String2Type<SIM::floatType>(nCellsString) );
            biasFactors.push_back( String2Type<SIM::floatType>(biasFactorString) );

            tempBoundsVector = ParseVectorString(boundsString, 2);
            boundsPair.first = tempBoundsVector[0];
            boundsPair.second = tempBoundsVector[1];
            segmentBounds.push_back( boundsPair );
        }
    }


    void ReadMesh(InputData &inputData, const pt::ptree &tree)
    {
        const pt::ptree &meshTree = tree.get_child("Mesh");

        // Domain
        const std::string &domainSizeString = meshTree.get<std::string>("domain");
        std::vector<SIM::floatType> domainSize = ParseVectorString(domainSizeString, 3);
        inputData.domainSizeX = domainSize[0];
        inputData.domainSizeY = domainSize[1];
        inputData.domainSizeZ = domainSize[2];

        // Grid X
        ReadGrid(meshTree, inputData.mesh.xCells, inputData.mesh.xBiasFactors, inputData.mesh.xSegmentBounds, "GridX");

        // Grid Y
        ReadGrid(meshTree, inputData.mesh.yCells, inputData.mesh.yBiasFactors, inputData.mesh.ySegmentBounds, "GridY");

        // Grid Z
        ReadGrid(meshTree, inputData.mesh.zCells, inputData.mesh.zBiasFactors, inputData.mesh.zSegmentBounds, "GridZ");
        
    }


}   // end anonymous namepsace


std::optional<InputData> ReadInputData(const std::string &inputFileName) 
{
    std::optional<pt::ptree> tree_optional = ParseFile(inputFileName);
    if (!tree_optional)
        return {};
    pt::ptree tree = tree_optional.value();

    InputData inputData;
    ReadMesh(inputData, tree);
    
    return inputData;

}




