#include "InputProcessing.h"

#define NDEBUG
#include "InputParser.h"
#include "SimulationParameters.h"
#include "boost/property_tree/ptree.hpp"
#include <utility>
#include <optional>
#include <iostream>
#include <algorithm>

#define VECTOR_START_CHAR       '('
#define VECTOR_END_CHAR         ')'
#define VECTOR_DELIMITER_CHAR   ','

namespace pt = boost::property_tree;

namespace
{

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
        std::vector<SIM::floatType> tempBoundsVector;
        for (auto segment : gridTree) {
            if (segment.first != "Segment") {
                // throw ERROR invalid child name
            }
            nCellsString     = segment.second.get<std::string>("nCells");
            boundsString     = segment.second.get<std::string>("bounds");
            biasFactorString = segment.second.get<std::string>("biasFactor");

            tempBoundsVector = ParseVectorString(boundsString, 2);
            meshSegments.emplace_back( tempBoundsVector[0],
                                       tempBoundsVector[1],
                                       String2Type<SIM::floatType>(nCellsString),
                                       String2Type<SIM::floatType>(biasFactorString) 
                                     );
        }
    }

    void ReadMesh(InputData &inputData, const pt::ptree &tree)
    {
        const pt::ptree &meshTree = tree.get_child("Mesh");

        // Domain
        const std::string &domainSizeString = meshTree.get<std::string>("domain");
        std::vector<SIM::floatType> domainSize = ParseVectorString(domainSizeString, 3);
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




