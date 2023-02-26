#include "InputProcessing.h"

#define NDEBUG
#include "boost/property_tree/ptree.hpp"
#include "InputParser.h"
#include <utility>
#include <optional>

namespace pt = boost::property_tree;


/*-------------------------------------------------------------------------------------*\
                                     InputData
\*-------------------------------------------------------------------------------------*/

namespace
{
    void ReadMesh(InputData &inputData, const pt::ptree &tree)
    {
        
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




