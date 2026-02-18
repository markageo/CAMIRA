#ifndef CAMIRA_INPUT_PARSER
#define CAMIRA_INPUT_PARSER

#include "boost/property_tree/ptree.hpp"
#include <string>
#include <stdexcept>

// Parsing function declaration
namespace INP
{

boost::property_tree::ptree ParseFile(const std::string &inputFileName) ;

}   // end namespace INP



#endif  // CAMIRA_INPUT_PARSER