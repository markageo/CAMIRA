#ifndef INPUT_PARSER
#define INPUT_PARSER

#include "boost/property_tree/ptree.hpp"
#include <string>
#include <optional>

// Parsing function declaration
std::optional<boost::property_tree::ptree> ReadInput(const std::string &inputFileName) ;


#endif  // INPUT_PARSER