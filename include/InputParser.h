#ifndef INPUT_PARSER
#define INPUT_PARSER

#include "boost/property_tree/ptree.hpp"
#include <string>

// Parsing function declaration
int ReadInput(boost::property_tree::ptree &, const std::string &); 


#endif  // INPUT_PARSER