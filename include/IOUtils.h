#ifndef IO_UTILS
#define IO_UTILS

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

namespace IO
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


// Return string with no whitespace
std::string RemoveWhitespace( std::string str )
{
    str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
    return str;
}


}   // end namespace UTIL


#endif // IO_UTILS