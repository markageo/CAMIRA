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


// Return relative path to given filename (LINUX FILESYSTEMS ONLY)
std::string RelativePath( const std::string &filename ) 
{ 
    std::string::const_iterator stringIterator = filename.end();
    for ( /* NULL */ ; *stringIterator != '/' && stringIterator != filename.begin(); stringIterator--) {};
    return std::string( filename.begin(), stringIterator );
}


}   // end namespace IO


#endif // IO_UTILS