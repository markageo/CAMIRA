#ifndef CAMIRA_PTREE_TRANSLATORS
#define CAMIRA_PTREE_TRANSLATORS

#include "Core/IO/IOTools.h"
#include "Boost/boost/property_tree/ptree.hpp"
#include <vector>
#include <string>


#define VECTOR_START_CHAR       '('
#define VECTOR_END_CHAR         ')'
#define VECTOR_DELIMITER_CHAR   ','
#define MULTI_DELIMITER_CHAR    ','


namespace pt = boost::property_tree;

// Parse vector string into an std::vector
template<typename T>
std::vector<T> ParseVectorString( const std::string &vecString )
{
    std::vector<T> vec;
    std::string::const_iterator stringIterator = vecString.begin();
    std::string valueString;

    if (*stringIterator != VECTOR_START_CHAR) {
        throw std::runtime_error( "'" + vecString + "' is not a valid input. Expected a vector." );
    }
    ++stringIterator;
    
    while ( stringIterator != vecString.end() ) {
        if ( *stringIterator == VECTOR_END_CHAR ) {
            vec.push_back( CAMIRA::CORE::IOTOOLS::String2Type<T>(valueString) );
            break;
        }

        if ( *stringIterator == VECTOR_DELIMITER_CHAR ) {
            vec.push_back( CAMIRA::CORE::IOTOOLS::String2Type<T>(valueString) );
            valueString.clear();
        } else {
            valueString += *stringIterator;
        }
        stringIterator++;
    }

    if (*stringIterator != VECTOR_END_CHAR) {
        throw std::runtime_error( "'" + vecString + "' vector not closed." );
    }

    return vec;
}


template< typename T >
struct VectorTranslator
{
    typedef std::string     internal_type;
    typedef std::vector<T>  external_type;

    external_type get_value( internal_type const &s ) {
        return ParseVectorString<T>( s );
    }
};


 // Specialization allows the translator to be used with ptree internally
namespace boost { namespace property_tree {

    template< typename T > 
    struct translator_between< std::string, std::vector< T > > 
    {
        typedef VectorTranslator< T > type;
    };

}   // end namespace property_tree  
}   // end namespace boost



#endif  // CAMIRA_PTREE_TRANSLATORS