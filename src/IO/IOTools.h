#ifndef CFD_IOTOOLS
#define CFD_IOTOOLS

#include <string>

namespace CFD 
{

namespace IOTOOLS
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
    inline std::string RemoveWhitespace( std::string str )
    {
        str.erase(remove_if(str.begin(), str.end(), isspace), str.end());
        return str;
    }


    // Return relative path to given filename (LINUX FILESYSTEMS ONLY)
    inline std::string RelativePath( const std::string &filename ) 
    { 
        std::string::const_iterator stringIterator = filename.end();
        for ( /* NULL */ ; *stringIterator != '/' && stringIterator != filename.begin(); stringIterator--) {};
        return std::string( filename.begin(), stringIterator );
    }


    // Remove file extension
    inline std::string RemoveFileExtension( const std::string &filename,
                                            const std::string &fileExtension )
    {
        // Check if there is a .vtk extension and remove it
        size_t lastPointPosition = filename.find_last_of(".");
        if ( lastPointPosition == std::string::npos ) {
            return filename;
        }
        if ( filename.substr( lastPointPosition ) == fileExtension ) {
            return filename.substr( 0, lastPointPosition );
        }
        return filename;
    }


}   // end namespace IOTOOLS

}   // end namespace CFD





#endif // CFD_IOTOOLS
