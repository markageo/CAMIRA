#ifndef UTILS
#define UTILS

#include "Types.h"
#include <unsupported/Eigen/CXX11/Tensor>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#ifdef PROFILING
#   if !defined(TIC) || !defined(TOC)
#       include "profiler/profiler.h"
#       define TIC(name) PROF::prof.tic(name);
#       define TOC(name) PROF::prof.toc(name);
#   endif

    namespace PROF {
        inline profiler<perf_counter::clock<time_units::SECONDS>> prof;
    }

#else
#   ifndef TIC
#       define TIC(name)
#   endif
#   ifndef TOC
#       define TOC(name)
#   endif
#endif


namespace UTIL
{

// Write eigen tensor to a file for debugging
template<typename T>
void WriteArray(const std::string &filename, const T &array, const int precision = 6)
{
    static_assert(std::is_same< T, CFD::array1D >::value ||
                  std::is_same< T, CFD::array2D >::value ||
                  std::is_same< T, CFD::array3D >::value,
                  "Array type invalid.");

    std::ofstream fileStream(filename);
    fileStream << std::fixed;
    fileStream << std::setprecision(precision);

    if (array.size() > 0) {

        const auto& dims = array.dimensions();

        if (dims.size() == 1) {
            fileStream << "dimensions" << " "
                       << dims[0]      << "\n\n";

            fileStream << array;

        } else if (dims.size() == 2) {
            fileStream << "dimensions" << " "
                       << dims[0]      << " "
                       << dims[1]      << "\n\n";

            fileStream << array;

        } else if (dims.size() == 3) {
            fileStream << "dimensions" << " "
                       << dims[0]      << " "
                       << dims[1]      << " "
                       << dims[2]      << "\n\n";

            for (int k = 0; k != dims[2]; k++) {
                fileStream << array.chip(k , 2);
                fileStream << "\n\n";
            }

        }

    }
}


// Functions used for constructing arrays in ReadArray function
namespace 
{

    // Template specialization for constructing array
    template<typename T> inline
    T ConstructArray(const std::array<CFD::intType, 3> &dims) = delete;

    template<> inline
    CFD::array1D ConstructArray<CFD::array1D>(const std::array<CFD::intType, 3> &dims)
    {
        return CFD::array1D(dims[0]);
    }

    template<> inline
    CFD::array2D ConstructArray<CFD::array2D>(const std::array<CFD::intType, 3> &dims)
    {
        return CFD::array2D(dims[0], dims[1]);
    }

    template<> inline
    CFD::array3D ConstructArray<CFD::array3D>(const std::array<CFD::intType, 3> &dims)
    {
        return CFD::array3D(dims[0], dims[1], dims[2]);
    }

}


// Read Eigen::Tensor (1D, 2D, or 3D) from file, file must be in the format created by WriteArray
template<typename T>
T ReadArray(const std::string &filename)
{

    static_assert(std::is_same< T, CFD::array1D >::value ||
                  std::is_same< T, CFD::array2D >::value ||
                  std::is_same< T, CFD::array3D >::value,
                  "Array type invalid.");

    std::ifstream fileStream(filename);

    // Get dimensions from first line
    std::string line, word;                  // Some temporary variables
    CFD::intType dim, ndims = 0;
    std::array<CFD::intType, 3> dims = {1, 1, 1}; // Unused dimensions have size 1   
    std::getline(fileStream, line);          // First line of file contains dimension information
    std::istringstream lineStream(line);     // For tokenizing line and casting 
    lineStream >> word;
    if (word != "dimensions") {              // File is invalid
        T emptyArray;
        return emptyArray;   
    }
    while (lineStream >> dim ) {             // Read in dimensions
        dims[ static_cast<size_t>( ndims ) ] = dim;
        ndims++;
    }

    // Array to return
    T array( ConstructArray<T>(dims) );
 
    // Read into the tensor, tensor must be column major, stride by number of columns to account 
    // for reading order not being continguous.
    auto *dataPointer = array.data();
    CFD::intType idx;

    for (CFD::intType k = 0; k != dims[2]; k++) {
        for (CFD::intType i = 0; i != dims[0]; i++) {
            for (CFD::intType j = 0; j != dims[1]; j++) {

                idx = i + j*dims[0] + k*dims[0]*dims[1];
                fileStream >> dataPointer[idx];

            }
        }
    }

    return array;
}



}   // end namespace UTIL


#endif // UTILS