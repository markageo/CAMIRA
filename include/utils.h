#ifndef UTILS
#define UTILS

#include "Types.h"
#include "Eigen/unsupported/CXX11/Tensor"

#include <fstream>
#include <iomanip>
#include <iostream>

#ifdef PROFILING
#   if !defined(TIC) || !defined(TOC)
#       include <profiler/profiler.h>
#       define TIC(name) PROF::prof.tic(name);
#       define TOC(name) PROF::prof.toc(name);
#   endif
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


// Read eigen tensor from file, file must be in the format created by WriteArray
template<typename T>
T ReadArray(const std::string &filename)
{
    static_assert(std::is_same< T, CFD::array1D >::value ||
                  std::is_same< T, CFD::array2D >::value ||
                  std::is_same< T, CFD::array3D >::value,
                  "Array type invalid.");

    

}


}   // end namespace UTIL


#endif // UTILS