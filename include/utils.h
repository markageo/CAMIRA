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
void writeArray(const std::string &filename, const T &array, const int precision = 6)
{
    std::ofstream fileStream(filename);
    fileStream << std::fixed;
    fileStream << std::setprecision(precision);

    if (array.size() > 0) {

        const auto& dims = array.dimensions();

        if (dims.size() < 3) {
            fileStream << array;
        } else if (dims.size() == 3) {
            for (int k = 0; k != dims[2]; k++) {
                fileStream << array.chip(k , 2);
                fileStream << "\n\n";
            }
        }
    }
}

}


#endif // UTILS