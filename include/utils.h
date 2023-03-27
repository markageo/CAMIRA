#ifndef UTILS
#define UTILS

#include "Types.h"
#include "Tensor"

#include <fstream>
#include <iomanip>
#include <iostream>

namespace UTIL
{

// Write eigen tensor to a file for debugging
template<typename T>
void writeArray(const std::string &filename, const T &array, const int precision = 4)
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