#include "utils.h"

#include "Types.h"
#include <fstream>
#include <iomanip>

namespace UTIL
{

using namespace CFD;

// Write an an array to file in an easy to read way, used for debugging
void writeArray(const std::string &filename, const array3D &array, const int precision)
{

    // Open file and set precision
    std::ofstream fileStream(filename);
    fileStream << std::fixed;
    fileStream << std::setprecision(precision);


    // Write each plane in a seperate block
    for (int k = 0; k != array.dimension(2); k++) {

        fileStream << "k = " << k << "\n";

        for (int i = 0; i != array.dimension(0); i++) {
            
            for (int j = 0; j != array.dimension(1); j++) {

                fileStream << array(i, j, k);
                if (j != array.dimension(1)-1)
                    fileStream << ", ";
                
            }
            fileStream << "\n";

        }
        fileStream << "\n";

    }

}

}