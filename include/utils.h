#ifndef UTILS
#define UTILS

#include "Types.h"

namespace UTIL
{

// Write an an array to file in an easy to read way, used for debugging
void writeArray(const std::string &, const CFD::array3D &, const int = 4);

}


#endif // UTILS