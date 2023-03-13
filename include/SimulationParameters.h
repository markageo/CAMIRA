#ifndef SIMULATION_PARAMETERS
#define SIMULATION_PARAMETERS

#include "Tensor"

namespace CFD 
{
    using floatType = double;
    using intType = int;
    using array1D = Eigen::Tensor<floatType, 1>;    // Column major
    using array2D = Eigen::Tensor<floatType, 2>;    // Column major
    using array3D = Eigen::Tensor<floatType, 3>;    // Column major
}


#endif // SIMULATION_PARAMETERS 