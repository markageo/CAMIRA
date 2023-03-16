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

    namespace Fields
    {
        enum ENUMDATA
        {
            U, // x velocity
            V, // y velocity
            W, // z velocity
            P, // Pressure
            count
        };
    }

    namespace BoundaryConditions 
    {
        enum ENUMDATA 
        {
            zeroGradient,
            uniform,
            extrapolated,
            count
        };
    };

    namespace BoundaryPatches 
    {   
        enum ENUMDATA
        {
            xPositive, xNegative,
            yPositive, yNegative,
            zPositive, zNegative,
            count
        };
    };

    namespace TransportCoefficients
    {
        enum ENUMDATA
        {
            p,  // (i  , j  , k  )
            n,  // (i  , j+1, k  )
            e,  // (i+1, j  , k  )
            s,  // (i  , j-1, k  )
            w,  // (i-1, j  , k  )
            t,  // (i  , j  , k+1)
            b,  // (i  , j  , k-1)
            nn, // (i  , j+2, k  )
            ee, // (i+2, j  , k  )
            ss, // (i  , j-2, k  )
            ww, // (i-2, j  , k  )
            tt, // (i  , j  , k+2)
            bb, // (i  , j  , k-2)
            count
        };
    }

}


#endif // SIMULATION_PARAMETERS 