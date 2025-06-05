#ifndef SOLVER
#define SOLVER

#include "../IO/InputProcessing.h"
#include "../Core/Types.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"

namespace CAMIRA
{

    template< MomentumInterpolation >
    void SolveSteady( const InputData &, const AxisTransformationMap &);

    template< MomentumInterpolation >
    void SolveTransient( const InputData &, const AxisTransformationMap &);

}   // end namespace CAMIRA    


#endif // SOLVER