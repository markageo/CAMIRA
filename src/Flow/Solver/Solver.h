#ifndef CAMIRA_SOLVER
#define CAMIRA_SOLVER

#include "Core/Types.h"
#include "Flow/InputProcessing/InputProcessing.h"
#include "Flow/FiniteVolume/FiniteVolume.h"
#include "Flow/CoordinateTransformations/AxisTransformationFunctions.h"

namespace CAMIRA
{

namespace FLOW
{

    template< MomentumInterpolation >
    void SolveSteady( const InputData &, const AxisTransformationMap &);

    template< MomentumInterpolation >
    void SolveTransient( const InputData &, const AxisTransformationMap &);

}   // end namespace FLOW

}   // end namespace CAMIRA    


#endif // CAMIRA_SOLVER