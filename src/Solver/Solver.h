#ifndef SOLVER
#define SOLVER

#include "../IO/InputProcessing.h"
#include "../Core/Types.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"

namespace CFD
{

    template< MomentumInterpolation, Linearisation >
    void SolveSteady( const InputData &, const AxisTransformationMap &);

    template< MomentumInterpolation, Linearisation >
    void SolveTransient( const InputData &, const AxisTransformationMap &);

}   // end namespace CFD    


#endif // SOLVER