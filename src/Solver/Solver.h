#ifndef SOLVER
#define SOLVER

#include "../IO/InputProcessing.h"
#include "../Types.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Tools/SweepTransformations.h"

namespace CFD
{

    template< MomentumInterpolation, Linearisation >
    void SolveSteady( const InputData &, const AxisTransformationMap &);

    template< MomentumInterpolation, Linearisation >
    void SolveTransient( const InputData &, const AxisTransformationMap &);

}   // end namespace CFD    


#endif // SOLVER