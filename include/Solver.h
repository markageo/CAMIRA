#ifndef SOLVER
#define SOLVER

#include "InputProcessing.h"
#include "Types.h"
#include "FiniteVolume.h"
#include "SweepTransformations.h"

namespace CFD
{

    template< MomentumInterpolation MI >
    void SweepSolve(FieldData<array3D> &, const Mesh &, const FieldData< BoundaryConditionData > &, const InputData &, const AxisTransformationMap &);

}   // end namespace CFD    


#endif // SOLVER