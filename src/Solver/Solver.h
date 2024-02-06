#ifndef SOLVER
#define SOLVER

#include "../IO/InputProcessing.h"
#include "../Types.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Tools/SweepTransformations.h"

namespace CFD
{

    template< MomentumInterpolation, Linearisation >
    void SweepSolve(FieldData<Tensor3D> &, const Mesh &, const BoundaryConditionData &, const InputData &, const AxisTransformationMap &);

}   // end namespace CFD    


#endif // SOLVER