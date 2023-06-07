#ifndef SOLVER
#define SOLVER

#include "InputProcessing.h"
#include "Types.h"
#include "FiniteVolume.h"
#include "SweepTransformations.h"

namespace CFD
{

    void SweepSolve(CFD::ArrayAllocator<CFD::Fields, CFD::array3D> &, const CFD::Mesh &, const CFD::InputData &, const AxisTransformationMap &);

}   // end namespace CFD    


#endif // SOLVER