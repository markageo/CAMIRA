#ifndef SOLVER
#define SOLVER

#include "InputProcessing.h"
#include "Types.h"
#include "FiniteVolume.h"

namespace CFD
{

    void SweepSolve(CFD::ArrayAllocator<CFD::Fields::ENUMDATA, CFD::array3D> &, const CFD::Mesh &, const CFD::InputData &);

}   // end namespace CFD    


#endif // SOLVER