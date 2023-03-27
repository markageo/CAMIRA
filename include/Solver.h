#ifndef SOLVER
#define SOLVER

#include "InputProcessing.h"
#include "Types.h"
#include "FiniteVolumeStructures.h"

namespace CFD
{

    void SweepSolve(CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &, const CFD::Mesh &, const CFD::InputData &);

}   // end namespace CFD    


#endif // SOLVER