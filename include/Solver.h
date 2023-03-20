#ifndef SOLVER
#define SOLVER

#include "InputProcessing.h"
#include "Types.h"
#include "FiniteVolumeStructures.h"

namespace CFD
{
    using namespace CFD;

    void SweepSolve(ArrayAllocator<CFD::Fields::ENUMDATA> &, const Mesh &, const InputData &);

}   // end namespace CFD    


#endif // SOLVER