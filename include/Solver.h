#ifndef SOLVER
#define SOLVER

#include "InputProcessing.h"
#include "Types.h"
#include "FiniteVolume.h"

namespace CFD
{

    // struct ConvergenceData
    // {
    //     std::vector< EnumVector<Fields, floatType> > residuals;

    //     void reserve( std::size_t capacity )
    //     {
    //         residuals.reserve( capacity );
    //     }

    //     void shrink_to_fit()
    //     {
    //         residuals.shrink_to_fit();
    //     }
    // };

    // ConvergenceData SweepSolve(CFD::ArrayAllocator<CFD::Fields, CFD::array3D> &, const CFD::Mesh &, const CFD::InputData &);
    void SweepSolve(CFD::ArrayAllocator<CFD::Fields, CFD::array3D> &, const CFD::Mesh &, const CFD::InputData &);

}   // end namespace CFD    


#endif // SOLVER