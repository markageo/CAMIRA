#ifndef PLANE_SOLVER
#define PLANE_SOLVER

#include "StaggerIndexing.h"
#include "StencilConstants.h"
#include "LineSolver.h"

#include "../Core/Types.h"
#include "../Core/Macros.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/FiniteVolume.h"

namespace CFD
{

using namespace FVT;

template <TransportCoefficients::ENUMDATA Wstag, 
          MomentumInterpolation MI >
class PlaneSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    // Staggering must be valid
    static_assert( (Wstag == TC::t) || (Wstag == TC::b) || (Wstag == TC::p), "Invalid W momentum staggering" );

public:
    PlaneSolver( FieldData<Tensor3D> &fields,
                 const FieldData<Tensor3D> &fieldsOld,
                 const Tensor3D &mask,
                 const FVCoefficients &fvCoeffs,
                 const InputData::SmootherSettings &smootherSettings) : 
                    m_fields( fields ),
                    m_fvCoeffs( fvCoeffs ),
                    m_planeConstants( Tensor2D( fvCoeffs.nCells(Axis::X) + 2*nGhost, fvCoeffs.nCells(Axis::Y) + 2*nGhost ).setZero() ),
                    m_ni( fvCoeffs.nCells(Axis::X) ),
                    m_nj( fvCoeffs.nCells(Axis::Y) )
    {
        if (m_nj == 1) {
            m_lineSolverCenter = std::make_unique<LineSolver<TC::p, Wstag, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            SolutionUpdater    = &PlaneSolver::Sweep2D;
            StateUpdater       = &PlaneSolver::UpdateState2D;
        } else {
            m_lineSolverNorth  = std::make_unique<LineSolver<TC::n, Wstag, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            m_lineSolverSouth  = std::make_unique<LineSolver<TC::s, Wstag, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            SolutionUpdater    = &PlaneSolver::Sweep3D;
            StateUpdater       = &PlaneSolver::UpdateState3D;
        }
    }


    void SolvePlane(const intType k)
    { (this->*SolutionUpdater)(k); }


    void UpdateState()
    { (this->*StateUpdater)(); }


private:

    FieldData<Tensor3D> &m_fields;
    const FVCoefficients &m_fvCoeffs;

    std::unique_ptr< LineSolver<TC::n, Wstag, MI > > m_lineSolverNorth;
    std::unique_ptr< LineSolver<TC::s, Wstag, MI > > m_lineSolverSouth;
    std::unique_ptr< LineSolver<TC::p, Wstag, MI > > m_lineSolverCenter;

    FieldData<Tensor2D> m_planeConstants;

    void (PlaneSolver::*SolutionUpdater)(const intType);
    void (PlaneSolver::*StateUpdater)();

    intType m_ni, m_nj;

    // For 3D simulations
    void Sweep3D(intType k)
    {
        UpdatePlaneConstants(k);

        for (intType j = 0; j != m_nj - 1; j++) { // Forward sweep
            m_lineSolverNorth->SolveLine(j, k, m_planeConstants);
        }

        for (intType j = m_nj - 1; j != 0; j--) { // Backward sweep
            m_lineSolverSouth->SolveLine(j, k, m_planeConstants);
        }
    }

    void UpdateState3D()
    {
        m_lineSolverNorth->UpdateState();
        m_lineSolverSouth->UpdateState();
    }


    // For 2D simulations
    void Sweep2D(intType k)
    { 
        UpdatePlaneConstants(k);
        m_lineSolverCenter->SolveLine(0, k, m_planeConstants); 
    }

    void UpdateState2D()
    { m_lineSolverCenter->UpdateState(); }


    // Precalculate parts of stencil that are constant along a plane
    void UpdatePlaneConstants(intType k)
    {
        m_planeConstants = CalculatePlaneConstants<Wstag, MI >(k, m_fvCoeffs, m_fields);
    }
};

}   // end namespace CFD    


#endif // PLANE_SOLVER