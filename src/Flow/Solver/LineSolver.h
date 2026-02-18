#ifndef CAMIRA_LINE_SOLVER
#define CAMIRA_LINE_SOLVER

#include "StaggerIndexing.h"
#include "StencilConstants.h"
#include "TriadSolver.h"

#include "Core/Types.h"
#include "Core/Macros.h"
#include "Core/FVTools.h"
#include "Core/FVLookups.h"
#include "Flow/InputProcessing/InputProcessing.h"
#include "Flow/CoordinateTransformations/AxisTransformationFunctions.h"
#include "Flow/FiniteVolume/FiniteVolume.h"

namespace CAMIRA
{

using namespace FVT;

template < TransportCoefficients::ENUMDATA Vstag,
           TransportCoefficients::ENUMDATA Wstag,
           MomentumInterpolation MI >
class LineSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    // Staggering must be valid
    static_assert( (Vstag == TC::n) || (Vstag == TC::s) || (Vstag == TC::p), "Invalid V momentum staggering" );
    static_assert( (Wstag == TC::t) || (Wstag == TC::b) || (Wstag == TC::p), "Invalid W momentum staggering" );

public:
    LineSolver( FieldData<Tensor3D> &fields,
                const FieldData<Tensor3D> &fieldsOld,
                const Tensor3D &mask,
                const FVCoefficients &fvCoeffs, 
                const InputData::SmootherSettings &smootherSettings) : 
                    m_fields( fields ),
                    m_fvCoeffs( fvCoeffs ),
                    m_lineConstants( Tensor1D( fvCoeffs.nCells(Axis::X) + 2*nGhost ).setZero() ),
                    m_ni( fvCoeffs.nCells(Axis::X) )
    {
        if (m_ni == 1) {
            m_triadSolverCenter = std::make_unique<TriadSolver<TC::p, Vstag, Wstag, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            SolutionUpdater     = &LineSolver::Sweep2D;
            StateUpdater        = &LineSolver::UpdateState2D;
        } else {
            m_triadSolverEast   = std::make_unique<TriadSolver<TC::e, Vstag, Wstag, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            m_triadSolverWest   = std::make_unique<TriadSolver<TC::w, Vstag, Wstag, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            SolutionUpdater     = &LineSolver::Sweep3D;
            StateUpdater        = &LineSolver::UpdateState3D;
        }
    }

    void SolveLine( const intType j, 
                    const intType k, 
                    const FieldData<Tensor2D> &planeConstants)
    { (this->*SolutionUpdater)(j, k, planeConstants); }

    void UpdateState()
    { (this->*StateUpdater)(); }


private:
    FieldData<Tensor3D> &m_fields;
    const FVCoefficients &m_fvCoeffs;

    std::unique_ptr<TriadSolver<TC::e, Vstag, Wstag, MI >> m_triadSolverEast;
    std::unique_ptr<TriadSolver<TC::w, Vstag, Wstag, MI >> m_triadSolverWest;
    std::unique_ptr<TriadSolver<TC::p, Vstag, Wstag, MI >> m_triadSolverCenter;

    FieldData<Tensor1D> m_lineConstants;

    void (LineSolver::*SolutionUpdater)(const intType, const intType, const FieldData<Tensor2D> &);
    void (LineSolver::*StateUpdater)(void);

    intType m_ni;

    // For 3D simulations
    void Sweep3D( const intType j, 
                  const intType k, 
                  const FieldData<Tensor2D> &planeConstants)
    {
        UpdateLineConstants(j, k, planeConstants);

        for (intType i = 0; i != m_ni - 1; i++) { // Forward sweep
            m_triadSolverEast->UpdateTriad(i, j, k, m_lineConstants);
            // m_triadSolverEast->UpdateTriad(i, j, k);
        }

        for (intType i = m_ni - 1; i != 0; i--) { // Backward sweep
            m_triadSolverWest->UpdateTriad(i, j, k, m_lineConstants);
            // m_triadSolverWest->UpdateTriad(i, j, k);
        }
    }

    void UpdateState3D()
    {
        m_triadSolverEast->UpdateGlobalConstants();
        m_triadSolverWest->UpdateGlobalConstants();
    }

    // For 2D simulations
    void Sweep2D( const intType j, 
                  const intType k,
                  const FieldData<Tensor2D> &planeConstants )
    { 
        UpdateLineConstants(j, k, planeConstants);
        m_triadSolverCenter->UpdateTriad(0, j, k, m_lineConstants); 
    }

    void UpdateState2D()
    { m_triadSolverCenter->UpdateGlobalConstants(); }


    // Precalculate parts of stencil that are constant along a line
    void UpdateLineConstants( const intType j, 
                              const intType k, 
                              const FieldData<Tensor2D> &planeConstants )
    {
        m_lineConstants = CalculateLineConstants<Vstag, Wstag, MI >(j, k, planeConstants, m_fvCoeffs, m_fields);
    }

};

}   // end namespace CAMIRA    


#endif // LINE_SOLVER