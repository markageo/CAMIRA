#ifndef LINE_SOLVER
#define LINE_SOLVER

#include "StaggerIndexing.h"
#include "TriadSolver.h"

#include "../Types.h"
#include "../Macros.h"
#include "../IO/InputProcessing.h"
#include "../Tools/SweepTransformations.h"
#include "../Tools/FVTools.h"
#include "../Tools/FVLookups.h"
#include "../FiniteVolume/FiniteVolume.h"

namespace CFD
{

using namespace FVT;

template < TransportCoefficients::ENUMDATA Vstag,
           TransportCoefficients::ENUMDATA Wstag,
           MomentumInterpolation MI,
           Linearisation LI >
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
                const FVCoefficients &fvCoeffs) : 
                    m_fields( fields ),
                    m_fvCoeffs( fvCoeffs ),
                    m_lineConstants( Tensor1D( fvCoeffs.nCells(Axis::X) ) ),
                    m_ni( fvCoeffs.nCells(Axis::X) )
    {
        if (m_ni == 1) {
            m_triadSolverCenter = std::make_unique<TriadSolver<TC::p, Vstag, Wstag, MI, LI>>(fields, fieldsOld, mask, fvCoeffs);
            SolutionUpdater = &LineSolver::Sweep2D;
            StateUpdater = &LineSolver::UpdateState2D;
        } else {
            m_triadSolverEast = std::make_unique<TriadSolver<TC::e, Vstag, Wstag, MI, LI>>(fields, fieldsOld, mask, fvCoeffs);
            m_triadSolverWest = std::make_unique<TriadSolver<TC::w, Vstag, Wstag, MI, LI>>(fields, fieldsOld, mask, fvCoeffs);
            SolutionUpdater = &LineSolver::Sweep3D;
            StateUpdater = &LineSolver::UpdateState3D;
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

    std::unique_ptr<TriadSolver<TC::e, Vstag, Wstag, MI, LI>> m_triadSolverEast;
    std::unique_ptr<TriadSolver<TC::w, Vstag, Wstag, MI, LI>> m_triadSolverWest;
    std::unique_ptr<TriadSolver<TC::p, Vstag, Wstag, MI, LI>> m_triadSolverCenter;

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
        }

        for (intType i = m_ni - 1; i != 0; i--) { // Backward sweep
            m_triadSolverWest->UpdateTriad(i, j, k, m_lineConstants);
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
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        using sCV = typename StaggerIndexing< Axis::Y, Vstag >::ContinuityVelocity;
        using sVP = typename StaggerIndexing< Axis::Y, Vstag >::MomentumPressure;
        using sCW = typename StaggerIndexing< Axis::Z, Wstag >::ContinuityVelocity;

        // Staggered indices, U momenutm is not staggered wrt to the line
        intType jV{ j + sCV::iCoupled }, kV{ k                 }; // V momentum
        intType jW{ j                 }, kW{ k + sCW::iCoupled }; // W momentum

        // Ghost cells
        intType jgV{ G(jV) }, kgV{ G(kV) };
        intType jgW{ G(jW) }, kgW{ G(kW) };
        intType   jg{ G(j) },   kg{ G(k) };

        CFD_PRAGMA_VECTORIZE
        for ( intType i = 0; i != m_ni; i++ ) {
            intType ig{ G(i) };

            // U momentum
            floatType newtonStencilX = 0.0f;
            if constexpr ( LI == Linearisation::Newton ) {
                newtonStencilX = - m_fvCoeffs.Mom[X].AU[Y][sCV::cLeft ]( i, j, k ) * m_fields.U[Y]( ig, jg+sCV::iLeft , kg )
                                 - m_fvCoeffs.Mom[X].AU[Y][sCV::cRight]( i, j, k ) * m_fields.U[Y]( ig, jg+sCV::iRight, kg );
            }
            m_lineConstants.U[X](i) = planeConstants.U[X](i, j)
                                    + ( 
                                      - m_fvCoeffs.Mom[X].AU[X][n](i, j, k) * m_fields.U[X]( ig  , jg+1, kg  )
                                      - m_fvCoeffs.Mom[X].AU[X][s](i, j, k) * m_fields.U[X]( ig  , jg-1, kg  )

                                      + newtonStencilX
                                      );

            // V momentum
            floatType newtonStencilY = 0.0f;
            if constexpr ( LI == Linearisation::Newton ) {
                newtonStencilY = - m_fvCoeffs.Mom[Y].AU[X][e]( i, jV, kV ) * m_fields.U[X]( ig+1, jgV  , kgV  )
                                 - m_fvCoeffs.Mom[Y].AU[X][p]( i, jV, kV ) * m_fields.U[X]( ig  , jgV  , kgV  )
                                 - m_fvCoeffs.Mom[Y].AU[X][w]( i, jV, kV ) * m_fields.U[X]( ig-1, jgV  , kgV  )

                                 - m_fvCoeffs.Mom[Y].AU[Z][sCW::cCoupled]( i, jV, kV ) * m_fields.U[Z]( ig, jgV, kgV+sCW::iCoupled );
            }
            m_lineConstants.U[Y](i) = planeConstants.U[Y](i, jV)
                                    + (
                                      - m_fvCoeffs.Mom[Y].AU[Y][n](i, jV, kV) * m_fields.U[Y]( ig  , jgV+1, kgV  )  
                                      - m_fvCoeffs.Mom[Y].AU[Y][s](i, jV, kV) * m_fields.U[Y]( ig  , jgV-1, kgV  ) 

                                      - m_fvCoeffs.Mom[Y].AP[sVP::cLeft ](jV) * m_fields.P( ig, jgV + sVP::iLeft , kgV)
                                      - m_fvCoeffs.Mom[Y].AP[sVP::cRight](jV) * m_fields.P( ig, jgV + sVP::iRight, kgV)

                                      + newtonStencilY
                                      );

            // W momentum
            floatType newtonStencilZ = 0.0f;
            m_lineConstants.U[Z](i) = planeConstants.U[Z](i, jW)
                                    + ( 
                                      - m_fvCoeffs.Mom[Z].AU[Z][n](i, jW, kW) * m_fields.U[Z]( ig  , jgW+1, kgW  ) 
                                      - m_fvCoeffs.Mom[Z].AU[Z][s](i, jW, kW) * m_fields.U[Z]( ig  , jgW-1, kgW  ) 
                                      
                                      + newtonStencilZ
                                      );

            // Continuity equation
            floatType pressureWideStencil = 0.0f;
            if constexpr ( MI == MomentumInterpolation::Implicit ) {
                pressureWideStencil = - m_fvCoeffs.Cont.AP[nn](i, j, k) * m_fields.P( ig  , jg+2, kg  )
                                      - m_fvCoeffs.Cont.AP[ss](i, j, k) * m_fields.P( ig  , jg-2, kg  );
            }

            m_lineConstants.P(i) = planeConstants.P(i, j)

                                 - m_fvCoeffs.Cont.AU[Y][sCV::cLeft ](j) * m_fields.U[Y]( ig, jg + sCV::iLeft , kg)
                                 - m_fvCoeffs.Cont.AU[Y][sCV::cRight](j) * m_fields.U[Y]( ig, jg + sCV::iRight, kg)

                                 - m_fvCoeffs.Cont.AP[n](i, j, k) * m_fields.P( ig  , jg+1, kg  ) 
                                 - m_fvCoeffs.Cont.AP[s](i, j, k) * m_fields.P( ig  , jg-1, kg  )
                                 
                                 + pressureWideStencil;

        }

    }

};

}   // end namespace CFD    


#endif // LINE_SOLVER