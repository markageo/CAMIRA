#ifndef PLANE_SOLVER
#define PLANE_SOLVER

#include "StaggerIndexing.h"
#include "LineSolver.h"

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

template <TransportCoefficients::ENUMDATA Wstag, 
          MomentumInterpolation MI>
class PlaneSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    // Staggering must be valid
    static_assert( (Wstag == TC::t) || (Wstag == TC::b) || (Wstag == TC::p), "Invalid W momentum staggering" );

public:
    PlaneSolver( FieldData<array3D> &fields,
                 const FieldData<array3D> &fieldsOld,
                 const FVCoefficients &fvCoeffs) : 
                    m_fields( fields ),
                    m_fvCoeffs( fvCoeffs ),
                    m_planeConstants( array2D( fvCoeffs.nCells(Axis::X), fvCoeffs.nCells(Axis::Y) ) ),
                    m_ni( fvCoeffs.nCells(Axis::X) ),
                    m_nj( fvCoeffs.nCells(Axis::Y) )
    {
        if (m_nj == 1) {
            m_lineSolverCenter = std::make_unique<LineSolver<TC::p, Wstag, MI>>(fields, fieldsOld, fvCoeffs);
            SolutionUpdater = &PlaneSolver::Sweep2D;
            StateUpdater = &PlaneSolver::UpdateState2D;
        } else {
            m_lineSolverNorth = std::make_unique<LineSolver<TC::n, Wstag, MI>>(fields, fieldsOld, fvCoeffs);
            m_lineSolverSouth = std::make_unique<LineSolver<TC::s, Wstag, MI>>(fields, fieldsOld, fvCoeffs);
            SolutionUpdater = &PlaneSolver::Sweep3D;
            StateUpdater = &PlaneSolver::UpdateState3D;
        }
    }


    void SolvePlane(const intType k)
    { (this->*SolutionUpdater)(k); }


    void UpdateState()
    { (this->*StateUpdater)(); }


private:

    FieldData<array3D> &m_fields;
    const FVCoefficients &m_fvCoeffs;

    std::unique_ptr< LineSolver<TC::n, Wstag, MI> > m_lineSolverNorth;
    std::unique_ptr< LineSolver<TC::s, Wstag, MI> > m_lineSolverSouth;
    std::unique_ptr< LineSolver<TC::p, Wstag, MI> > m_lineSolverCenter;

    FieldData<array2D> m_planeConstants;

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

        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        using sCW = typename StaggerIndexing< Axis::Z, Wstag >::ContinuityVelocity;
        using sWP = typename StaggerIndexing< Axis::Z, Wstag >::MomentumPressure;

        // Staggered indices, U and V momentum is not staggered with respect to a plane
        intType kW{ k + sCW::iCoupled }; // W momentum

        // Ghost cells
        intType kgW{ G(kW) };
        intType  kg{ G(k) };

        for ( intType j = 0; j != m_nj; j++ ) {

            CFD_PRAGMA_VECTORIZE
            for ( intType i = 0; i != m_ni; i++ ) {
                intType ig{ G(i) }, jg{ G(j) };

                // U momentum
                m_planeConstants.U[X](i, j) = ( m_fvCoeffs.Mom[X].B(i, j, k)

                                              - m_fvCoeffs.Mom[X].AU[X][t](i, j, k) * m_fields.U[X]( ig  , jg  , kg+1) 
                                              - m_fvCoeffs.Mom[X].AU[X][b](i, j, k) * m_fields.U[X]( ig  , jg  , kg-1)

                                              );

                // V momentum
                m_planeConstants.U[Y](i, j) = ( m_fvCoeffs.Mom[Y].B(i, j, k)

                                              - m_fvCoeffs.Mom[Y].AU[Y][t](i, j, k) * m_fields.U[Y]( ig  , jg  , kg+1) 
                                              - m_fvCoeffs.Mom[Y].AU[Y][b](i, j, k) * m_fields.U[Y]( ig  , jg  , kg-1)

                                              );

                // W momentum 
                m_planeConstants.U[Z](i, j) = ( m_fvCoeffs.Mom[Z].B(i, j, kW)
                                
                                            - m_fvCoeffs.Mom[Z].AU[Z][t](i, j, kW) * m_fields.U[Z]( ig  , jg  , kgW+1) 
                                            - m_fvCoeffs.Mom[Z].AU[Z][b](i, j, kW) * m_fields.U[Z]( ig  , jg  , kgW-1)

                                            - m_fvCoeffs.Mom[Z].AP[sWP::cLeft ](kW) * m_fields.P( ig, jg, kgW + sWP::iLeft ) 
                                            - m_fvCoeffs.Mom[Z].AP[sWP::cRight](kW) * m_fields.P( ig, jg, kgW + sWP::iRight)

                                            );


                // Continuity equation
                floatType pressureWideStencil = 0.0f;
                if constexpr ( MI == MomentumInterpolation::Implicit ) {
                    pressureWideStencil = - m_fvCoeffs.Cont.AP[tt](i, j, k) * m_fields.P( ig  , jg  , kg+2) 
                                          - m_fvCoeffs.Cont.AP[bb](i, j, k) * m_fields.P( ig  , jg  , kg-2);
                }

                m_planeConstants.P(i, j) = m_fvCoeffs.Cont.B(i, j, k)

                                         - m_fvCoeffs.Cont.AU[Z][sCW::cLeft ](k) * m_fields.U[Z]( ig, jg, kg + sCW::iLeft )
                                         - m_fvCoeffs.Cont.AU[Z][sCW::cRight](k) * m_fields.U[Z]( ig, jg, kg + sCW::iRight)

                                         - m_fvCoeffs.Cont.AP[t](i, j, k) * m_fields.P( ig  , jg  , kg+1) 
                                         - m_fvCoeffs.Cont.AP[b](i, j, k) * m_fields.P( ig  , jg  , kg-1)
        
                                         + pressureWideStencil;
                
            }
        }

    }
};

}   // end namespace CFD    


#endif // PLANE_SOLVER