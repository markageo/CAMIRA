#ifndef LINEAR_SOLVER
#define LINEAR_SOLVER

#include "PlaneSolver.h"
#include "TriadSolver.h"
#include "StencilConstants.h"
#include "ResidualFunctions.h"

#include "../Core/Types.h"
#include "../Core/FVTools.h"
#include "../Core/FVLookups.h"
#include "../IO/InputProcessing.h"
#include "../FiniteVolume/FiniteVolume.h"

namespace CFD
{

using namespace FVT;

// Linear solver interface class
template< MomentumInterpolation MI,
          Linearisation LI >
class LinearSolverInterface
{
    public:
        virtual void Solve() = 0;
        virtual void UpdateState() = 0;
};


// Two orientation symmetric sweeping
template< MomentumInterpolation MI,
          Linearisation LI >
class domainSymmetricSolver : public LinearSolverInterface< MI, LI >
{
    using TC = TransportCoefficients::ENUMDATA;
    using A = Axis::ENUMDATA;

public:
    domainSymmetricSolver( FieldData<Tensor3D> &fields,
                           const FieldData<Tensor3D> &fieldsOld,
                           const Tensor3D &mask,
                           const FVCoefficients &fvCoeffs, 
                           const InputData::SmootherSettings &smootherSettings) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_fvCoeffs( fvCoeffs ),
                    m_maxIterations( smootherSettings.maxIterations ),
                    m_maxResiduals( smootherSettings.maxResiduals ),
                    m_relaxation( smootherSettings.relaxation ),

                    m_ni( fvCoeffs.nCells(A::X) ),
                    m_nj( fvCoeffs.nCells(A::Y) ),
                    m_nk( fvCoeffs.nCells(A::Z) )
    {
        m_triadSolverForward  = std::make_unique<TriadSolver<TC::e, TC::n, TC::t, MI, LI>>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
        m_triadSolverBackward = std::make_unique<TriadSolver<TC::w, TC::s, TC::b, MI, LI>>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
    }


    void Solve()
    {
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        for ( intType nIterations = 1 ; nIterations <= m_maxIterations; nIterations++ )
        {
            // Reset residuals
            ForAllFieldData( [&] (intType f) { m_residuals[f] = 0.0f; });

            // Sweep domain
            Sweep3D();

            // Normalise residuals
            ForAllFieldData( [&] (intType f) { m_residuals[f] /= static_cast<floatType>(m_ni * m_nj * m_nk); });
            if ( nIterations == 1 ) {
                ForAllFieldData( [&] (intType f) { m_residualsInitialInv[f] = 1.0f / m_residuals[f]; });
            }
            NormaliseResiduals(m_residuals, m_residualsInitialInv);

            // Check residual tolerence
            if ( MetResidualTolerence(m_residuals, m_maxResiduals) ) {
                break;
            }
        }

    }


    // Update any precomputed values
    void UpdateState()
    {
        m_triadSolverForward->UpdateGlobalConstants();
        m_triadSolverBackward->UpdateGlobalConstants();
    }


private:

    FieldData<Tensor3D> &m_fields;
    const FieldData<Tensor3D> &m_fieldsOld;
    const FVCoefficients &m_fvCoeffs;
    const intType m_maxIterations;
    const FieldData<floatType> m_maxResiduals;
    const FieldData<floatType> m_relaxation;

    std::unique_ptr<TriadSolver<TC::e, TC::n, TC::t, MI, LI>> m_triadSolverForward;
    std::unique_ptr<TriadSolver<TC::w, TC::s, TC::b, MI, LI>> m_triadSolverBackward;

    FieldData<floatType> m_residuals, m_residualsInitialInv;

    intType m_ni, m_nj, m_nk;

    void Sweep3D()
    {

        TIC("Sweeping")
        // Triad starting on lo side
        for ( intType k = 0; k != m_nk; k++ ) {

            // FieldData<Tensor2D> planeConstants = CalculatePlaneConstants<TC::t, MI, LI>(k, m_fvCoeffs, m_fields);

            for ( intType j = 0; j != m_nj; j++ ) {

                // FieldData<Tensor1D> lineConstants = CalculateLineConstants<TC::n, TC::t, MI, LI>(j, k, planeConstants, m_fvCoeffs, m_fields);

                for ( intType i = 0; i != m_ni; i++ ) {

                    FieldData<floatType> oldValues;
                    oldValues.P    = m_fields.P( G(i, j, k) );
                    oldValues.U[0] = m_fields.U[0]( G(i+1, j  , k  ) );
                    oldValues.U[1] = m_fields.U[1]( G(i  , j+1, k  ) );
                    oldValues.U[2] = m_fields.U[2]( G(i  , j  , k+1) );

                    // m_triadSolverForward->UpdateTriad( i, j, k, lineConstants );
                    m_triadSolverForward->UpdateTriad( i, j, k );

                    m_residuals.P    += abs( oldValues.P    - m_fields.P( G(i, j, k) ) );
                    m_residuals.U[0] += abs( oldValues.U[0] - m_fields.U[0]( G(i+1, j  , k  ) ) );
                    m_residuals.U[1] += abs( oldValues.U[1] - m_fields.U[1]( G(i  , j+1, k  ) ) );
                    m_residuals.U[2] += abs( oldValues.U[2] - m_fields.U[2]( G(i  , j  , k+1) ) );

                }
            }
        }


        // Triad starting on hi side
        for ( intType k = m_nk-1; k != -1; k-- ) {

            // FieldData<Tensor2D> planeConstants = CalculatePlaneConstants<TC::b, MI, LI>(k, m_fvCoeffs, m_fields);

            for ( intType j = m_nj-1; j != -1; j-- ) {

                // FieldData<Tensor1D> lineConstants = CalculateLineConstants<TC::s, TC::b, MI, LI>(j, k, planeConstants, m_fvCoeffs, m_fields);

                for ( intType i = m_ni-1; i != -1; i-- ) {

                    FieldData<floatType> oldValues;
                    oldValues.P    = m_fields.P( G(i, j, k) );
                    oldValues.U[0] = m_fields.U[0]( G(i-1, j  , k  ) );
                    oldValues.U[1] = m_fields.U[1]( G(i  , j-1, k  ) );
                    oldValues.U[2] = m_fields.U[2]( G(i  , j  , k-1) );

                    // m_triadSolverBackward->UpdateTriad( i, j, k, lineConstants );
                    m_triadSolverBackward->UpdateTriad( i, j, k );

                    m_residuals.P    += abs( oldValues.P    - m_fields.P( G(i, j, k) ) );
                    m_residuals.U[0] += abs( oldValues.U[0] - m_fields.U[0]( G(i-1, j  , k  ) ) );
                    m_residuals.U[1] += abs( oldValues.U[1] - m_fields.U[1]( G(i  , j-1, k  ) ) );
                    m_residuals.U[2] += abs( oldValues.U[2] - m_fields.U[2]( G(i  , j  , k-1) ) );

                }
            }
        }

        TOC()
    }


};





// Nested symmetric sweeping
template< MomentumInterpolation MI,
          Linearisation LI >
class nestedLineSymmetricSolver : public LinearSolverInterface< MI, LI >
{
    using TC = TransportCoefficients::ENUMDATA;
    using A = Axis::ENUMDATA;

public:
    nestedLineSymmetricSolver( FieldData<Tensor3D> &fields,
                               const FieldData<Tensor3D> &fieldsOld,
                               const Tensor3D &mask,
                               const FVCoefficients &fvCoeffs, 
                               const InputData::SmootherSettings &smootherSettings) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_maxIterations( smootherSettings.maxIterations ),
                    m_maxResiduals( smootherSettings.maxResiduals ),
                    m_oldPlane( Tensor2D( m_fields.P.dimension(A::X), m_fields.P.dimension(A::Y) ) ),
                    m_ni( fvCoeffs.nCells(A::X) ),
                    m_nj( fvCoeffs.nCells(A::Y) ),
                    m_nk( fvCoeffs.nCells(A::Z) )
    {
        if (m_nk == 1) {
            m_planeSolverCenter = std::make_unique<PlaneSolver<TC::p, MI, LI>>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            SolutionUpdater     = &nestedLineSymmetricSolver::Sweep2D;
            StateUpdater        = &nestedLineSymmetricSolver::UpdateState2D;
        } else {
            m_planeSolverTop    = std::make_unique<PlaneSolver<TC::t, MI, LI>>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            m_planeSolverBottom = std::make_unique<PlaneSolver<TC::b, MI, LI>>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            SolutionUpdater     = &nestedLineSymmetricSolver::Sweep3D;
            StateUpdater        = &nestedLineSymmetricSolver::UpdateState3D;
        }
    }


    void Solve()
    {
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        for ( intType nIterations = 1 ; nIterations <= m_maxIterations; nIterations++ )
        {
            // Reset residuals
            ForAllFieldData( [&] (intType f) { m_residuals[f] = 0.0f; });

            // Update plane
            (this->*SolutionUpdater)();

            // Normalise residuals
            ForAllFieldData( [&] (intType f) { m_residuals[f] /= static_cast<floatType>(m_ni * m_nj * m_nk); });
            if ( nIterations == 1 ) {
                ForAllFieldData( [&] (intType f) { m_residualsInitialInv[f] = 1.0f / m_residuals[f]; });
            }
            NormaliseResiduals(m_residuals, m_residualsInitialInv);

            // Check residual tolerence
            if ( MetResidualTolerence(m_residuals, m_maxResiduals) ) {
                break;
            }
        }

    }


    // Update any precomputed values
    void UpdateState()
    {
        (this->*StateUpdater)();
    }


private:

    FieldData<Tensor3D> &m_fields;
    const FieldData<Tensor3D> &m_fieldsOld;
    const intType m_maxIterations;
    const FieldData<floatType> m_maxResiduals;

    std::unique_ptr< PlaneSolver<TC::t, MI, LI> > m_planeSolverTop;
    std::unique_ptr< PlaneSolver<TC::b, MI, LI> > m_planeSolverBottom;
    std::unique_ptr< PlaneSolver<TC::p, MI, LI> > m_planeSolverCenter;

    void (nestedLineSymmetricSolver::*SolutionUpdater)(void);
    void (nestedLineSymmetricSolver::*StateUpdater)(void);

    FieldData<Tensor2D> m_oldPlane;
    FieldData<floatType> m_residuals, m_residualsInitialInv;
    FieldData<intType> m_kS;

    intType m_ni, m_nj, m_nk;

    // For 3D simulations
    void Sweep3D()
    {
        for (intType k = 0; k != m_nk - 1; k++) { // Forward sweep
            Update(m_planeSolverTop, k);
        }

        for (intType k = m_nk - 1; k != 0; k--) { // Backward sweep
            Update(m_planeSolverBottom, k);
        }
    }

    void UpdateState3D()
    {
        m_planeSolverTop->UpdateState();
        m_planeSolverBottom->UpdateState();
    }


    // For 2D simulations
    void Sweep2D()
    {
        Update(m_planeSolverCenter, 0);
    }

    void UpdateState2D()
    {
        m_planeSolverCenter->UpdateState();
    }


    template <TC Wstag>
    void Update( std::unique_ptr<PlaneSolver<Wstag, MI, LI>> &planeSolver, intType k )
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Axis::ENUMDATA;

        ForAllFieldData( [&] (intType f) { m_kS[f] = k; });  // Set iterating coefficient
        m_kS.U[Z] += LUT::CoeffIndex[Wstag];                 // W momentum is staggered

        ForAllFieldData( [&] (intType f) { m_oldPlane[f] = m_fields[f].chip( G(m_kS[f]), Z ); } ); // Set old plane

        planeSolver->SolvePlane(k);

        ForAllFieldData( [&] (intType f) {
            auto fieldPlane = m_fields[f].chip( G(m_kS[f]), Z );
            auto delta = fieldPlane - m_oldPlane[f];
            m_residuals[f] += static_cast<Tensor0D>( delta.abs().sum() )(0);
        } );
    }
};

}   // end namespace CFD    


#endif // LINEAR_SOLVER