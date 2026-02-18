#ifndef CAMIRA_LINEAR_SOLVER
#define CAMIRA_LINEAR_SOLVER

#include "PlaneSolver.h"
#include "TriadSolver.h"
#include "StencilConstants.h"
#include "ResidualFunctions.h"

#include "Core/Types.h"
#include "Core/FVTools.h"
#include "Core/FVLookups.h"
#include "Core/ArrayIndexConversions.h"
#include "Flow/InputProcessing/InputProcessing.h"
#include "Flow/FiniteVolume/FiniteVolume.h"
#include "Flow/FiniteVolume/GhostCells.h"
#include "Flow/Parallel/Parallel.h"

#include <omp.h>

#include <cassert>

namespace CAMIRA
{

using namespace FVT;

// Linear solver interface class
template< MomentumInterpolation MI >
class LinearSolverInterface
{
    public:
        virtual void Solve() = 0;
        virtual void UpdateState() = 0;
        virtual ~LinearSolverInterface() {};
};


// Two orientation symmetric sweeping (Serial)
template< MomentumInterpolation MI >
class domainSymmetricSolverSerial : public LinearSolverInterface< MI >
{
    using TC = TransportCoefficients::ENUMDATA;
    using A = Axis::ENUMDATA;

public:
    domainSymmetricSolverSerial( FieldData<Tensor3D> &fields,
                                 const FieldData<Tensor3D> &fieldsOld,
                                 const Tensor3D &mask,
                                 const FVCoefficients &fvCoeffs, 
                                 const Mesh &mesh,
                                 const BoundaryConditionData &bcData,
                                 const InputData::SmootherSettings &smootherSettings) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_fvCoeffs( fvCoeffs ),
                    m_mesh( mesh ),
                    m_bcData( bcData ),
                    m_maxIterations( smootherSettings.maxIterations ),
                    m_maxResiduals( smootherSettings.maxResiduals ),
                    m_relaxation( smootherSettings.relaxation ),

                    m_ni( fvCoeffs.nCells(A::X) ),
                    m_nj( fvCoeffs.nCells(A::Y) ),
                    m_nk( fvCoeffs.nCells(A::Z) )
    {
        m_triadSolverForward  = std::make_unique<TriadSolver<TC::e, TC::n, TC::t, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
        m_triadSolverBackward = std::make_unique<TriadSolver<TC::w, TC::s, TC::b, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
    }


    void Solve()
    {
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        for ( intType nIterations = 1 ; nIterations <= m_maxIterations; nIterations++ )
        {
            // Reset residuals
            ForAllFieldData( [&] (intType f) { m_residuals[f] = 1.0f; });

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
    const Mesh &m_mesh;
    const BoundaryConditionData &m_bcData;
    const intType m_maxIterations;
    const FieldData<floatType> m_maxResiduals;
    const FieldData<floatType> m_relaxation;

    std::unique_ptr<TriadSolver<TC::e, TC::n, TC::t, MI >> m_triadSolverForward;
    std::unique_ptr<TriadSolver<TC::w, TC::s, TC::b, MI >> m_triadSolverBackward;

    FieldData<floatType> m_residuals, m_residualsInitialInv;

    intType m_ni, m_nj, m_nk;

    void Sweep3D()
    {
        
        TIC("Setting Ghost Cells")
        SetGhostCells(m_fields, m_mesh, m_bcData);
        TOC()

        // Triad starting on lo side
        TIC("Sweeping")
        for ( intType k = 0; k != m_nk; k++ ) {

            // FieldData<Tensor2D> planeConstants = CalculatePlaneConstants<TC::t, MI >(k, m_fvCoeffs, m_fields);

            for ( intType j = 0; j != m_nj; j++ ) {

                // FieldData<Tensor1D> lineConstants = CalculateLineConstants<TC::n, TC::t, MI >(j, k, planeConstants, m_fvCoeffs, m_fields);

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
        TOC()

        TIC("Setting Ghost Cells")
        SetGhostCells(m_fields, m_mesh, m_bcData);
        TOC()


        // Triad starting on hi side
        TIC("Sweeping")
        for ( intType k = m_nk-1; k != -1; k-- ) {

            // FieldData<Tensor2D> planeConstants = CalculatePlaneConstants<TC::b, MI >(k, m_fvCoeffs, m_fields);

            for ( intType j = m_nj-1; j != -1; j-- ) {

                // FieldData<Tensor1D> lineConstants = CalculateLineConstants<TC::s, TC::b, MI >(j, k, planeConstants, m_fvCoeffs, m_fields);

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





// Nested symmetric sweeping (Serial)
template< MomentumInterpolation MI >
class nestedLineSymmetricSolverSerial : public LinearSolverInterface< MI >
{
    using TC = TransportCoefficients::ENUMDATA;
    using A = Axis::ENUMDATA;

public:
    nestedLineSymmetricSolverSerial( FieldData<Tensor3D> &fields,
                                     const FieldData<Tensor3D> &fieldsOld,
                                     const Tensor3D &mask,
                                     const FVCoefficients &fvCoeffs, 
                                     const Mesh &mesh,
                                     const BoundaryConditionData &bcData,
                                     const InputData::SmootherSettings &smootherSettings) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_mesh( mesh ),
                    m_bcData( bcData ),
                    m_maxIterations( smootherSettings.maxIterations ),
                    m_maxResiduals( smootherSettings.maxResiduals ),
                    m_oldPlane( Tensor2D( m_fields.P.dimension(A::X), m_fields.P.dimension(A::Y) ) ),
                    m_ni( fvCoeffs.nCells(A::X) ),
                    m_nj( fvCoeffs.nCells(A::Y) ),
                    m_nk( fvCoeffs.nCells(A::Z) )
    {
        if (m_nk == 1) {
            m_planeSolverCenter = std::make_unique<PlaneSolver<TC::p, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            SolutionUpdater     = &nestedLineSymmetricSolverSerial::Sweep2D;
            StateUpdater        = &nestedLineSymmetricSolverSerial::UpdateState2D;
        } else {
            m_planeSolverTop    = std::make_unique<PlaneSolver<TC::t, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            m_planeSolverBottom = std::make_unique<PlaneSolver<TC::b, MI >>(fields, fieldsOld, mask, fvCoeffs, smootherSettings);
            SolutionUpdater     = &nestedLineSymmetricSolverSerial::Sweep3D;
            StateUpdater        = &nestedLineSymmetricSolverSerial::UpdateState3D;
        }
    }


    void Solve()
    {
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        for ( intType nIterations = 1 ; nIterations <= m_maxIterations; nIterations++ )
        {
            // Reset residuals
            ForAllFieldData( [&] (intType f) { m_residuals[f] = 1.0f; });

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
    const Mesh &m_mesh;
    const BoundaryConditionData &m_bcData;
    const intType m_maxIterations;
    const FieldData<floatType> m_maxResiduals;

    std::unique_ptr< PlaneSolver<TC::t, MI > > m_planeSolverTop;
    std::unique_ptr< PlaneSolver<TC::b, MI > > m_planeSolverBottom;
    std::unique_ptr< PlaneSolver<TC::p, MI > > m_planeSolverCenter;

    void (nestedLineSymmetricSolverSerial::*SolutionUpdater)(void);
    void (nestedLineSymmetricSolverSerial::*StateUpdater)(void);

    FieldData<Tensor2D> m_oldPlane;
    FieldData<floatType> m_residuals, m_residualsInitialInv;
    FieldData<intType> m_kS;

    intType m_ni, m_nj, m_nk;

    // For 3D simulations
    void Sweep3D()
    {
        TIC("Setting Ghost Cells")
        SetGhostCells(m_fields, m_mesh, m_bcData);
        TOC()

        TIC("Sweeping")
        for (intType k = 0; k != m_nk - 1; k++) { // Forward sweep
            Update(m_planeSolverTop, k);
        }
        TOC()

        TIC("Setting Ghost Cells")
        SetGhostCells(m_fields, m_mesh, m_bcData);
        TOC()

        TIC("Sweeping")
        for (intType k = m_nk - 1; k != 0; k--) { // Backward sweep
            Update(m_planeSolverBottom, k);
        }
        TOC()

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
    void Update( std::unique_ptr<PlaneSolver<Wstag, MI >> &planeSolver, intType k )
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




// Two orientation symmetric sweeping (Parallel)
template< MomentumInterpolation MI >
class domainSymmetricSolverParallel : public LinearSolverInterface< MI >
{
    using TC = TransportCoefficients::ENUMDATA;
    using A = Axis::ENUMDATA;    

public:
    domainSymmetricSolverParallel( FieldData<Tensor3D> &fields,
                                   const FieldData<Tensor3D> &fieldsOld,
                                   const Tensor3D &mask,
                                   const FVCoefficients &fvCoeffs, 
                                   const Mesh &mesh,
                                   const BoundaryConditionData &bcData,
                                   const InputData::SmootherSettings &smootherSettings ) : 

                    m_fields( fields ),
                    m_pressureOld( fields.P ),
                    m_fvCoeffs( fvCoeffs ),
                    m_mesh( mesh ),
                    m_bcData( bcData ),
                    m_maxIterations( smootherSettings.maxIterations ),
                    m_maxResiduals( smootherSettings.maxResiduals ),
                    m_relaxation( smootherSettings.relaxation ),

                    #if   defined( CAMIRA_PARALLEL_COLOURING_POINTS )

                        m_forwardColorSet( CreateForward3DColorSet(m_fvCoeffs.nCells) ),
                        m_reverseColorSet( CreateReverse3DColorSet(m_fvCoeffs.nCells) ),

                    #elif defined( CAMIRA_PARALLEL_COLOURING_LINES )

                        m_forwardColorSet( CreateForward2DColorSet(m_fvCoeffs.nCells) ),
                        m_reverseColorSet( CreateReverse2DColorSet(m_fvCoeffs.nCells) ),

                    #elif defined( CAMIRA_PARALLEL_COLOURING_PLANES )

                        m_forwardColorSet( CreateForward1DColorSet(m_fvCoeffs.nCells(2)) ),
                        m_reverseColorSet( CreateReverse1DColorSet(m_fvCoeffs.nCells(2)) ),

                    #endif

                    m_triadSolverForward(fields, fieldsOld, mask, fvCoeffs, smootherSettings),
                    m_triadSolverBackward(fields, fieldsOld, mask, fvCoeffs, smootherSettings),

                    m_nCells( fvCoeffs.nCells )
    { }


    void Solve()
    {
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        for ( intType nIterations = 1 ; nIterations <= m_maxIterations; nIterations++ )
        {
            m_pressureOld = m_fields.P;     // Only look at the pressure field to save time

            // Sweep domain
            Sweep3D();

            // Normalise residuals
            m_residual = L1Diff(m_pressureOld, m_fields.P);
            if ( nIterations == 1 ) {
                m_residualInitialInv = 1.0f / m_residual;

            }
            m_residual *= m_residualInitialInv;

            if ( m_residual <= m_maxResiduals.P ) {
                break;
            }
        }

    }


    // Update any precomputed values
    void UpdateState()
    {
        m_triadSolverForward.UpdateGlobalConstants();
        m_triadSolverBackward.UpdateGlobalConstants();
    }


private:

    FieldData<Tensor3D> &m_fields;
    Tensor3D m_pressureOld;
    const FVCoefficients &m_fvCoeffs;
    const Mesh &m_mesh;
    const BoundaryConditionData &m_bcData;
    const intType m_maxIterations;
    const FieldData<floatType> m_maxResiduals;
    const FieldData<floatType> m_relaxation;

    std::vector< std::vector< intType > > m_forwardColorSet,
                                          m_reverseColorSet;

    TriadSolver<TC::e, TC::n, TC::t, MI > m_triadSolverForward;
    TriadSolver<TC::w, TC::s, TC::b, MI > m_triadSolverBackward;

    iArray3 m_nCells;
    floatType m_residual, m_residualInitialInv;

   
    __attribute__((flatten))
    void Sweep3D()
    {
        #if defined( CAMIRA_PARALLEL_COLOURING_LINES )
            const iArray2 nCellsPlane = { m_nCells(1), m_nCells(2) };
        #endif

        TIC("Setting Ghost Cells")
        SetGhostCells(m_fields, m_mesh, m_bcData);
        TOC()

        // Forward sweep
        TIC("Sweeping")
        #if   defined( CAMIRA_PARALLEL_COLOURING_POINTS )

            for ( const std::vector<intType> &color : m_forwardColorSet ) {

                #pragma omp parallel for
                for ( intType idx : color ) {

                    const auto subs = Ind2Sub(  m_nCells, idx );
                    const intType i = subs[0],
                                  j = subs[1],
                                  k = subs[2];

                    m_triadSolverForward.UpdateTriad( i, j, k );

                }
            }


        #elif defined( CAMIRA_PARALLEL_COLOURING_LINES )

            for ( const std::vector<intType> &color : m_forwardColorSet ) {

                #pragma omp parallel for
                for ( intType idx : color ) {

                    const auto subs = Ind2Sub( nCellsPlane, idx );
                    const intType j = subs[0],
                                  k = subs[1];

                    for ( intType i = 0; i != m_nCells(0); i++ ) {
                        m_triadSolverForward.UpdateTriad( i, j, k );
                    }

                }
            }


        #elif defined( CAMIRA_PARALLEL_COLOURING_PLANES )

            for ( const std::vector<intType> &color : m_forwardColorSet ) {

                #pragma omp parallel for
                for ( intType k : color ) {

                    for ( intType j = 0; j != m_nCells(1); j++ ) {
                        for ( intType i = 0; i != m_nCells(0); i++ ) {
                            m_triadSolverForward.UpdateTriad( i, j, k );
                        }
                    }

                }
            }

        #endif
        TOC()


        TIC("Setting Ghost Cells")
        SetGhostCells(m_fields, m_mesh, m_bcData);
        TOC()


        // Reverse plane sweep
        TIC("Sweeping")
        #if defined( CAMIRA_PARALLEL_COLOURING_POINTS )

            for ( const std::vector<intType> &color : m_reverseColorSet ) {

                #pragma omp parallel for
                for ( intType idx : color ) {

                    auto subs = Ind2Sub(  m_nCells, idx );
                    const intType i = subs[0],
                                  j = subs[1],
                                  k = subs[2];

                    m_triadSolverBackward.UpdateTriad( i, j, k );

                }
            }

        #elif defined( CAMIRA_PARALLEL_COLOURING_LINES )

            for ( const std::vector<intType> &color : m_reverseColorSet ) {

                #pragma omp parallel for
                for ( intType idx : color ) {

                    const auto subs = Ind2Sub( nCellsPlane, idx );
                    const intType j = subs[0],
                                  k = subs[1];

                    for ( intType i = m_nCells(0)-1; i != -1; i-- ) {
                        m_triadSolverBackward.UpdateTriad( i, j, k );
                    }

                }
            }

        #elif defined( CAMIRA_PARALLEL_COLOURING_PLANES )

            for ( const std::vector<intType> &color : m_reverseColorSet ) {

                #pragma omp parallel for
                for ( intType k : color ) {

                    for ( intType j = m_nCells(1)-1; j != -1; j-- ) {
                        for ( intType i = m_nCells(0)-1; i != -1; i-- ) {
                            m_triadSolverBackward.UpdateTriad( i, j, k );
                        }
                    }

                }
            } 

        #endif

        TOC()

    }

};


}   // end namespace CAMIRA    


#endif // LINEAR_SOLVER