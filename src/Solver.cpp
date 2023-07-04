#include "Solver.h"
#include "Types.h"
#include "InputProcessing.h"
#include "SweepTransformations.h"
#include "FiniteVolume.h"
#include "ConvergenceLogging.h"
#include "FieldProbe.h"
#include "StaggerIndexing.h"
#include "Utils.h"

#include <type_traits>
#include <iostream>
#include <memory>

namespace CFD
{

namespace
{

    // Calculated as the L1 norm of the difference between two arrays.
    [[ maybe_unused ]]
    FieldData<floatType> L1DiffResiduals( const FieldData<array3D> &fields1,
                                          const FieldData<array3D> &fields2 )
    {
        FieldData<floatType> result;
        ForAllFieldData( [&] (intType i) { 

            auto fieldDiff = fields2[i] - fields1[i];  // auto lazily evaluates
            result[i] = static_cast<array0D>( fieldDiff.abs().mean() )(0); 

        });

        return result;
    }



    // Calculate the absolute residual of each equation from the finite volume stencil
    template< MomentumInterpolation MI > [[ maybe_unused ]]
    FieldData<floatType> StencilResiduals( const FieldData<array3D> &fields,
                                           const FVCoefficients<MI> &fvCoeffs )
    {
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        FieldData<floatType> residuals{0};
        FieldData<floatType> scalingFactor{0};


        for ( intType k = 0; k != fvCoeffs.nCells[Z]; k++ ) {
            for ( intType j = 0; j != fvCoeffs.nCells[Y]; j++ ) {
                for ( intType i = 0; i != fvCoeffs.nCells[X]; i++ ) {

                    intType ig{ G(i) }, jg{ G(j) }, kg{ G(k) };

                    // U momentum
                    residuals.U[X] += abs( 
                                            fvCoeffs.Mom[X].AU[X][p](i, j, k) * fields.U[X]( ig  , jg  , kg  ) 
                                          + fvCoeffs.Mom[X].AU[X][n](i, j, k) * fields.U[X]( ig  , jg+1, kg  ) 
                                          + fvCoeffs.Mom[X].AU[X][e](i, j, k) * fields.U[X]( ig+1, jg  , kg  ) 
                                          + fvCoeffs.Mom[X].AU[X][s](i, j, k) * fields.U[X]( ig  , jg-1, kg  ) 
                                          + fvCoeffs.Mom[X].AU[X][w](i, j, k) * fields.U[X]( ig-1, jg  , kg  ) 
                                          + fvCoeffs.Mom[X].AU[X][t](i, j, k) * fields.U[X]( ig  , jg  , kg+1) 
                                          + fvCoeffs.Mom[X].AU[X][b](i, j, k) * fields.U[X]( ig  , jg  , kg-1) 

                                          + fvCoeffs.Mom[X].AP[e](i) * fields.P( ig+1, jg  , kg  )
                                          + fvCoeffs.Mom[X].AP[p](i) * fields.P( ig  , jg  , kg  )
                                          + fvCoeffs.Mom[X].AP[w](i) * fields.P( ig-1, jg  , kg  )

                                          - fvCoeffs.Mom[X].B(i, j, k)  );

                    scalingFactor.U[X] += abs( fvCoeffs.Mom[X].AU[X][p](i, j, k) * fields.U[X]( ig  , jg  , kg  ) );


                    // V momentum
                    residuals.U[Y] += abs( 
                                            fvCoeffs.Mom[Y].AU[Y][p](i, j, k) * fields.U[Y]( ig  , jg  , kg  ) 
                                          + fvCoeffs.Mom[Y].AU[Y][n](i, j, k) * fields.U[Y]( ig  , jg+1, kg  ) 
                                          + fvCoeffs.Mom[Y].AU[Y][e](i, j, k) * fields.U[Y]( ig+1, jg  , kg  ) 
                                          + fvCoeffs.Mom[Y].AU[Y][s](i, j, k) * fields.U[Y]( ig  , jg-1, kg  ) 
                                          + fvCoeffs.Mom[Y].AU[Y][w](i, j, k) * fields.U[Y]( ig-1, jg  , kg  ) 
                                          + fvCoeffs.Mom[Y].AU[Y][t](i, j, k) * fields.U[Y]( ig  , jg  , kg+1) 
                                          + fvCoeffs.Mom[Y].AU[Y][b](i, j, k) * fields.U[Y]( ig  , jg  , kg-1) 

                                          + fvCoeffs.Mom[Y].AP[n](j) * fields.P( ig  , jg+1, kg  )
                                          + fvCoeffs.Mom[Y].AP[p](j) * fields.P( ig  , jg  , kg  )
                                          + fvCoeffs.Mom[Y].AP[s](j) * fields.P( ig  , jg-1, kg  )

                                          - fvCoeffs.Mom[Y].B(i, j, k)  );

                    scalingFactor.U[Y] += abs( fvCoeffs.Mom[Y].AU[Y][p](i, j, k) * fields.U[Y]( ig  , jg  , kg  ) );


                    // W momentm
                    residuals.U[Z] += abs( 
                                            fvCoeffs.Mom[Z].AU[Z][p](i, j, k) * fields.U[Z]( ig  , jg  , kg  ) 
                                          + fvCoeffs.Mom[Z].AU[Z][n](i, j, k) * fields.U[Z]( ig  , jg+1, kg  ) 
                                          + fvCoeffs.Mom[Z].AU[Z][e](i, j, k) * fields.U[Z]( ig+1, jg  , kg  ) 
                                          + fvCoeffs.Mom[Z].AU[Z][s](i, j, k) * fields.U[Z]( ig  , jg-1, kg  ) 
                                          + fvCoeffs.Mom[Z].AU[Z][w](i, j, k) * fields.U[Z]( ig-1, jg  , kg  ) 
                                          + fvCoeffs.Mom[Z].AU[Z][t](i, j, k) * fields.U[Z]( ig  , jg  , kg+1) 
                                          + fvCoeffs.Mom[Z].AU[Z][b](i, j, k) * fields.U[Z]( ig  , jg  , kg-1) 

                                          + fvCoeffs.Mom[Z].AP[t](k) * fields.P( ig  , jg  , kg+1)
                                          + fvCoeffs.Mom[Z].AP[p](k) * fields.P( ig  , jg  , kg  )
                                          + fvCoeffs.Mom[Z].AP[b](k) * fields.P( ig  , jg  , kg-1)

                                          - fvCoeffs.Mom[Z].B(i, j, k)  );

                    scalingFactor.U[Z] += abs( fvCoeffs.Mom[Z].AU[Z][p](i, j, k) * fields.U[Z]( ig  , jg  , kg  ) );


                    // Continuity
                    residuals.P += abs( 
                                      fvCoeffs.Cont.AU[X][e](i) * fields.U[X]( ig+1, jg  , kg  )
                                    + fvCoeffs.Cont.AU[X][p](i) * fields.U[X]( ig  , jg  , kg  )
                                    + fvCoeffs.Cont.AU[X][w](i) * fields.U[X]( ig-1, jg  , kg  )

                                    + fvCoeffs.Cont.AU[Y][n](j) * fields.U[Y]( ig  , jg+1, kg  )
                                    + fvCoeffs.Cont.AU[Y][p](j) * fields.U[Y]( ig  , jg  , kg  )
                                    + fvCoeffs.Cont.AU[Y][s](j) * fields.U[Y]( ig  , jg-1, kg  )

                                    + fvCoeffs.Cont.AU[Z][t](k) * fields.U[Z]( ig  , jg  , kg+1)
                                    + fvCoeffs.Cont.AU[Z][p](k) * fields.U[Z]( ig  , jg  , kg  )
                                    + fvCoeffs.Cont.AU[Z][b](k) * fields.U[Z]( ig  , jg  , kg-1)

                                    + fvCoeffs.Cont.AP[p](i, j, k) * fields.P( ig  , jg  , kg  )
                                    + fvCoeffs.Cont.AP[n](i, j, k) * fields.P( ig  , jg+1, kg  ) 
                                    + fvCoeffs.Cont.AP[e](i, j, k) * fields.P( ig+1, jg  , kg  ) 
                                    + fvCoeffs.Cont.AP[s](i, j, k) * fields.P( ig  , jg-1, kg  ) 
                                    + fvCoeffs.Cont.AP[w](i, j, k) * fields.P( ig-1, jg  , kg  ) 
                                    + fvCoeffs.Cont.AP[t](i, j, k) * fields.P( ig  , jg  , kg+1) 
                                    + fvCoeffs.Cont.AP[b](i, j, k) * fields.P( ig  , jg  , kg-1)

                                    + fvCoeffs.Cont.AP[nn](i, j, k) * fields.P( ig  , jg+2, kg  ) 
                                    + fvCoeffs.Cont.AP[ee](i, j, k) * fields.P( ig+2, jg  , kg  ) 
                                    + fvCoeffs.Cont.AP[ss](i, j, k) * fields.P( ig  , jg-2, kg  ) 
                                    + fvCoeffs.Cont.AP[ww](i, j, k) * fields.P( ig-2, jg  , kg  ) 
                                    + fvCoeffs.Cont.AP[tt](i, j, k) * fields.P( ig  , jg  , kg+2) 
                                    + fvCoeffs.Cont.AP[bb](i, j, k) * fields.P( ig  , jg  , kg-2)
                                    
                                    - fvCoeffs.Cont.B(i, j, k) );

                }
            }
        }

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
            residuals.U[axis] /= scalingFactor.U[axis];
        } );

        // ForAllFieldData( [&] (intType f) {
        //     residuals[f] /= static_cast<floatType>( fvCoeffs.nCells[X] * fvCoeffs.nCells[Y] * fvCoeffs.nCells[Z] ); 
        // } );

        return residuals;
    }


    // Calculate global mass flux residual at the domain boundary
    [[ maybe_unused ]]
    floatType BoundaryMassFluxResidual( const EnumVector<Axis, array3D> &faceFluxes,
                                        const Mesh &mesh )
    {
        floatType massFluxResidual = 0.0f;

        EnumFor<Axis>([&](Axis::ENUMDATA axis) {

            // Positive face, area normal is in positive direction
            auto faceFluxesPositive = faceFluxes[axis].chip( mesh.nCells(axis), axis ) * mesh.cellFaceAreas[axis];
            massFluxResidual += static_cast<array0D>( faceFluxesPositive.sum() )(0);

            // Negative face, area normal is in negative direction
            auto faceFluxesNegative = - faceFluxes[axis].chip( 0, axis ) * mesh.cellFaceAreas[axis];
            massFluxResidual += static_cast<array0D>( faceFluxesNegative.sum() )(0);

        });

        return massFluxResidual;
    }



    // Normalise the residual by the first iteration
    [[ maybe_unused ]]
    void NormaliseResiduals( FieldData<floatType> &residuals,
                         FieldData<floatType> &residualsScaleFactor,
                         const intType nIterations )
    {
        if (nIterations == 1) {

            // Momentum equations
            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                residualsScaleFactor.U[axis] = 1.0f;
                // if ( residuals.U[axis] != 0 ) { // Division by zero
                //     residualsScaleFactor.U[axis] = 1.0f / residuals.U[axis];
                // } 

            } );

            // Continuity equation
            residualsScaleFactor.P = 1.0f;
            if ( residuals.P != 0 ) { // Division by zero
                residualsScaleFactor.P = 1.0f / residuals.P;
            } 
        }

        ForAllFieldData( [&] (intType i) { residuals[i] *= residualsScaleFactor[i]; } );
    }



    // Check if residual tolerence is met
    [[ maybe_unused ]]
    bool MetResidualTolerence( const FieldData<floatType> &residuals,
                               const FieldData<floatType> &residualsTarget )
    {
        for ( intType i = 0; i != FieldData<floatType>::nData; i++ ) {  // Can't use ForAllFieldData since returning inside loop

            if ( residuals[i] > residualsTarget[i] ) 
                return false;

        }
        return true;
    }



    // Update a vector of field probes
    std::vector< FieldData<floatType> > FieldProbeValues( const FieldData<array3D> &fields,
                                                          const std::vector<FieldProbe> &fieldProbes )
    {
        std::vector< FieldData<floatType> > probeValues( fieldProbes.size() );

        for ( size_t p = 0; p != fieldProbes.size(); p++ ) {
            ForAllFieldData( [&] (intType f) { 
                probeValues[p][f] =  fieldProbes[p].GetFieldValue( fields[f] );
            } );
        }

        return probeValues;
    }


} // end anonymous namespace





// --------------------------------------------------------- TriadSolver --------------------------------------------------------- //
template < TransportCoefficients::ENUMDATA Ustag,
           TransportCoefficients::ENUMDATA Vstag,
           TransportCoefficients::ENUMDATA Wstag, 
           MomentumInterpolation MI>
class TriadSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    // Staggering must be valid
    static_assert( (Ustag == TC::e) || (Ustag == TC::w) || (Ustag == TC::p), "Invalid U momentum staggering" );
    static_assert( (Vstag == TC::n) || (Vstag == TC::s) || (Vstag == TC::p), "Invalid V momentum staggering" );
    static_assert( (Wstag == TC::t) || (Wstag == TC::b) || (Wstag == TC::p), "Invalid W momentum staggering" );

    // Aliases for the staggering offsets
    using sCU = typename StaggerIndexing< Axis::X, Ustag >::ContinuityVelocity;
    using sUP = typename StaggerIndexing< Axis::X, Ustag >::MomentumPressure;

    using sCV = typename StaggerIndexing< Axis::Y, Vstag >::ContinuityVelocity;
    using sVP = typename StaggerIndexing< Axis::Y, Vstag >::MomentumPressure;

    using sCW = typename StaggerIndexing< Axis::Z, Wstag >::ContinuityVelocity;
    using sWP = typename StaggerIndexing< Axis::Z, Wstag >::MomentumPressure;

public:
    TriadSolver( FieldData<array3D> &fields,
                 const FieldData<array3D> &fieldsOld,
                 const FVCoefficients<MI> &fvCoeffs ) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_fvCoeffs( fvCoeffs ),
                    m_ni( fvCoeffs.nCells(0) ),
                    m_nj( fvCoeffs.nCells(1) ),
                    m_nk (fvCoeffs.nCells(2) ),
                    m_K( array3D(m_ni, m_nj, m_nk).setZero() )
    { UpdateGlobalConstants(); };



    // Core function which updates the local coupled system. Templated by staggering direction.
    // This makes use of precomputed line constants
    __attribute__((always_inline)) 
    inline void UpdateTriad( const intType i, 
                              const intType j, 
                              const intType k,
                              const FieldData<array1D> &lineConstants )
    {
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        // For indexing the staggered cells
        intType iU{ i + sCU::iCoupled }, jU{ j                 }, kU{ k                 }; // U momentum
        intType iV{ i                 }, jV{ j + sCV::iCoupled }, kV{ k                 }; // V momentum
        intType iW{ i                 }, jW{ j                 }, kW{ k + sCW::iCoupled }; // W momentum

        // With ghost cells, this is faster than using the G() function inline every time
        intType igU{ G(iU) }, jgU{ G(jU) }, kgU{ G(kU) };
        intType igV{ G(iV) }, jgV{ G(jV) }, kgV{ G(kV) };
        intType igW{ G(iW) }, jgW{ G(jW) }, kgW{ G(kW) };
        intType   ig{ G(i) },   jg{ G(j) },   kg{ G(k) };


        // Precompute momentum RHS divided by AP coefficients
        // U momentum
        floatType bU = ( lineConstants.U[X](iU)  

                       - m_fvCoeffs.Mom[X].AU[X][e](iU, jU, kU) * m_fields.U[X]( igU+1, jgU  , kgU  )
                       - m_fvCoeffs.Mom[X].AU[X][w](iU, jU, kU) * m_fields.U[X]( igU-1, jgU  , kgU  )

                       - m_fvCoeffs.Mom[X].AP[sUP::cLeft ](iU) * m_fields.P( igU + sUP::iLeft , jgU, kgU)
                       - m_fvCoeffs.Mom[X].AP[sUP::cRight](iU) * m_fields.P( igU + sUP::iRight, jgU, kgU) 

                       ) * m_fvCoeffs.Mom[X].diagCoeffInv(iU, jU, kU);


        // V momentum
        floatType bV = ( lineConstants.U[Y](iV)

                       - m_fvCoeffs.Mom[Y].AU[Y][e](iV, jV, kV) * m_fields.U[Y]( igV+1, jgV  , kgV  ) 
                       - m_fvCoeffs.Mom[Y].AU[Y][w](iV, jV, kV) * m_fields.U[Y]( igV-1, jgV  , kgV  ) 

                       ) * m_fvCoeffs.Mom[Y].diagCoeffInv(iV, jV, kV);


        // W momentum
        floatType bW = ( lineConstants.U[Z](iW)

                       - m_fvCoeffs.Mom[Z].AU[Z][e](iW, jW, kW) * m_fields.U[Z]( igW+1, jgW  , kgW  ) 
                       - m_fvCoeffs.Mom[Z].AU[Z][w](iW, jW, kW) * m_fields.U[Z]( igW-1, jgW  , kgW  ) 

                       ) * m_fvCoeffs.Mom[Z].diagCoeffInv(iW, jW, kW);


        // Continuity for pressure
        floatType bP = lineConstants.P(i)

                     - m_fvCoeffs.Cont.AU[X][sCU::cLeft ](i) * m_fields.U[X]( ig + sCU::iLeft , jg, kg)
                     - m_fvCoeffs.Cont.AU[X][sCU::cRight](i) * m_fields.U[X]( ig + sCU::iRight, jg, kg)

                     - m_fvCoeffs.Cont.AP[e](i, j, k) * m_fields.P( ig+1, jg  , kg  ) 
                     - m_fvCoeffs.Cont.AP[w](i, j, k) * m_fields.P( ig-1, jg  , kg  ) 

                     - m_fvCoeffs.Cont.AP[ee](i, j, k) * m_fields.P( ig+2, jg  , kg  ) 
                     - m_fvCoeffs.Cont.AP[ww](i, j, k) * m_fields.P( ig-2, jg  , kg  );


        // Update P from continuity
        m_fields.P( ig, jg, kg ) = ( 1 - m_fvCoeffs.Cont.relaxation ) * m_fieldsOld.P( ig, jg, kg )
                                 + m_fvCoeffs.Cont.relaxation * 
                                   ( bP 
                                   - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](i) * bU 
                                   - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](j) * bV 
                                   - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](k) * bW 
                                   ) * m_K(i, j, k);


        // Update U from momentum
        m_fields.U[X]( igU, jgU, kgU ) = ( 1 - m_fvCoeffs.Mom[X].relaxation ) * m_fieldsOld.U[X]( igU, jgU, kgU )
                                       + m_fvCoeffs.Mom[X].relaxation * ( bU - m_fvCoeffs.Mom[X].AP[sUP::cCoupled](iU) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[X].diagCoeffInv(iU, jU, kU) );

        // Update V from momentum
        m_fields.U[Y]( igV, jgV, kgV ) = ( 1 - m_fvCoeffs.Mom[Y].relaxation ) * m_fieldsOld.U[Y]( igV, jgV, kgV )
                                       + m_fvCoeffs.Mom[Y].relaxation * ( bV - m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jV) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Y].diagCoeffInv(iV, jV, kV) );

        // Update W from momentum
        m_fields.U[Z]( igW, jgW, kgW ) = ( 1 - m_fvCoeffs.Mom[Z].relaxation ) * m_fieldsOld.U[Z]( igW, jgW, kgW ) 
                                       + m_fvCoeffs.Mom[Z].relaxation * ( bW - m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kW) * m_fields.P( ig, jg, kg ) * m_fvCoeffs.Mom[Z].diagCoeffInv(iW, jW, kW) );

    }


    // Constants which are global to the linear solver
    void UpdateGlobalConstants()
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Axis::ENUMDATA;

        // Staggered indexing for fields
        intType iU, jU, kU,
                iV, jV, kV,
                iW, jW, kW;

        // Starting and ending indices, since K cannot be calculated on some boundaries due to the staggering
        intType iStart = 1 + sCU::iLeft,
                iLength = m_ni - 1 + sCU::iRight,

                jStart = 1 + sCV::iLeft,
                jLength = m_nj - 1 + sCV::iRight,

                kStart = 1 + sCW::iLeft,
                kLength = m_nk - 1 + sCW::iRight;

        for (intType k = kStart; k != kLength; k++) {

            kU = k;
            kV = k;
            kW = k + sCW::iCoupled;

            for (intType j = jStart; j != jLength; j++) {

                jU = j;
                jV = j + sCV::iCoupled;
                jW = j;

                #pragma GCC ivdep
                for (intType i = iStart; i != iLength; i++) {

                    iU = i + sCU::iCoupled;
                    iV = i;
                    iW = i;

                    m_K(i, j, k) = m_fvCoeffs.Cont.AP[p](i, j, k) 
                                 - m_fvCoeffs.Cont.AU[X][sCU::cCoupled](i) * m_fvCoeffs.Mom[X].AP[sUP::cCoupled](iU) * m_fvCoeffs.Mom[X].diagCoeffInv(iU, jU, kU) 
                                 - m_fvCoeffs.Cont.AU[Y][sCV::cCoupled](j) * m_fvCoeffs.Mom[Y].AP[sVP::cCoupled](jV) * m_fvCoeffs.Mom[Y].diagCoeffInv(iV, jV, kV) 
                                 - m_fvCoeffs.Cont.AU[Z][sCW::cCoupled](k) * m_fvCoeffs.Mom[Z].AP[sWP::cCoupled](kW) * m_fvCoeffs.Mom[Z].diagCoeffInv(iW, jW, kW);
                    m_K(i, j, k) = 1.0f / m_K(i, j, k);
                }
            }
        }
    }


private:
    FieldData<array3D> &m_fields;
    const FieldData<array3D> &m_fieldsOld;
    const FVCoefficients<MI> &m_fvCoeffs;
    const intType m_ni, m_nj, m_nk;
    array3D m_K;

};





// --------------------------------------------------------- LineSolver --------------------------------------------------------- //
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
    LineSolver( FieldData<array3D> &fields,
                const FieldData<array3D> &fieldsOld,
                const FVCoefficients<MI> &fvCoeffs) : 
                    m_fields( fields ),
                    m_fvCoeffs( fvCoeffs ),
                    m_lineConstants( array1D( fvCoeffs.nCells(Axis::X) ) ),
                    m_ni( fvCoeffs.nCells(Axis::X) )
    {
        if (m_ni == 1) {
            m_triadSolverCenter = std::make_unique<TriadSolver<TC::p, Vstag, Wstag, MI>>(fields, fieldsOld, fvCoeffs);
            SolutionUpdater = &LineSolver::Sweep2D;
            StateUpdater = &LineSolver::UpdateState2D;
        } else {
            m_triadSolverEast = std::make_unique<TriadSolver<TC::e, Vstag, Wstag, MI>>(fields, fieldsOld, fvCoeffs);
            m_triadSolverWest = std::make_unique<TriadSolver<TC::w, Vstag, Wstag, MI>>(fields, fieldsOld, fvCoeffs);
            SolutionUpdater = &LineSolver::Sweep3D;
            StateUpdater = &LineSolver::UpdateState3D;
        }
    }

    void SolveLine( const intType j, 
                    const intType k, 
                    const FieldData<array2D> &planeConstants)
    { (this->*SolutionUpdater)(j, k, planeConstants); }

    void UpdateState()
    { (this->*StateUpdater)(); }


private:
    FieldData<array3D> &m_fields;
    const FVCoefficients<MI> &m_fvCoeffs;

    std::unique_ptr<TriadSolver<TC::e, Vstag, Wstag, MI>> m_triadSolverEast;
    std::unique_ptr<TriadSolver<TC::w, Vstag, Wstag, MI>> m_triadSolverWest;
    std::unique_ptr<TriadSolver<TC::p, Vstag, Wstag, MI>> m_triadSolverCenter;

    FieldData<array1D> m_lineConstants;

    void (LineSolver::*SolutionUpdater)(const intType, const intType, const FieldData<array2D> &);
    void (LineSolver::*StateUpdater)(void);

    intType m_ni;

    // For 3D simulations
    void Sweep3D( const intType j, 
                  const intType k, 
                  const FieldData<array2D> &planeConstants)
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
                  const FieldData<array2D> &planeConstants )
    { 
        UpdateLineConstants(j, k, planeConstants);
        m_triadSolverCenter->UpdateTriad(0, j, k, m_lineConstants); 
    }

    void UpdateState2D()
    { m_triadSolverCenter->UpdateGlobalConstants(); }


    // Precalculate parts of stencil that are constant along a line
    void UpdateLineConstants( const intType j, 
                              const intType k, 
                              const FieldData<array2D> &planeConstants )
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

        // U momentum
        #pragma GCC ivdep
        for ( intType i = 0; i != m_ni; i++ ) {
            intType ig{ G(i) };

            // U momentum
            m_lineConstants.U[X](i) = planeConstants.U[X](i, j)
                                    + ( 
                                      - m_fvCoeffs.Mom[X].AU[X][n](i, j, k) * m_fields.U[X]( ig  , jg+1, kg  )
                                      - m_fvCoeffs.Mom[X].AU[X][s](i, j, k) * m_fields.U[X]( ig  , jg-1, kg  )
                                      );

            // V momentum
            m_lineConstants.U[Y](i) = planeConstants.U[Y](i, jV)
                                    + (
                                      - m_fvCoeffs.Mom[Y].AU[Y][n](i, jV, kV) * m_fields.U[Y]( ig  , jgV+1, kgV  )  
                                      - m_fvCoeffs.Mom[Y].AU[Y][s](i, jV, kV) * m_fields.U[Y]( ig  , jgV-1, kgV  ) 

                                      - m_fvCoeffs.Mom[Y].AP[sVP::cLeft ](jV) * m_fields.P( ig, jgV + sVP::iLeft , kgV)
                                      - m_fvCoeffs.Mom[Y].AP[sVP::cRight](jV) * m_fields.P( ig, jgV + sVP::iRight, kgV)
                                      );

            // W momentum
            m_lineConstants.U[Z](i) = planeConstants.U[Z](i, jW)
                                    + ( 
                                      - m_fvCoeffs.Mom[Z].AU[Z][n](i, jW, kW) * m_fields.U[Z]( ig  , jgW+1, kgW  ) 
                                      - m_fvCoeffs.Mom[Z].AU[Z][s](i, jW, kW) * m_fields.U[Z]( ig  , jgW-1, kgW  ) 
                                      );

            // Continuity equation
            m_lineConstants.P(i) = planeConstants.P(i, j)

                                 - m_fvCoeffs.Cont.AU[Y][sCV::cLeft ](j) * m_fields.U[Y]( ig, jg + sCV::iLeft , kg)
                                 - m_fvCoeffs.Cont.AU[Y][sCV::cRight](j) * m_fields.U[Y]( ig, jg + sCV::iRight, kg)

                                 - m_fvCoeffs.Cont.AP[n](i, j, k) * m_fields.P( ig  , jg+1, kg  ) 
                                 - m_fvCoeffs.Cont.AP[s](i, j, k) * m_fields.P( ig  , jg-1, kg  ) 
 
                                 - m_fvCoeffs.Cont.AP[nn](i, j, k) * m_fields.P( ig  , jg+2, kg  )
                                 - m_fvCoeffs.Cont.AP[ss](i, j, k) * m_fields.P( ig  , jg-2, kg  );
        }

    }

};





// --------------------------------------------------------- PlaneSolver --------------------------------------------------------- //
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
                 const FVCoefficients<MI> &fvCoeffs) : 
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
    const FVCoefficients<MI> &m_fvCoeffs;

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

        #pragma GCC ivdep
        for ( intType j = 0; j != m_nj; j++ ) {

            #pragma GCC ivdep
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
                m_planeConstants.P(i, j) = m_fvCoeffs.Cont.B(i, j, k)

                                         - m_fvCoeffs.Cont.AU[Z][sCW::cLeft ](k) * m_fields.U[Z]( ig, jg, kg + sCW::iLeft )
                                         - m_fvCoeffs.Cont.AU[Z][sCW::cRight](k) * m_fields.U[Z]( ig, jg, kg + sCW::iRight)

                                         - m_fvCoeffs.Cont.AP[t](i, j, k) * m_fields.P( ig  , jg  , kg+1) 
                                         - m_fvCoeffs.Cont.AP[b](i, j, k) * m_fields.P( ig  , jg  , kg-1)
        
                                         - m_fvCoeffs.Cont.AP[tt](i, j, k) * m_fields.P( ig  , jg  , kg+2) 
                                         - m_fvCoeffs.Cont.AP[bb](i, j, k) * m_fields.P( ig  , jg  , kg-2);
            }
        }

    }
};





// --------------------------------------------------------- LinearSolver --------------------------------------------------------- //
template< MomentumInterpolation MI >
class LinearSolver
{
    using TC = TransportCoefficients::ENUMDATA;
    using A = Axis::ENUMDATA;

public:
    LinearSolver( FieldData<array3D> &fields,
                  const FieldData<array3D> &fieldsOld,
                  const FVCoefficients<MI> &fvCoeffs, 
                  const InputData::LinearSolverSettings &linearSolverSettings) : 
                    m_fields( fields ),
                    m_maxIterations( linearSolverSettings.maxIterations ),
                    m_maxResiduals( linearSolverSettings.maxResiduals ),
                    m_relaxation( linearSolverSettings.relaxation ),

                    m_delta( array2D( m_fields.P.dimension(A::X), m_fields.P.dimension(A::Y) ) ),
                    m_oldPlane( array2D( m_fields.P.dimension(A::X), m_fields.P.dimension(A::Y) ) ),

                    m_ni( fvCoeffs.nCells(A::X) ),
                    m_nj( fvCoeffs.nCells(A::Y) ),
                    m_nk( fvCoeffs.nCells(A::Z) )
    {
        if (m_nk == 1) {
            m_planeSolverCenter = std::make_unique<PlaneSolver<TC::p, MI>>(fields, fieldsOld, fvCoeffs);
            SolutionUpdater = &LinearSolver::Sweep2D;
            StateUpdater = &LinearSolver::UpdateState2D;
        } else {
            m_planeSolverTop = std::make_unique<PlaneSolver<TC::t, MI>>(fields, fieldsOld, fvCoeffs);
            m_planeSolverBottom = std::make_unique<PlaneSolver<TC::b, MI>>(fields, fieldsOld, fvCoeffs);
            SolutionUpdater = &LinearSolver::Sweep3D;
            StateUpdater = &LinearSolver::UpdateState3D;
        }
    }


    void Solve()
    {
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        for ( intType nIterations = 1; nIterations <= m_maxIterations; nIterations++ )
        {
            // Reset residuals
            ForAllFieldData( [&] (intType f) { m_residuals[f] = 0.0f; });

            // Update plane
            (this->*SolutionUpdater)();

            // Normalise residuals
            ForAllFieldData( [&] (intType f) { m_residuals[f] /= static_cast<floatType>(m_ni * m_nj * m_nk); });
            NormaliseResiduals(m_residuals, m_residualsInitialInv, nIterations);

            // Check residual tolerence
            if ( MetResidualTolerence(m_residuals, m_maxResiduals) ) {
                std::cout << "*** INNER ITERATIONS CONVERGED ***"
                            << "\n\n";
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

    FieldData<array3D> &m_fields;
    const intType m_maxIterations;
    const FieldData<floatType> m_maxResiduals;
    const FieldData<floatType> m_relaxation;

    std::unique_ptr<PlaneSolver<TC::t, MI>> m_planeSolverTop;
    std::unique_ptr<PlaneSolver<TC::b, MI>> m_planeSolverBottom;
    std::unique_ptr<PlaneSolver<TC::p, MI>> m_planeSolverCenter;

    void (LinearSolver::*SolutionUpdater)(void);
    void (LinearSolver::*StateUpdater)(void);

    FieldData<array2D> m_delta, m_oldPlane;
    FieldData<floatType> m_residuals, m_residualsInitialInv;
    FieldData<intType> m_kS;

    intType m_ni, m_nj, m_nk;

    // For 3D simulations
    void Sweep3D()
    {
        for (intType k = 0; k != m_nk - 1; k++) { // Forward sweep
            UpdateAndRelax(m_planeSolverTop, k);
        }

        for (intType k = m_nk - 1; k != 0; k--) { // Backward sweep
            UpdateAndRelax(m_planeSolverBottom, k);
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
        UpdateAndRelax(m_planeSolverCenter, 0);
    }

    void UpdateState2D()
    {
        m_planeSolverCenter->UpdateState();
    }


    template <TC Wstag>
    void UpdateAndRelax( std::unique_ptr<PlaneSolver<Wstag, MI>> &planeSolver, intType k )
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Axis::ENUMDATA;

        ForAllFieldData( [&] (intType f) { m_kS[f] = k; });  // Set iterating coefficient
        m_kS.U[Z] += LUT::CoeffIndex[Wstag];                 // W momentum is staggered

        ForAllFieldData( [&] (intType f) { m_oldPlane[f] = m_fields[f].chip( G(m_kS[f]), Z ); } ); // Set old plane

        planeSolver->SolvePlane(k);

        ForAllFieldData( [&] (intType f) {
            auto fieldPlane = m_fields[f].chip( G(m_kS[f]), Z );
            m_delta[f] = m_delta[f].constant( m_relaxation[f] ) * (fieldPlane - m_oldPlane[f]); // Relaxed change in plane
            fieldPlane = m_oldPlane[f] + m_delta[f];                                            // Relax
            m_residuals[f] += static_cast<array0D>( m_delta[f].abs().sum() )(0);                // Add to residual count
        } );
    }
};





// --------------------------------------------------------- SweepSolve --------------------------------------------------------- //

template< MomentumInterpolation MI >
void SweepSolve( FieldData<array3D> &fields,
                 const Mesh &mesh,
                 const InputData &inputData,
                 const AxisTransformationMap &axisTransformation)
{
    using enum Axis::ENUMDATA;

    // Extract from input data
    const InputData::LinearSolverSettings linearSolverSettings = inputData.linearSolverSettings;
    const intType maxOuterIterations = inputData.schemes.maxOuterIterations;
    const FieldData<floatType> maxOuterResiduals = inputData.schemes.maxOuterResiduals;

    // Initialise
    EnumVector<Axis, array3D> faceFluxes = InitialiseFaceFluxes(mesh, fields.U, inputData);
    FieldData<array3D> fieldsOld( fields );
    FVCoefficients<MI> fvCoeffs = InitialiseFVCoefficients<MI>(mesh, faceFluxes, inputData);

    // Initialise residuals
    FieldData<floatType> residualsOuter, residualsScaleFactor;
    floatType massFluxResidual;

    // Logging objects
    std::vector< FieldProbe > fieldProbes;
    std::vector< ProbeLogFile > probeLogFiles;
    for ( const auto &probeData : inputData.probes ) {
        fieldProbes.emplace_back( mesh, probeData.location );
        probeLogFiles.emplace_back( probeData.filename, axisTransformation, fieldProbes.back() );
    }
    std::vector< FieldData<floatType> > probeValues( fieldProbes.size() );
    
    ResidualLogFile residualsLogFile( inputData.residualHistoryFilename, axisTransformation );
    ConsoleLog consoleLog( axisTransformation );

    // Instantiate linear solver, this holds references to the fields
    LinearSolver<MI> linearSolver(fields, fieldsOld, fvCoeffs, linearSolverSettings);


    // Outer iterations
    for ( intType nOuterIterations = 1; nOuterIterations <= maxOuterIterations; nOuterIterations++ )
    {
        linearSolver.UpdateState();
        linearSolver.Solve();

        UpdateFaceFluxes(faceFluxes, mesh, fields.U, inputData);
        UpdateFVCoefficients(fvCoeffs, mesh, faceFluxes, inputData);

        // residualsOuter   = L1DiffResiduals(fields, fieldsOld);
        residualsOuter   = StencilResiduals(fields, fvCoeffs); 
        NormaliseResiduals( residualsOuter, residualsScaleFactor, nOuterIterations );

        massFluxResidual = BoundaryMassFluxResidual(faceFluxes, mesh);

        probeValues      = FieldProbeValues(fields, fieldProbes); 
        

        fieldsOld = fields;

        consoleLog.WriteResiduals( residualsOuter, massFluxResidual, nOuterIterations );
        residualsLogFile.WriteData( residualsOuter, massFluxResidual, nOuterIterations );
        for ( size_t p = 0; p != fieldProbes.size(); p++ ) {
            probeLogFiles[p].WriteData( probeValues[p], nOuterIterations );
        }

        if ( MetResidualTolerence(residualsOuter, maxOuterResiduals) ) {
            std::cout << "*** OUTER ITERATIONS CONVERGED ***"
                        << "\n\n";
            break;
        }
    }
}
template void SweepSolve<MomentumInterpolation::Implicit>( FieldData<array3D> &, const Mesh &, const InputData &, const AxisTransformationMap &);
template void SweepSolve<MomentumInterpolation::SemiExplicit>( FieldData<array3D> &, const Mesh &, const InputData &, const AxisTransformationMap &);

} // end namespace CFD