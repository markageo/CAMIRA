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
    template <typename arrayType, typename enumStruct> [[ maybe_unused ]]
    EnumVector<enumStruct, floatType> L1ArrayDiff( const ArrayAllocator<enumStruct, arrayType> &array1,
                                                   const ArrayAllocator<enumStruct, arrayType> &array2)
    {
        static_assert( std::is_same<arrayType, array1D>::value ||
                       std::is_same<arrayType, array2D>::value ||
                       std::is_same<arrayType, array3D>::value);

        EnumVector<enumStruct, floatType> result;

        EnumFor<enumStruct>([&](typename enumStruct::ENUMDATA enumName) { 
            auto fieldDiff = array1[enumName] - array2[enumName];  // auto lazily evaluates
            result[enumName] = static_cast<array0D>( fieldDiff.abs().mean() )(0); 
        });

        return result;
    }



    // Calculate the absolute residual of each equation from the finite volume stencil
    [[ maybe_unused ]]
    EnumVector<Fields, floatType> StencilResiduals( const ArrayAllocator<Fields, array3D> &fields,
                                                    const FVCoefficients &fvCoeffs )
    {
        using enum Fields::ENUMDATA;
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        EnumVector<CFD::Fields, floatType> residuals;

        for ( intType k = 0; k != fvCoeffs.nCells[Z]; k++ ) {
            for ( intType j = 0; j != fvCoeffs.nCells[Y]; j++ ) {
                for ( intType i = 0; i != fvCoeffs.nCells[Z]; i++ ) {

                    // U momentum
                    residuals[U] += abs(
                                      fvCoeffs.Umom.AU[p](i, j, k) * fields[U]( G(i  , j  , k  ) )
                                    + fvCoeffs.Umom.AU[n](i, j, k) * fields[U]( G(i  , j+1, k  ) )
                                    + fvCoeffs.Umom.AU[e](i, j, k) * fields[U]( G(i+1, j  , k  ) )
                                    + fvCoeffs.Umom.AU[s](i, j, k) * fields[U]( G(i  , j-1, k  ) )
                                    + fvCoeffs.Umom.AU[w](i, j, k) * fields[U]( G(i-1, j  , k  ) )
                                    + fvCoeffs.Umom.AU[t](i, j, k) * fields[U]( G(i  , j  , k+1) )
                                    + fvCoeffs.Umom.AU[b](i, j, k) * fields[U]( G(i  , j  , k-1) )

                                    + fvCoeffs.Umom.AP[e](i) * fields[P]( G(i+1, j  , k  ) )
                                    + fvCoeffs.Umom.AP[p](i) * fields[P]( G(i  , j  , k  ) )
                                    + fvCoeffs.Umom.AP[w](i) * fields[P]( G(i-1, j  , k  ) )

                                    - fvCoeffs.Umom.B(i, j, k) );

                    // V momentum
                    residuals[V] += abs( 
                                      fvCoeffs.Vmom.AV[p](i, j, k) * fields[V]( G(i  , j  , k  ) )
                                    + fvCoeffs.Vmom.AV[n](i, j, k) * fields[V]( G(i  , j+1, k  ) )
                                    + fvCoeffs.Vmom.AV[e](i, j, k) * fields[V]( G(i+1, j  , k  ) )
                                    + fvCoeffs.Vmom.AV[s](i, j, k) * fields[V]( G(i  , j-1, k  ) )
                                    + fvCoeffs.Vmom.AV[w](i, j, k) * fields[V]( G(i-1, j  , k  ) )
                                    + fvCoeffs.Vmom.AV[t](i, j, k) * fields[V]( G(i  , j  , k+1) )
                                    + fvCoeffs.Vmom.AV[b](i, j, k) * fields[V]( G(i  , j  , k-1) )

                                    + fvCoeffs.Vmom.AP[n](j) * fields[P]( G(i  , j+1, k  ) )
                                    + fvCoeffs.Vmom.AP[p](j) * fields[P]( G(i  , j  , k  ) )
                                    + fvCoeffs.Vmom.AP[s](j) * fields[P]( G(i  , j-1, k  ) )

                                    - fvCoeffs.Vmom.B(i, j, k) );

                    // W momentum
                    residuals[W] += abs( 
                                      fvCoeffs.Wmom.AW[p](i, j, k) * fields[W]( G(i  , j  , k  ) )
                                    + fvCoeffs.Wmom.AW[n](i, j, k) * fields[W]( G(i  , j+1, k  ) )
                                    + fvCoeffs.Wmom.AW[e](i, j, k) * fields[W]( G(i+1, j  , k  ) )
                                    + fvCoeffs.Wmom.AW[s](i, j, k) * fields[W]( G(i  , j-1, k  ) )
                                    + fvCoeffs.Wmom.AW[w](i, j, k) * fields[W]( G(i-1, j  , k  ) )
                                    + fvCoeffs.Wmom.AW[t](i, j, k) * fields[W]( G(i  , j  , k+1) )
                                    + fvCoeffs.Wmom.AW[b](i, j, k) * fields[W]( G(i  , j  , k-1) )

                                    + fvCoeffs.Wmom.AP[t](k) * fields[P]( G(i  , j  , k+1) )
                                    + fvCoeffs.Wmom.AP[p](k) * fields[P]( G(i  , j  , k  ) )
                                    + fvCoeffs.Wmom.AP[b](k) * fields[P]( G(i  , j  , k-1) )

                                    - fvCoeffs.Wmom.B(i, j, k) );

                    // Continuity
                    residuals[P] += abs( 
                                      fvCoeffs.Cont.AU[e](i) * fields[U]( G(i+1, j  , k  ) )
                                    + fvCoeffs.Cont.AU[p](i) * fields[U]( G(i  , j  , k  ) )
                                    + fvCoeffs.Cont.AU[w](i) * fields[U]( G(i-1, j  , k  ) )

                                    + fvCoeffs.Cont.AV[n](j) * fields[V]( G(i  , j+1, k  ) )
                                    + fvCoeffs.Cont.AV[p](j) * fields[V]( G(i  , j  , k  ) )
                                    + fvCoeffs.Cont.AV[s](j) * fields[V]( G(i  , j-1, k  ) )

                                    + fvCoeffs.Cont.AW[t](k) * fields[W]( G(i  , j  , k+1) )
                                    + fvCoeffs.Cont.AW[p](k) * fields[W]( G(i  , j  , k  ) )
                                    + fvCoeffs.Cont.AW[b](k) * fields[W]( G(i  , j  , k-1) )

                                    + fvCoeffs.Cont.AP[p](i, j, k) * fields[P]( G(i  , j  , k  ) )
                                    + fvCoeffs.Cont.AP[n](i, j, k) * fields[P]( G(i  , j+1, k  ) ) 
                                    + fvCoeffs.Cont.AP[e](i, j, k) * fields[P]( G(i+1, j  , k  ) ) 
                                    + fvCoeffs.Cont.AP[s](i, j, k) * fields[P]( G(i  , j-1, k  ) ) 
                                    + fvCoeffs.Cont.AP[w](i, j, k) * fields[P]( G(i-1, j  , k  ) ) 
                                    + fvCoeffs.Cont.AP[t](i, j, k) * fields[P]( G(i  , j  , k+1) ) 
                                    + fvCoeffs.Cont.AP[b](i, j, k) * fields[P]( G(i  , j  , k-1) )

                                    + fvCoeffs.Cont.AP[nn](i, j, k) * fields[P]( G(i  , j+2, k  ) ) 
                                    + fvCoeffs.Cont.AP[ee](i, j, k) * fields[P]( G(i+2, j  , k  ) ) 
                                    + fvCoeffs.Cont.AP[ss](i, j, k) * fields[P]( G(i  , j-2, k  ) ) 
                                    + fvCoeffs.Cont.AP[ww](i, j, k) * fields[P]( G(i-2, j  , k  ) ) 
                                    + fvCoeffs.Cont.AP[tt](i, j, k) * fields[P]( G(i  , j  , k+2) ) 
                                    + fvCoeffs.Cont.AP[bb](i, j, k) * fields[P]( G(i  , j  , k-2) )
                                    
                                    - fvCoeffs.Cont.B(i, j, k) );

                }
            }
        }

        EnumFor<Fields>( [&] (Fields::ENUMDATA field) {
            residuals[field] /= static_cast<floatType>( fvCoeffs.nCells[X] * fvCoeffs.nCells[Y] * fvCoeffs.nCells[Z] ); 
        } );

        return residuals;
    }



    // Calculate global mass flux residual at the domain boundary
    [[ maybe_unused ]]
    floatType BoundaryMassFluxResidual( const ArrayAllocator<Fields, array3D> &faceVelocities,
                                        const Mesh &mesh )
    {
        floatType massFluxResidual = 0.0f;

        EnumFor<Axis>([&](Axis::ENUMDATA axis) {

            // Positive face, area normal is in positive direction
            auto faceFluxesPositive = faceVelocities[AxisVelocity[axis]].chip( mesh.nCells(axis), axis ) * mesh.cellFaceAreas[axis];
            massFluxResidual += static_cast<array0D>( faceFluxesPositive.sum() )(0);

            // Negative face, area normal is in negative direction
            auto faceFluxesNegative = -faceVelocities[AxisVelocity[axis]].chip( 0, axis ) * mesh.cellFaceAreas[axis];
            massFluxResidual += static_cast<array0D>( faceFluxesNegative.sum() )(0);

        });

        return massFluxResidual;
    }



    // Turn the residual into a relative residual
    [[ maybe_unused ]]
    void RelativeResidual( EnumVector<Fields, floatType> &residuals,
                           EnumVector<Fields, floatType> &residualsInitialInv,
                           const intType nIterations )
    {
        if (nIterations == 1) {

            EnumFor<Fields>([&](Fields::ENUMDATA field) {
                if ( residuals[field] != 0 ) { // Division by zero
                    residualsInitialInv[field] = 1.0f / residuals[field];
                } else {
                    residualsInitialInv[field] = 1.0f;
                }
            });

        }

        EnumFor<Fields>([&](Fields::ENUMDATA field) { residuals[field] *= residualsInitialInv[field]; });
    }



    // Check if residual tolerence is met
    bool MetResidualTolerence( const EnumVector<Fields, floatType> &residuals,
                               const EnumVector<Fields, floatType> &residualsTarget )
    {
        bool met = true;
        EnumFor<Fields>( [&] (Fields::ENUMDATA field) {

            if (residuals[field] > residualsTarget[field]) 
                met = false;
                
        } );
        return met;
    }

} // end anonymous namespace





// Performs a single local update of block coupled equations
template < TransportCoefficients::ENUMDATA Ustag,
           TransportCoefficients::ENUMDATA Vstag,
           TransportCoefficients::ENUMDATA Wstag>
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
    TriadSolver( ArrayAllocator<Fields, array3D> &fields,
                 const ArrayAllocator<Fields, array3D> &fieldsOld,
                 const FVCoefficients &fvCoeffs ) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_fvCoeffs( fvCoeffs ),
                    m_ni( fvCoeffs.nCells(0) ),
                    m_nj( fvCoeffs.nCells(1) ),
                    m_nk (fvCoeffs.nCells(2) ),
                    m_K( array3D(m_ni, m_nj, m_nk).setZero() )
    {
        UpdateGlobalConstants();
    };


    // Core function which updates the local coupled system. Templated by staggering direction.
    void UpdateTriad( const intType i, 
                      const intType j, 
                      const intType k )
    {
        using enum Fields::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        // For indexing the staggered cells
        intType iU( i + sCU::iCoupled ), jU( j                 ), kU( k                 ); // U momentum
        intType iV( i                 ), jV( j + sCV::iCoupled ), kV( k                 ); // V momentum
        intType iW( i                 ), jW( j                 ), kW( k + sCW::iCoupled ); // W momentum


        // Precompute momentum RHS divided by AP coefficients
        // U momentum
        floatType bU = ( m_fvCoeffs.Umom.B(iU, jU, kU)

                       - m_fvCoeffs.Umom.AU[n](iU, jU, kU) * m_fields[U]( G(iU  , jU+1, kU  ) ) 
                       - m_fvCoeffs.Umom.AU[e](iU, jU, kU) * m_fields[U]( G(iU+1, jU  , kU  ) ) 
                       - m_fvCoeffs.Umom.AU[s](iU, jU, kU) * m_fields[U]( G(iU  , jU-1, kU  ) ) 
                       - m_fvCoeffs.Umom.AU[w](iU, jU, kU) * m_fields[U]( G(iU-1, jU  , kU  ) ) 
                       - m_fvCoeffs.Umom.AU[t](iU, jU, kU) * m_fields[U]( G(iU  , jU  , kU+1) ) 
                       - m_fvCoeffs.Umom.AU[b](iU, jU, kU) * m_fields[U]( G(iU  , jU  , kU-1) )

                       - m_fvCoeffs.Umom.AP[sUP::cLeft ](iU) * m_fields[P]( G(iU + sUP::iLeft , jU, kU) ) 
                       - m_fvCoeffs.Umom.AP[sUP::cRight](iU) * m_fields[P]( G(iU + sUP::iRight, jU, kU) ) 

                       ) * m_fvCoeffs.Umom.diagCoeffInv(iU, jU, kU);


        // V momentum
        floatType bV = ( m_fvCoeffs.Vmom.B(iV, jV, kV)

                       - m_fvCoeffs.Vmom.AV[n](iV, jV, kV) * m_fields[V]( G(iV  , jV+1, kV  ) ) 
                       - m_fvCoeffs.Vmom.AV[e](iV, jV, kV) * m_fields[V]( G(iV+1, jV  , kV  ) ) 
                       - m_fvCoeffs.Vmom.AV[s](iV, jV, kV) * m_fields[V]( G(iV  , jV-1, kV  ) ) 
                       - m_fvCoeffs.Vmom.AV[w](iV, jV, kV) * m_fields[V]( G(iV-1, jV  , kV  ) ) 
                       - m_fvCoeffs.Vmom.AV[t](iV, jV, kV) * m_fields[V]( G(iV  , jV  , kV+1) ) 
                       - m_fvCoeffs.Vmom.AV[b](iV, jV, kV) * m_fields[V]( G(iV  , jV  , kV-1) )

                       - m_fvCoeffs.Vmom.AP[sVP::cLeft ](jV) * m_fields[P]( G(iV, jV + sVP::iLeft , kV) ) 
                       - m_fvCoeffs.Vmom.AP[sVP::cRight](jV) * m_fields[P]( G(iV, jV + sVP::iRight, kV) )

                       ) * m_fvCoeffs.Vmom.diagCoeffInv(iV, jV, kV);


        // W momentum
        floatType bW = ( m_fvCoeffs.Wmom.B(iW, jW, kW)

                       - m_fvCoeffs.Wmom.AW[n](iW, jW, kW) * m_fields[W]( G(iW  , jW+1, kW  ) ) 
                       - m_fvCoeffs.Wmom.AW[e](iW, jW, kW) * m_fields[W]( G(iW+1, jW  , kW  ) ) 
                       - m_fvCoeffs.Wmom.AW[s](iW, jW, kW) * m_fields[W]( G(iW  , jW-1, kW  ) ) 
                       - m_fvCoeffs.Wmom.AW[w](iW, jW, kW) * m_fields[W]( G(iW-1, jW  , kW  ) ) 
                       - m_fvCoeffs.Wmom.AW[t](iW, jW, kW) * m_fields[W]( G(iW  , jW  , kW+1) ) 
                       - m_fvCoeffs.Wmom.AW[b](iW, jW, kW) * m_fields[W]( G(iW  , jW  , kW-1) )

                       - m_fvCoeffs.Wmom.AP[sWP::cLeft ](kW) * m_fields[P]( G(iW, jW, kW + sWP::iLeft) ) 
                       - m_fvCoeffs.Wmom.AP[sWP::cRight](kW) * m_fields[P]( G(iW, jW, kW + sWP::iRight) )

                       ) * m_fvCoeffs.Wmom.diagCoeffInv(iW, jW, kW);


        // Continuity for pressure
        floatType bP = m_fvCoeffs.Cont.B(i, j, k)

                     - m_fvCoeffs.Cont.AU[sCU::cLeft ](i) * m_fields[U]( G(i + sCU::iLeft , j, k) ) 
                     - m_fvCoeffs.Cont.AU[sCU::cRight](i) * m_fields[U]( G(i + sCU::iRight, j, k) )

                     - m_fvCoeffs.Cont.AV[sCV::cLeft ](j) * m_fields[V]( G(i, j + sCV::iLeft , k) ) 
                     - m_fvCoeffs.Cont.AV[sCV::cRight](j) * m_fields[V]( G(i, j + sCV::iRight, k) )

                     - m_fvCoeffs.Cont.AW[sCW::cLeft ](k) * m_fields[W]( G(i, j, k + sCW::iLeft ) )
                     - m_fvCoeffs.Cont.AW[sCW::cRight](k) * m_fields[W]( G(i, j, k + sCW::iRight) )

                     - m_fvCoeffs.Cont.AP[n](i, j, k) * m_fields[P]( G(i  , j+1, k  )) 
                     - m_fvCoeffs.Cont.AP[e](i, j, k) * m_fields[P]( G(i+1, j  , k  )) 
                     - m_fvCoeffs.Cont.AP[s](i, j, k) * m_fields[P]( G(i  , j-1, k  )) 
                     - m_fvCoeffs.Cont.AP[w](i, j, k) * m_fields[P]( G(i-1, j  , k  )) 
                     - m_fvCoeffs.Cont.AP[t](i, j, k) * m_fields[P]( G(i  , j  , k+1)) 
                     - m_fvCoeffs.Cont.AP[b](i, j, k) * m_fields[P]( G(i  , j  , k-1))

                     - m_fvCoeffs.Cont.AP[nn](i, j, k) * m_fields[P]( G(i  , j+2, k  ) ) 
                     - m_fvCoeffs.Cont.AP[ee](i, j, k) * m_fields[P]( G(i+2, j  , k  ) ) 
                     - m_fvCoeffs.Cont.AP[ss](i, j, k) * m_fields[P]( G(i  , j-2, k  ) ) 
                     - m_fvCoeffs.Cont.AP[ww](i, j, k) * m_fields[P]( G(i-2, j  , k  ) ) 
                     - m_fvCoeffs.Cont.AP[tt](i, j, k) * m_fields[P]( G(i  , j  , k+2) ) 
                     - m_fvCoeffs.Cont.AP[bb](i, j, k) * m_fields[P]( G(i  , j  , k-2) );


        // Update P from continuity
        m_fields[P]( G(i, j, k) ) = ( 1 - m_fvCoeffs.Cont.relaxation ) * m_fieldsOld[P]( G(i, j, k) )
                                  + m_fvCoeffs.Cont.relaxation * 
                                    ( bP 
                                    - m_fvCoeffs.Cont.AU[sCU::cCoupled](i) * bU 
                                    - m_fvCoeffs.Cont.AV[sCV::cCoupled](j) * bV 
                                    - m_fvCoeffs.Cont.AW[sCW::cCoupled](k) * bW 
                                    ) * m_K(i, j, k);


        // Update U from momentum
        m_fields[U]( G(iU, jU, kU) ) = ( 1 - m_fvCoeffs.Umom.relaxation ) * m_fieldsOld[U]( G(iU, jU, kU) )
                                     + m_fvCoeffs.Umom.relaxation * ( bU - m_fvCoeffs.Umom.AP[sUP::cCoupled](iU) * m_fields[P]( G(i, j, k) ) * m_fvCoeffs.Umom.diagCoeffInv(iU, jU, kU) );

        // Update V from momentum
        m_fields[V]( G(iV, jV, kV) ) = ( 1 - m_fvCoeffs.Vmom.relaxation ) * m_fieldsOld[V]( G(iV, jV, kV) )
                                     + m_fvCoeffs.Vmom.relaxation * ( bV - m_fvCoeffs.Vmom.AP[sVP::cCoupled](jV) * m_fields[P]( G(i, j, k) ) * m_fvCoeffs.Vmom.diagCoeffInv(iV, jV, kV) );

        // Update W from momentum
        m_fields[W]( G(iW, jW, kW) ) = ( 1 - m_fvCoeffs.Wmom.relaxation ) * m_fieldsOld[W]( G(iW, jW, kW) ) 
                                     + m_fvCoeffs.Wmom.relaxation * ( bW - m_fvCoeffs.Wmom.AP[sWP::cCoupled](kW) * m_fields[P]( G(i, j, k) ) * m_fvCoeffs.Wmom.diagCoeffInv(iW, jW, kW) );

    }

    // Constants which are global to the linear solver
    void UpdateGlobalConstants()
    {
        using enum TransportCoefficients::ENUMDATA;

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

                for (intType i = iStart; i != iLength; i++) {

                    iU = i + sCU::iCoupled;
                    iV = i;
                    iW = i;

                    m_K(i, j, k) = m_fvCoeffs.Cont.AP[p](i, j, k) 
                                 - m_fvCoeffs.Cont.AU[sCU::cCoupled](i) * m_fvCoeffs.Umom.AP[sUP::cCoupled](iU) * m_fvCoeffs.Umom.diagCoeffInv(iU, jU, kU) 
                                 - m_fvCoeffs.Cont.AV[sCV::cCoupled](j) * m_fvCoeffs.Vmom.AP[sVP::cCoupled](jV) * m_fvCoeffs.Vmom.diagCoeffInv(iV, jV, kV) 
                                 - m_fvCoeffs.Cont.AW[sCW::cCoupled](k) * m_fvCoeffs.Wmom.AP[sWP::cCoupled](kW) * m_fvCoeffs.Wmom.diagCoeffInv(iW, jW, kW);
                    m_K(i, j, k) = 1.0f / m_K(i, j, k);
                }
            }
        }
    }


private:
    ArrayAllocator<Fields, array3D> &m_fields;
    const ArrayAllocator<Fields, array3D> &m_fieldsOld;
    const FVCoefficients &m_fvCoeffs;
    intType m_ni, m_nj, m_nk;
    array3D m_K;

};





template <TransportCoefficients::ENUMDATA Vstag,
            TransportCoefficients::ENUMDATA Wstag>
class LineSolver
{
    using TC = TransportCoefficients::ENUMDATA;

    // Staggering must be valid
    static_assert( (Vstag == TC::n) || (Vstag == TC::s) || (Vstag == TC::p), "Invalid V momentum staggering" );
    static_assert( (Wstag == TC::t) || (Wstag == TC::b) || (Wstag == TC::p), "Invalid W momentum staggering" );

public:
    LineSolver( ArrayAllocator<Fields, array3D> &fields,
                const ArrayAllocator<Fields, array3D> &fieldsOld,
                const FVCoefficients &fvCoeffs,
                const InputData::LineSolverSettings &lineSolverSettings) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_maxIterations( lineSolverSettings.maxIterations ),
                    m_maxResiduals( lineSolverSettings.maxResiduals ),
                    m_relaxation( lineSolverSettings.relaxation ),
                    m_ni( fvCoeffs.nCells(0) )
    {
        if (m_ni == 1) {

            m_blockSolverCenter = std::make_unique<TriadSolver<TC::p, Vstag, Wstag>>(fields, fieldsOld, fvCoeffs);
            SolutionUpdater = &LineSolver::Sweep2D;
            StateUpdater = &LineSolver::UpdateState2D;

        } else {
            m_blockSolverEast = std::make_unique<TriadSolver<TC::e, Vstag, Wstag>>(fields, fieldsOld, fvCoeffs);
            m_blockSolverWest = std::make_unique<TriadSolver<TC::w, Vstag, Wstag>>(fields, fieldsOld, fvCoeffs);
            SolutionUpdater = &LineSolver::Sweep3D;
            StateUpdater = &LineSolver::UpdateState3D;

        }
    }


    void SolveLine(const intType j, const intType k)
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Fields::ENUMDATA;
        using enum Axis::ENUMDATA;

        // Staggered indexing
        EnumFor<Fields>( [&] (Fields::ENUMDATA f) {
            m_jS[f] = j; 
            m_kS[f] = k; 
        } );
        m_kS[V] += CoeffIndex[Vstag];
        m_kS[W] += CoeffIndex[Wstag];

        // Solver loop
        for ( intType nIterations = 1; nIterations <= m_maxIterations; nIterations++ )
        {
            EnumFor<Fields>( [&] (Fields::ENUMDATA field) { m_residuals[field] = 0.0f; });

            // Update block
            (this->*SolutionUpdater)(j, k);

            // Normalise residuals
            EnumFor<Fields>( [&] (Fields::ENUMDATA field) { m_residuals[field] /= static_cast<floatType>(m_ni); } );
            RelativeResidual( m_residuals, m_residualsInitialInv, nIterations );

            // Check residual
            if (MetResidualTolerence(m_residuals, m_maxResiduals)) {
                break;
            }
        }
    }

    void UpdateState()
    {
        (this->*StateUpdater)();
    }


private:
    ArrayAllocator<Fields, array3D> &m_fields;
    const ArrayAllocator<Fields, array3D> &m_fieldsOld;
    const intType m_maxIterations;
    const EnumVector<Fields, floatType> m_maxResiduals;
    const EnumVector<Fields, floatType> m_relaxation;

    std::unique_ptr<TriadSolver<TC::e, Vstag, Wstag>> m_blockSolverEast;
    std::unique_ptr<TriadSolver<TC::w, Vstag, Wstag>> m_blockSolverWest;
    std::unique_ptr<TriadSolver<TC::p, Vstag, Wstag>> m_blockSolverCenter;

    void (LineSolver::*SolutionUpdater)(intType, intType);
    void (LineSolver::*StateUpdater)(void);

    EnumVector<Fields, floatType> m_oldBlock, m_delta, m_residuals, m_residualsInitialInv;
    EnumVector<Fields, intType> m_iS, m_jS, m_kS;
    intType m_ni;

    // For 3D simulations
    void Sweep3D(intType j, intType k)
    {
        for (intType i = 0; i != m_ni - 1; i++) { // Forward sweep
            UpdateAndRelax(m_blockSolverEast, i, j, k);
        }

        for (intType i = m_ni - 1; i != 0; i--) { // Backward sweep
            UpdateAndRelax(m_blockSolverWest, i, j, k);
        }

    }

    void UpdateState3D()
    {
        m_blockSolverEast->UpdateGlobalConstants();
        m_blockSolverWest->UpdateGlobalConstants();
    }


    // For 2D simulations
    void Sweep2D(intType j, intType k)
    {
        UpdateAndRelax(m_blockSolverCenter, 0, j, k);
    }

    void UpdateState2D()
    {
        m_blockSolverCenter->UpdateGlobalConstants();
    }

    template <TC Ustag>
    void UpdateAndRelax(std::unique_ptr<TriadSolver<Ustag, Vstag, Wstag>> &blockSolver, intType i, intType j, intType k)
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Fields::ENUMDATA;
        using enum Axis::ENUMDATA;

        EnumFor<Fields>( [&] (Fields::ENUMDATA f) { m_iS[f] = i; } ); // Set iterating coefficient
        m_iS[U] += CoeffIndex[Ustag];      // U momentum is staggered

        EnumFor<Fields>([&](Fields::ENUMDATA f)
                        { m_oldBlock[f] = m_fields[f](G(m_iS[f], m_jS[f], m_kS[f])); }); // Set old block values

        blockSolver->UpdateTriad(i, j, k);

        EnumFor<Fields>( [&] (Fields::ENUMDATA f) {
            auto &fieldBlock = m_fields[f](G(m_iS[f], m_jS[f], m_kS[f]));
            m_delta[f] = m_relaxation[f] * (fieldBlock - m_oldBlock[f]); // Relaxed change in solution
            fieldBlock = m_oldBlock[f] + m_delta[f];                     // Apply relaxation
            m_residuals[f] += abs(m_delta[f]);                           // Add to residual count
        } );
    }
};





template <TransportCoefficients::ENUMDATA Wstag>
class PlaneSolver
{
    using TC = TransportCoefficients::ENUMDATA;
    using A = Axis::ENUMDATA;

    // Staggering must be valid
    static_assert( (Wstag == TC::t) || (Wstag == TC::b) || (Wstag == TC::p), "Invalid W momentum staggering" );

public:
    PlaneSolver( ArrayAllocator<Fields, array3D> &fields,
                 const ArrayAllocator<Fields, array3D> &fieldsOld,
                 const FVCoefficients &fvCoeffs,
                 const InputData::PlaneSolverSettings &planeSolverSettings) : 
                    m_fields( fields ),
                    m_fieldsOld( fieldsOld ),
                    m_maxIterations( planeSolverSettings.maxIterations ),
                    m_maxResiduals( planeSolverSettings.maxResiduals ),
                    m_relaxation( planeSolverSettings.relaxation ),

                    m_oldLine( array1D( m_fields[F::U].dimension(A::X) ) ),
                    m_delta( array1D( m_fields[F::U].dimension(A::X) ) ),

                    m_ni( fvCoeffs.nCells(A::X) ),
                    m_nj( fvCoeffs.nCells(A::Y) )
    {
        if (m_nj == 1) {

            m_lineSolverCenter = std::make_unique<LineSolver<TC::p, Wstag>>(fields, fieldsOld, fvCoeffs, planeSolverSettings.lineSolverSettings);
            SolutionUpdater = &PlaneSolver::Sweep2D;
            StateUpdater = &PlaneSolver::UpdateState2D;

        } else {

            m_lineSolverNorth = std::make_unique<LineSolver<TC::n, Wstag>>(fields, fieldsOld, fvCoeffs, planeSolverSettings.lineSolverSettings);
            m_lineSolverSouth = std::make_unique<LineSolver<TC::s, Wstag>>(fields, fieldsOld, fvCoeffs, planeSolverSettings.lineSolverSettings);
            SolutionUpdater = &PlaneSolver::Sweep3D;
            StateUpdater = &PlaneSolver::UpdateState3D;

        }
    }


    void SolvePlane(const intType k)
    {
        using enum Fields::ENUMDATA;
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        // Staggered indexing
        EnumFor<Fields>( [&] (Fields::ENUMDATA f) { m_kS[f] = k; } );
        m_kS[W] += CoeffIndex[Wstag];

        // Solver loop
        for ( intType nIterations = 1; nIterations <= m_maxIterations; nIterations++ )
        {
            EnumFor<Fields>( [&] (Fields::ENUMDATA field) { m_residuals[field] = 0.0f; } );

            // Update Line
            (this->*SolutionUpdater)(k);

            // Normalise residuals
            EnumFor<Fields>( [&] (Fields::ENUMDATA field) { m_residuals[field] /= static_cast<floatType>(m_ni * m_nj); } );
            RelativeResidual(m_residuals, m_residualsInitialInv, nIterations);

            // Check residual tolerence
            if ( MetResidualTolerence(m_residuals, m_maxResiduals) ) {
                break;
            }
        }
    }


    void UpdateState()
    {
        (this->*StateUpdater)();
    }


private:

    ArrayAllocator<Fields, array3D> &m_fields;
    const ArrayAllocator<Fields, array3D> &m_fieldsOld;
    const intType m_maxIterations;
    const EnumVector<Fields, floatType> m_maxResiduals;
    const EnumVector<Fields, floatType> m_relaxation;

    std::unique_ptr< LineSolver<TC::n, Wstag> > m_lineSolverNorth;
    std::unique_ptr< LineSolver<TC::s, Wstag> > m_lineSolverSouth;
    std::unique_ptr< LineSolver<TC::p, Wstag> > m_lineSolverCenter;

    void (PlaneSolver::*SolutionUpdater)(intType);
    void (PlaneSolver::*StateUpdater)();

    EnumVector<Fields, array1D> m_oldLine, m_delta;
    EnumVector<Fields, floatType> m_residuals, m_residualsInitialInv;
    EnumVector<Fields, intType> m_jS, m_kS;

    intType m_ni, m_nj;

    // For 3D simulations
    void Sweep3D(intType k)
    {
        for (intType j = 0; j != m_nj - 1; j++) { // Forward sweep
            UpdateAndRelax(m_lineSolverNorth, j, k);
        }

        for (intType j = m_nj - 1; j != 0; j--) { // Backward sweep
            UpdateAndRelax(m_lineSolverSouth, j, k);
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
        UpdateAndRelax(m_lineSolverCenter, 0, k);
    }

    void UpdateState2D()
    {
        m_lineSolverCenter->UpdateState();
    }


    template <TC Vstag>
    void UpdateAndRelax(std::unique_ptr<LineSolver<Vstag, Wstag>> &lineSolver, intType j, intType k)
    {
        using enum Fields::ENUMDATA;
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        EnumFor<Fields>( [&] (Fields::ENUMDATA f) { m_jS[f] = j; }); // Set iterating coefficient
        m_jS[V] += CoeffIndex[Vstag];      // V momentum is staggered

        EnumFor<Fields>( [&] (Fields::ENUMDATA f) { m_oldLine[f] = m_fields[f].chip( G(m_kS[f]), Z ).chip( G(m_jS[f]), Y ); } ); // Set old line

        lineSolver->SolveLine(j, k);

        EnumFor<Fields>( [&] (Fields::ENUMDATA f) {
            auto fieldLine = m_fields[f].chip( G(m_kS[f]), Z ).chip( G(m_jS[f]), Y );
            m_delta[f] = m_delta[f].constant( m_relaxation[f] ) * ( fieldLine - m_oldLine[f]) ; // Relaxed change in line
            fieldLine = m_oldLine[f] + m_delta[f];                                          // Relax
            m_residuals[f] += static_cast<array0D>( m_delta[f].abs().sum() )(0);              // Add to residual count
        } );
    }
};





class LinearSolver
{
    using TC = TransportCoefficients::ENUMDATA;
    using A = Axis::ENUMDATA;

public:
    LinearSolver( ArrayAllocator<Fields, array3D> &fields,
                  const ArrayAllocator<Fields, array3D> &fieldsOld,
                  const FVCoefficients &fvCoeffs, 
                  const InputData::LinearSolverSettings &linearSolverSettings) : 
                    m_fields( fields ),
                    m_maxIterations( linearSolverSettings.maxIterations ),
                    m_maxResiduals( linearSolverSettings.maxResiduals ),
                    m_relaxation( linearSolverSettings.relaxation ),

                    m_delta( array2D( m_fields[F::U].dimension(A::X), m_fields[F::U].dimension(A::Y) ) ),
                    m_oldPlane( array2D( m_fields[F::U].dimension(A::X), m_fields[F::U].dimension(A::Y) ) ),

                    m_ni( fvCoeffs.nCells(A::X) ),
                    m_nj( fvCoeffs.nCells(A::Y) ),
                    m_nk( fvCoeffs.nCells(A::Z) )
    {
        if (m_nk == 1) {

            m_planeSolverCenter = std::make_unique<PlaneSolver<TC::p>>(fields, fieldsOld, fvCoeffs, linearSolverSettings.planeSolverSettings);
            SolutionUpdater = &LinearSolver::Sweep2D;
            StateUpdater = &LinearSolver::UpdateState2D;

        } else {

            m_planeSolverTop = std::make_unique<PlaneSolver<TC::t>>(fields, fieldsOld, fvCoeffs, linearSolverSettings.planeSolverSettings);
            m_planeSolverBottom = std::make_unique<PlaneSolver<TC::b>>(fields, fieldsOld, fvCoeffs, linearSolverSettings.planeSolverSettings);
            SolutionUpdater = &LinearSolver::Sweep3D;
            StateUpdater = &LinearSolver::UpdateState3D;

        }
    }


    void Solve()
    {
        using enum Fields::ENUMDATA;
        using enum Axis::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;

        for ( intType nIterations = 1; nIterations <= m_maxIterations; nIterations++ )
        {

            // Reset residuals
            EnumFor<Fields>( [&] (Fields::ENUMDATA field) { m_residuals[field] = 0.0f; });

            // Update plane
            (this->*SolutionUpdater)();

            // Normalise residuals
            EnumFor<Fields>( [&] (Fields::ENUMDATA field) { m_residuals[field] /= static_cast<floatType>(m_ni * m_nj * m_nk); });
            RelativeResidual(m_residuals, m_residualsInitialInv, nIterations);

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

    ArrayAllocator<Fields, array3D> &m_fields;
    const intType m_maxIterations;
    const EnumVector<Fields, floatType> m_maxResiduals;
    const EnumVector<Fields, floatType> m_relaxation;

    std::unique_ptr<PlaneSolver<TC::t>> m_planeSolverTop;
    std::unique_ptr<PlaneSolver<TC::b>> m_planeSolverBottom;
    std::unique_ptr<PlaneSolver<TC::p>> m_planeSolverCenter;

    void (LinearSolver::*SolutionUpdater)(void);
    void (LinearSolver::*StateUpdater)(void);

    EnumVector<Fields, array2D> m_delta, m_oldPlane;
    EnumVector<Fields, floatType> m_residuals, m_residualsInitialInv;
    EnumVector<Fields, intType> m_kS;

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
    void UpdateAndRelax(std::unique_ptr<PlaneSolver<Wstag>> &planeSolver, intType k)
    {
        using enum Fields::ENUMDATA;
        using enum TransportCoefficients::ENUMDATA;
        using enum Axis::ENUMDATA;

        EnumFor<Fields>( [&] (Fields::ENUMDATA f) { m_kS[f] = k; });  // Set iterating coefficient
        m_kS[W] += CoeffIndex[Wstag];                                 // W momentum is staggered

        EnumFor<Fields>( [&] (Fields::ENUMDATA f) { m_oldPlane[f] = m_fields[f].chip( G(m_kS[f]), Z ); } ); // Set old plane

        planeSolver->SolvePlane(k);

        EnumFor<Fields>( [&] (Fields::ENUMDATA f) {
            auto fieldPlane = m_fields[f].chip( G(m_kS[f]), Z );
            m_delta[f] = m_delta[f].constant( m_relaxation[f] ) * (fieldPlane - m_oldPlane[f]); // Relaxed change in plane
            fieldPlane = m_oldPlane[f] + m_delta[f];                                            // Relax
            m_residuals[f] += static_cast<array0D>( m_delta[f].abs().sum() )(0);                // Add to residual count
        } );
    }
};





void SweepSolve( CellFields &fields,
                 const Mesh &mesh,
                 const InputData &inputData,
                 const AxisTransformationMap &axisTransformation)
{
    using enum Axis::ENUMDATA;
    using enum Fields::ENUMDATA;

    // Extract from input data
    const InputData::LinearSolverSettings linearSolverSettings = inputData.linearSolverSettings;
    const intType maxOuterIterations = inputData.schemes.maxOuterIterations;
    const EnumVector<Fields, floatType> maxOuterResiduals = inputData.schemes.maxOuterResiduals;

    // Initialise
    ArrayAllocator<Fields, array3D> faceVelocities = InitialiseFaceVelocities(mesh, fields, inputData);
    ArrayAllocator<Fields, array3D> fieldsOld(fields);
    FVCoefficients fvCoeffs = InitialiseFVCoefficients(mesh, faceVelocities, inputData);

    EnumVector<Fields, floatType> residualsOuter, residualsOuterInitialInv;
    floatType massFluxResidual;

    ResidualLogFile residualsLogFile("convergence_history.csv", axisTransformation);
    ConsoleLog consoleLog( axisTransformation );

    // Instantiate linear solver, this holds references to the fields
    LinearSolver linearSolver(fields, fieldsOld, fvCoeffs, linearSolverSettings);

    // Outer iterations
    for ( intType nOuterIterations = 1; nOuterIterations <= maxOuterIterations; nOuterIterations++ )
    {
        linearSolver.UpdateState();
        linearSolver.Solve();

        UpdateFaceVelocities(faceVelocities, mesh, fields, inputData);
        UpdateFVCoefficients(fvCoeffs, mesh, faceVelocities, inputData);
        
        residualsOuter = L1ArrayDiff(fields, fieldsOld);
        massFluxResidual = BoundaryMassFluxResidual(faceVelocities, mesh);
        RelativeResidual(residualsOuter, residualsOuterInitialInv, nOuterIterations);

        fieldsOld = fields;

        consoleLog.WriteResiduals(residualsOuter, massFluxResidual, nOuterIterations);
        residualsLogFile.WriteData(residualsOuter, massFluxResidual, nOuterIterations);

        if (MetResidualTolerence(residualsOuter, maxOuterResiduals)) {
            std::cout << "*** OUTER ITERATIONS CONVERGED ***"
                        << "\n\n";
            break;
        }
        
    }
}

} // end namespace CFD