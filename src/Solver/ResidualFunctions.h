#ifndef RESIDUAL_FUNCTIONS
#define RESIDUAL_FUNCTIONS

#include "../Core/Macros.h"
#include "../Core/Types.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../DerivedQuantities/DerivedQuantities.h"


#include <cmath>

namespace CFD
{

using namespace FVT;

// Calculated as the L1 norm of the difference between two arrays.
[[ maybe_unused ]]
inline FieldData<floatType> L1DiffResiduals( const FieldData<Tensor3D> &fields1,
                                             const FieldData<Tensor3D> &fields2 )
{
    FieldData<floatType> result;
    ForAllFieldData( [&] (intType i) { 

        auto fieldDiff = fields2[i] - fields1[i];  // auto lazily evaluates
        result[i] = static_cast<Tensor0D>( fieldDiff.abs().mean() )(0); 

    });

    return result;
}



// Cell residual from stencil for momentum equation
__attribute__((always_inline)) 
inline floatType CellMomentumResidual_x( const FieldData<Tensor3D> &fields,
                                         const FVCoefficients &fvCoeffs,
                                         const intType ig,
                                         const intType jg,
                                         const intType kg )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    auto & xMomCoeffs = fvCoeffs.Mom.coeffs[X];

    floatType residual = - xMomCoeffs.AU[p](ig, jg, kg) * fields.U[X]( ig  , jg  , kg  ) 
                         - xMomCoeffs.AU[n](ig, jg, kg) * fields.U[X]( ig  , jg+1, kg  ) 
                         - xMomCoeffs.AU[e](ig, jg, kg) * fields.U[X]( ig+1, jg  , kg  ) 
                         - xMomCoeffs.AU[s](ig, jg, kg) * fields.U[X]( ig  , jg-1, kg  ) 
                         - xMomCoeffs.AU[w](ig, jg, kg) * fields.U[X]( ig-1, jg  , kg  ) 
                         - xMomCoeffs.AU[t](ig, jg, kg) * fields.U[X]( ig  , jg  , kg+1) 
                         - xMomCoeffs.AU[b](ig, jg, kg) * fields.U[X]( ig  , jg  , kg-1) 

                         - xMomCoeffs.AP[e](ig) * fields.P( ig+1, jg  , kg  )
                         - xMomCoeffs.AP[p](ig) * fields.P( ig  , jg  , kg  )
                         - xMomCoeffs.AP[w](ig) * fields.P( ig-1, jg  , kg  )

                         - xMomCoeffs.B(ig, jg, kg) 
                         + xMomCoeffs.F(ig, jg, kg);

    return residual;
}



__attribute__((always_inline)) 
inline floatType CellMomentumResidual_y( const FieldData<Tensor3D> &fields,
                                         const FVCoefficients &fvCoeffs,
                                         const intType ig,
                                         const intType jg,
                                         const intType kg )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    auto & yMomCoeffs = fvCoeffs.Mom.coeffs[Y];

    floatType residual = - yMomCoeffs.AU[p](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg  ) 
                         - yMomCoeffs.AU[n](ig, jg, kg) * fields.U[Y]( ig  , jg+1, kg  ) 
                         - yMomCoeffs.AU[e](ig, jg, kg) * fields.U[Y]( ig+1, jg  , kg  ) 
                         - yMomCoeffs.AU[s](ig, jg, kg) * fields.U[Y]( ig  , jg-1, kg  ) 
                         - yMomCoeffs.AU[w](ig, jg, kg) * fields.U[Y]( ig-1, jg  , kg  ) 
                         - yMomCoeffs.AU[t](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg+1) 
                         - yMomCoeffs.AU[b](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg-1) 

                         - yMomCoeffs.AP[n](jg) * fields.P( ig  , jg+1, kg  )
                         - yMomCoeffs.AP[p](jg) * fields.P( ig  , jg  , kg  )
                         - yMomCoeffs.AP[s](jg) * fields.P( ig  , jg-1, kg  )

                         - yMomCoeffs.B(ig, jg, kg)
                         + yMomCoeffs.F(ig, jg, kg);

    return residual;
}



__attribute__((always_inline)) 
inline floatType CellMomentumResidual_z( const FieldData<Tensor3D> &fields,
                                         const FVCoefficients &fvCoeffs,
                                         const intType ig,
                                         const intType jg,
                                         const intType kg )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    auto & zMomCoeffs = fvCoeffs.Mom.coeffs[Z];

    floatType residual = - zMomCoeffs.AU[p](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg  ) 
                         - zMomCoeffs.AU[n](ig, jg, kg) * fields.U[Z]( ig  , jg+1, kg  ) 
                         - zMomCoeffs.AU[e](ig, jg, kg) * fields.U[Z]( ig+1, jg  , kg  ) 
                         - zMomCoeffs.AU[s](ig, jg, kg) * fields.U[Z]( ig  , jg-1, kg  ) 
                         - zMomCoeffs.AU[w](ig, jg, kg) * fields.U[Z]( ig-1, jg  , kg  ) 
                         - zMomCoeffs.AU[t](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg+1) 
                         - zMomCoeffs.AU[b](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg-1) 

                         - zMomCoeffs.AP[t](kg) * fields.P( ig  , jg  , kg+1)
                         - zMomCoeffs.AP[p](kg) * fields.P( ig  , jg  , kg  )
                         - zMomCoeffs.AP[b](kg) * fields.P( ig  , jg  , kg-1)

                         - zMomCoeffs.B(ig, jg, kg)
                         + zMomCoeffs.F(ig, jg, kg);

    return residual;
}



template< MomentumInterpolation MI > __attribute__((always_inline)) 
inline floatType CellContinuityResidual( const FieldData<Tensor3D> &fields,
                                         const FVCoefficients &fvCoeffs,
                                         const intType ig,
                                         const intType jg,
                                         const intType kg )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType pressureWideStencil = 0.0f;
    auto & contCoeffs = fvCoeffs.Cont.coeffs;
    if constexpr ( MI == MomentumInterpolation::Implicit ) {
        pressureWideStencil = contCoeffs.AP[nn](ig, jg, kg) * fields.P( ig  , jg+2, kg  ) 
                            + contCoeffs.AP[ee](ig, jg, kg) * fields.P( ig+2, jg  , kg  ) 
                            + contCoeffs.AP[ss](ig, jg, kg) * fields.P( ig  , jg-2, kg  ) 
                            + contCoeffs.AP[ww](ig, jg, kg) * fields.P( ig-2, jg  , kg  ) 
                            + contCoeffs.AP[tt](ig, jg, kg) * fields.P( ig  , jg  , kg+2) 
                            + contCoeffs.AP[bb](ig, jg, kg) * fields.P( ig  , jg  , kg-2);
    }
    floatType residual = - contCoeffs.AU[X][e](ig) * fields.U[X]( ig+1, jg  , kg  )
                         - contCoeffs.AU[X][p](ig) * fields.U[X]( ig  , jg  , kg  )
                         - contCoeffs.AU[X][w](ig) * fields.U[X]( ig-1, jg  , kg  )

                         - contCoeffs.AU[Y][n](jg) * fields.U[Y]( ig  , jg+1, kg  )
                         - contCoeffs.AU[Y][p](jg) * fields.U[Y]( ig  , jg  , kg  )
                         - contCoeffs.AU[Y][s](jg) * fields.U[Y]( ig  , jg-1, kg  )

                         - contCoeffs.AU[Z][t](kg) * fields.U[Z]( ig  , jg  , kg+1)
                         - contCoeffs.AU[Z][p](kg) * fields.U[Z]( ig  , jg  , kg  )
                         - contCoeffs.AU[Z][b](kg) * fields.U[Z]( ig  , jg  , kg-1)

                         - contCoeffs.AP[p](ig, jg, kg) * fields.P( ig  , jg  , kg  )
                         - contCoeffs.AP[n](ig, jg, kg) * fields.P( ig  , jg+1, kg  ) 
                         - contCoeffs.AP[e](ig, jg, kg) * fields.P( ig+1, jg  , kg  ) 
                         - contCoeffs.AP[s](ig, jg, kg) * fields.P( ig  , jg-1, kg  ) 
                         - contCoeffs.AP[w](ig, jg, kg) * fields.P( ig-1, jg  , kg  ) 
                         - contCoeffs.AP[t](ig, jg, kg) * fields.P( ig  , jg  , kg+1) 
                         - contCoeffs.AP[b](ig, jg, kg) * fields.P( ig  , jg  , kg-1)

                         - pressureWideStencil
                         
                         - contCoeffs.B(ig, jg, kg)
                         + contCoeffs.F(ig, jg, kg);

    return residual;
}



// Calculate the absolute residual of each equation from the finite volume stencil
template< MomentumInterpolation MI > [[ maybe_unused ]]
inline FieldData<floatType> ScaledL1NormResiduals( const FieldData<Tensor3D> &fields,
                                                   const FVCoefficients &fvCoeffs,
                                                   const Tensor3D &mask )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    FieldData<floatType> residuals{0};
    FieldData<floatType> scalingFactor{0};

    // Temporaries to allow vector reductions
    floatType residualU{0}     , residualV{0}     , residualW{0}     , residualP{0},
              scalingFactorU{0}, scalingFactorV{0}, scalingFactorW{0};

    for ( intType k = 0; k != fvCoeffs.nCells[Z]; k++ ) {
        for ( intType j = 0; j != fvCoeffs.nCells[Y]; j++ ) {
            for ( intType i = 0; i != fvCoeffs.nCells[X]; i++ ) {

                intType ig{ G(i) }, jg{ G(j) }, kg{ G(k) };

                // U momentum
                residualU      += mask(ig, jg, kg) * abs( CellMomentumResidual_x(fields, fvCoeffs, ig, jg, kg) );
                scalingFactorU += mask(ig, jg, kg) * abs( fvCoeffs.Mom.coeffs[X].AU[p](ig, jg, kg) * fields.U[X]( ig  , jg  , kg  ) );


                // V momentum
                residualV      += mask(ig, jg, kg) * abs( CellMomentumResidual_y(fields, fvCoeffs, ig, jg, kg) );
                scalingFactorV += mask(ig, jg, kg) * abs( fvCoeffs.Mom.coeffs[Y].AU[p](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg  ) );


                // W momentm
                residualW      += mask(ig, jg, kg) * abs( CellMomentumResidual_z(fields, fvCoeffs, ig, jg, kg) );             
                scalingFactorW += mask(ig, jg, kg) * abs( fvCoeffs.Mom.coeffs[Z].AU[p](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg  ) );


                // Continuity 
                residualP   += mask(ig, jg, kg) * abs( CellContinuityResidual<MI>(fields, fvCoeffs, ig, jg, kg) );

            }
        }
    }

    // U momentum
    residuals.U[X] = residualU;
    scalingFactor.U[X] = scalingFactorU;

    // V momentum
    residuals.U[Y] = residualV;
    scalingFactor.U[Y] = scalingFactorV;

    // W momentm
    residuals.U[Z] = residualW;
    scalingFactor.U[Z] = scalingFactorW;

    // Continuity 
    residuals.P = residualP;


    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        if ( scalingFactor.U[axis] != 0.0f ) { // Division by zero
            residuals.U[axis] /= scalingFactor.U[axis];
        }
    } );
    residuals.P /= static_cast<floatType>( fvCoeffs.nCells(X) * fvCoeffs.nCells(Y) * fvCoeffs.nCells(Z) );

    return residuals;
}



// Calculate the local discrete equation residual field
template< MomentumInterpolation MI >
inline FieldData<Tensor3D> ResidualsField( const FieldData<Tensor3D> &fields,
                                           const FVCoefficients &fvCoeffs,
                                           const Tensor3D &mask )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    FieldData<Tensor3D> residuals( Tensor3D( fvCoeffs.nCells[X] + 2*CFD::nGhost, 
                                             fvCoeffs.nCells[Y] + 2*CFD::nGhost, 
                                             fvCoeffs.nCells[Z] + 2*CFD::nGhost ).setZero() );

    for ( intType k = 0; k != fvCoeffs.nCells[Z]; k++ ) {
        for ( intType j = 0; j != fvCoeffs.nCells[Y]; j++ ) {
            for ( intType i = 0; i != fvCoeffs.nCells[X]; i++ ) {

                intType ig{ G(i) }, jg{ G(j) }, kg{ G(k) };

                // U momentum
                residuals.U[X](ig, jg, kg) = mask(ig, jg, kg) * CellMomentumResidual_x(fields, fvCoeffs, ig, jg, kg);

                // V momentum
                residuals.U[Y](ig, jg, kg) = mask(ig, jg, kg) * CellMomentumResidual_y(fields, fvCoeffs, ig, jg, kg);

                // W momentm
                residuals.U[Z](ig, jg, kg) = mask(ig, jg, kg) * CellMomentumResidual_z(fields, fvCoeffs, ig, jg, kg);

                // Continuity 
                residuals.P(ig, jg, kg)    = mask(ig, jg, kg) * CellContinuityResidual<MI>(fields, fvCoeffs, ig, jg, kg);

            }
        }
    }

    return residuals;
}



// Calculate global mass flux residual at the domain boundary
[[ maybe_unused ]]
inline floatType BoundaryMassFluxResidual( const EnumVector<Axis, Tensor3D> &faceFluxes,
                                           const Mesh &mesh )
{
    floatType massFluxResidual = 0.0f;

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        // Positive face, area normal is in positive direction
        auto faceFluxesPositive = faceFluxes[axis].chip( mesh.nCells(axis), axis ) * mesh.cellFaceAreas[axis];
        massFluxResidual += static_cast<Tensor0D>( faceFluxesPositive.sum() )(0);

        // Negative face, area normal is in negative direction
        auto faceFluxesNegative = - faceFluxes[axis].chip( 0, axis ) * mesh.cellFaceAreas[axis];
        massFluxResidual += static_cast<Tensor0D>( faceFluxesNegative.sum() )(0);

    });

    return massFluxResidual;
}

// Set the normalisation factor for the residuals
[[ maybe_unused ]]
inline void SetResidualsNormalisationFactor( FieldData<floatType> &residualsScaleFactor,
                                             const FieldData<floatType> &residuals ) {
    ForAllFieldData( [&] (intType f) {
        residualsScaleFactor[f] = 1.0f;
    } );
    if ( residuals.P != 0 ) { // Division by zero
        residualsScaleFactor.P = 1.0f / residuals.P;
    } 
}



// Apply normalisation to residuals
[[ maybe_unused ]]
inline void NormaliseResiduals( FieldData<floatType> &residuals,
                                FieldData<floatType> &residualsScaleFactor )
{
    ForAllFieldData( [&] (intType i) { residuals[i] *= residualsScaleFactor[i]; } );
}



// Check if residual tolerence is met
inline bool MetResidualTolerence( const FieldData<floatType> &residuals,
                                  const FieldData<floatType> &residualsTarget )
{
    for ( intType i = 0; i != FieldData<floatType>::nData; i++ ) {  // Can't use ForAllFieldData since returning inside loop

        if ( residuals[i] > residualsTarget[i] )
            return false;

    }
    return true;
}



// Check if any residuals have diverged
inline bool ResidualsDiverged( const FieldData<floatType> &residuals )
{
    for ( intType i = 0; i != FieldData<floatType>::nData; i++ ) {
        if ( !std::isfinite( residuals[i] ) ) {
            return true;
        }   
    }

    return false;
}

}   // end namespace CFD

#endif // RESIDUAL_FUNCTIONS