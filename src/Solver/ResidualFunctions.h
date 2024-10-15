#ifndef RESIDUAL_FUNCTIONS
#define RESIDUAL_FUNCTIONS

#include "../Types.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Tools/FieldProbe.h"
#include "../Tools/FVTools.h"
#include "../Macros.h"

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
template< Linearisation LI > __attribute__((always_inline)) 
inline floatType CellMomentumResidual_x( const FieldData<Tensor3D> &fields,
                                         const FVCoefficients &fvCoeffs,
                                         const intType ig,
                                         const intType jg,
                                         const intType kg )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType newtonStencilX = 0.0f;
    if constexpr ( LI == Linearisation::Newton ) {
        newtonStencilX = fvCoeffs.Mom[X].AU[Y][n](ig, jg, kg) * fields.U[Y]( ig  , jg+1, kg  )
                       + fvCoeffs.Mom[X].AU[Y][p](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg  )
                       + fvCoeffs.Mom[X].AU[Y][s](ig, jg, kg) * fields.U[Y]( ig  , jg-1, kg  )

                       + fvCoeffs.Mom[X].AU[Z][t](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg+1)
                       + fvCoeffs.Mom[X].AU[Z][p](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg  )
                       + fvCoeffs.Mom[X].AU[Z][b](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg-1);
    }
    floatType residual = - fvCoeffs.Mom[X].AU[X][p](ig, jg, kg) * fields.U[X]( ig  , jg  , kg  ) 
                         - fvCoeffs.Mom[X].AU[X][n](ig, jg, kg) * fields.U[X]( ig  , jg+1, kg  ) 
                         - fvCoeffs.Mom[X].AU[X][e](ig, jg, kg) * fields.U[X]( ig+1, jg  , kg  ) 
                         - fvCoeffs.Mom[X].AU[X][s](ig, jg, kg) * fields.U[X]( ig  , jg-1, kg  ) 
                         - fvCoeffs.Mom[X].AU[X][w](ig, jg, kg) * fields.U[X]( ig-1, jg  , kg  ) 
                         - fvCoeffs.Mom[X].AU[X][t](ig, jg, kg) * fields.U[X]( ig  , jg  , kg+1) 
                         - fvCoeffs.Mom[X].AU[X][b](ig, jg, kg) * fields.U[X]( ig  , jg  , kg-1) 

                         - fvCoeffs.Mom[X].AP[e](ig) * fields.P( ig+1, jg  , kg  )
                         - fvCoeffs.Mom[X].AP[p](ig) * fields.P( ig  , jg  , kg  )
                         - fvCoeffs.Mom[X].AP[w](ig) * fields.P( ig-1, jg  , kg  )

                         - newtonStencilX

                         - fvCoeffs.Mom[X].B(ig, jg, kg) 
                         + fvCoeffs.Mom[X].F(ig, jg, kg);

    return residual;
}



template< Linearisation LI > __attribute__((always_inline)) 
inline floatType CellMomentumResidual_y( const FieldData<Tensor3D> &fields,
                                         const FVCoefficients &fvCoeffs,
                                         const intType ig,
                                         const intType jg,
                                         const intType kg )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType newtonStencilY = 0.0f;
    if constexpr ( LI == Linearisation::Newton ) {
        newtonStencilY = fvCoeffs.Mom[Y].AU[X][e](ig, jg, kg) * fields.U[X]( ig+1, jg  , kg  )
                       + fvCoeffs.Mom[Y].AU[X][p](ig, jg, kg) * fields.U[X]( ig  , jg  , kg  )
                       + fvCoeffs.Mom[Y].AU[X][w](ig, jg, kg) * fields.U[X]( ig-1, jg  , kg  )
        
                       + fvCoeffs.Mom[Y].AU[Z][t](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg+1)
                       + fvCoeffs.Mom[Y].AU[Z][p](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg  )
                       + fvCoeffs.Mom[Y].AU[Z][b](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg-1);      
    }
    floatType residual = - fvCoeffs.Mom[Y].AU[Y][p](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg  ) 
                         - fvCoeffs.Mom[Y].AU[Y][n](ig, jg, kg) * fields.U[Y]( ig  , jg+1, kg  ) 
                         - fvCoeffs.Mom[Y].AU[Y][e](ig, jg, kg) * fields.U[Y]( ig+1, jg  , kg  ) 
                         - fvCoeffs.Mom[Y].AU[Y][s](ig, jg, kg) * fields.U[Y]( ig  , jg-1, kg  ) 
                         - fvCoeffs.Mom[Y].AU[Y][w](ig, jg, kg) * fields.U[Y]( ig-1, jg  , kg  ) 
                         - fvCoeffs.Mom[Y].AU[Y][t](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg+1) 
                         - fvCoeffs.Mom[Y].AU[Y][b](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg-1) 

                         - fvCoeffs.Mom[Y].AP[n](jg) * fields.P( ig  , jg+1, kg  )
                         - fvCoeffs.Mom[Y].AP[p](jg) * fields.P( ig  , jg  , kg  )
                         - fvCoeffs.Mom[Y].AP[s](jg) * fields.P( ig  , jg-1, kg  )

                         - newtonStencilY

                         - fvCoeffs.Mom[Y].B(ig, jg, kg)
                         + fvCoeffs.Mom[Y].F(ig, jg, kg);

    return residual;
}



template< Linearisation LI > __attribute__((always_inline)) 
inline floatType CellMomentumResidual_z( const FieldData<Tensor3D> &fields,
                                         const FVCoefficients &fvCoeffs,
                                         const intType ig,
                                         const intType jg,
                                         const intType kg )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    floatType newtonStencilZ = 0.0f;
    if constexpr ( LI == Linearisation::Newton ) {
        newtonStencilZ =  fvCoeffs.Mom[Z].AU[X][e](ig, jg, kg) * fields.U[X]( ig+1, jg  , kg  )
                        + fvCoeffs.Mom[Z].AU[X][p](ig, jg, kg) * fields.U[X]( ig  , jg  , kg  )
                        + fvCoeffs.Mom[Z].AU[X][w](ig, jg, kg) * fields.U[X]( ig-1, jg  , kg  )
        
                        + fvCoeffs.Mom[Z].AU[Y][n](ig, jg, kg) * fields.U[Y]( ig  , jg+1, kg  )
                        + fvCoeffs.Mom[Z].AU[Y][p](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg  )
                        + fvCoeffs.Mom[Z].AU[Y][s](ig, jg, kg) * fields.U[Y]( ig  , jg-1, kg  );    
    }
    floatType residual = - fvCoeffs.Mom[Z].AU[Z][p](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg  ) 
                         - fvCoeffs.Mom[Z].AU[Z][n](ig, jg, kg) * fields.U[Z]( ig  , jg+1, kg  ) 
                         - fvCoeffs.Mom[Z].AU[Z][e](ig, jg, kg) * fields.U[Z]( ig+1, jg  , kg  ) 
                         - fvCoeffs.Mom[Z].AU[Z][s](ig, jg, kg) * fields.U[Z]( ig  , jg-1, kg  ) 
                         - fvCoeffs.Mom[Z].AU[Z][w](ig, jg, kg) * fields.U[Z]( ig-1, jg  , kg  ) 
                         - fvCoeffs.Mom[Z].AU[Z][t](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg+1) 
                         - fvCoeffs.Mom[Z].AU[Z][b](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg-1) 

                         - fvCoeffs.Mom[Z].AP[t](kg) * fields.P( ig  , jg  , kg+1)
                         - fvCoeffs.Mom[Z].AP[p](kg) * fields.P( ig  , jg  , kg  )
                         - fvCoeffs.Mom[Z].AP[b](kg) * fields.P( ig  , jg  , kg-1)

                         - newtonStencilZ

                         - fvCoeffs.Mom[Z].B(ig, jg, kg)
                         + fvCoeffs.Mom[Z].F(ig, jg, kg);

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
    if constexpr ( MI == MomentumInterpolation::Implicit ) {
        pressureWideStencil = fvCoeffs.Cont.AP[nn](ig, jg, kg) * fields.P( ig  , jg+2, kg  ) 
                            + fvCoeffs.Cont.AP[ee](ig, jg, kg) * fields.P( ig+2, jg  , kg  ) 
                            + fvCoeffs.Cont.AP[ss](ig, jg, kg) * fields.P( ig  , jg-2, kg  ) 
                            + fvCoeffs.Cont.AP[ww](ig, jg, kg) * fields.P( ig-2, jg  , kg  ) 
                            + fvCoeffs.Cont.AP[tt](ig, jg, kg) * fields.P( ig  , jg  , kg+2) 
                            + fvCoeffs.Cont.AP[bb](ig, jg, kg) * fields.P( ig  , jg  , kg-2);
    }
    floatType residual = - fvCoeffs.Cont.AU[X][e](ig) * fields.U[X]( ig+1, jg  , kg  )
                         - fvCoeffs.Cont.AU[X][p](ig) * fields.U[X]( ig  , jg  , kg  )
                         - fvCoeffs.Cont.AU[X][w](ig) * fields.U[X]( ig-1, jg  , kg  )

                         - fvCoeffs.Cont.AU[Y][n](jg) * fields.U[Y]( ig  , jg+1, kg  )
                         - fvCoeffs.Cont.AU[Y][p](jg) * fields.U[Y]( ig  , jg  , kg  )
                         - fvCoeffs.Cont.AU[Y][s](jg) * fields.U[Y]( ig  , jg-1, kg  )

                         - fvCoeffs.Cont.AU[Z][t](kg) * fields.U[Z]( ig  , jg  , kg+1)
                         - fvCoeffs.Cont.AU[Z][p](kg) * fields.U[Z]( ig  , jg  , kg  )
                         - fvCoeffs.Cont.AU[Z][b](kg) * fields.U[Z]( ig  , jg  , kg-1)

                         - fvCoeffs.Cont.AP[p](ig, jg, kg) * fields.P( ig  , jg  , kg  )
                         - fvCoeffs.Cont.AP[n](ig, jg, kg) * fields.P( ig  , jg+1, kg  ) 
                         - fvCoeffs.Cont.AP[e](ig, jg, kg) * fields.P( ig+1, jg  , kg  ) 
                         - fvCoeffs.Cont.AP[s](ig, jg, kg) * fields.P( ig  , jg-1, kg  ) 
                         - fvCoeffs.Cont.AP[w](ig, jg, kg) * fields.P( ig-1, jg  , kg  ) 
                         - fvCoeffs.Cont.AP[t](ig, jg, kg) * fields.P( ig  , jg  , kg+1) 
                         - fvCoeffs.Cont.AP[b](ig, jg, kg) * fields.P( ig  , jg  , kg-1)

                         - pressureWideStencil
                         
                         - fvCoeffs.Cont.B(ig, jg, kg)
                         + fvCoeffs.Cont.F(ig, jg, kg);

    return residual;
}



// Calculate the absolute residual of each equation from the finite volume stencil
template< MomentumInterpolation MI,
          Linearisation LI > [[ maybe_unused ]]
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
                residualU      += mask(ig, jg, kg) * abs( CellMomentumResidual_x<LI>(fields, fvCoeffs, ig, jg, kg) );
                scalingFactorU += mask(ig, jg, kg) * abs( fvCoeffs.Mom[X].AU[X][p](ig, jg, kg) * fields.U[X]( ig  , jg  , kg  ) );


                // V momentum
                residualV      += mask(ig, jg, kg) * abs( CellMomentumResidual_y<LI>(fields, fvCoeffs, ig, jg, kg) );
                scalingFactorV += mask(ig, jg, kg) * abs( fvCoeffs.Mom[Y].AU[Y][p](ig, jg, kg) * fields.U[Y]( ig  , jg  , kg  ) );


                // W momentm
                residualW      += mask(ig, jg, kg) * abs( CellMomentumResidual_z<LI>(fields, fvCoeffs, ig, jg, kg) );             
                scalingFactorW += mask(ig, jg, kg) * abs( fvCoeffs.Mom[Z].AU[Z][p](ig, jg, kg) * fields.U[Z]( ig  , jg  , kg  ) );


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
template< MomentumInterpolation MI,
          Linearisation LI >
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
                residuals.U[X](ig, jg, kg) = mask(ig, jg, kg) * CellMomentumResidual_x<LI>(fields, fvCoeffs, ig, jg, kg);

                // V momentum
                residuals.U[Y](ig, jg, kg) = mask(ig, jg, kg) * CellMomentumResidual_y<LI>(fields, fvCoeffs, ig, jg, kg);

                // W momentm
                residuals.U[Z](ig, jg, kg) = mask(ig, jg, kg) * CellMomentumResidual_z<LI>(fields, fvCoeffs, ig, jg, kg);

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
[[ maybe_unused ]]
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
[[ maybe_unused ]]
inline bool ResidualsDiverged( const FieldData<floatType> &residuals )
{
    for ( intType i = 0; i != FieldData<floatType>::nData; i++ ) {

        // Check if the compiler honours nans and infinities
        #ifdef CFD_HONOR_INFINITIES_AND_NANS
            if ( !std::isfinite( residuals[i] ) ) 
                return true;
        #else
            #define CFD_MAX_RESIDUAL_BEFORE_DIVERGENCE 1e8f
            if ( residuals[i] > CFD_MAX_RESIDUAL_BEFORE_DIVERGENCE )
                return true;
        #endif

    }
    return false;
}



// Update a vector of field probes
inline std::vector< FieldData<floatType> > SetFieldProbeValues( const FieldData<Tensor3D> &fields,
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


}   // end namespace CFD

#endif // RESIDUAL_FUNCTIONS