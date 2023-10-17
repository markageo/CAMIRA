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



// Calculate the absolute residual of each equation from the finite volume stencil
template< MomentumInterpolation MI,
          Linearisation LI > [[ maybe_unused ]]
inline FieldData<floatType> StencilResiduals( const FieldData<Tensor3D> &fields,
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
                floatType newtonStencilX = 0.0f;
                if constexpr ( LI == Linearisation::Newton ) {
                    newtonStencilX = fvCoeffs.Mom[X].AU[Y][n](i, j, k) * fields.U[Y]( ig  , jg+1, kg  )
                                   + fvCoeffs.Mom[X].AU[Y][p](i, j, k) * fields.U[Y]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[X].AU[Y][s](i, j, k) * fields.U[Y]( ig  , jg-1, kg  )

                                   + fvCoeffs.Mom[X].AU[Z][t](i, j, k) * fields.U[Z]( ig  , jg  , kg+1)
                                   + fvCoeffs.Mom[X].AU[Z][p](i, j, k) * fields.U[Z]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[X].AU[Z][b](i, j, k) * fields.U[Z]( ig  , jg  , kg-1);
                }
                residualU  += mask(i, j, k)
                            * abs( 
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

                                        + newtonStencilX

                                        - fvCoeffs.Mom[X].B(i, j, k)  );

                scalingFactorU += mask(i, j, k) * abs( fvCoeffs.Mom[X].AU[X][p](i, j, k) * fields.U[X]( ig  , jg  , kg  ) );


                // V momentum
                floatType newtonStencilY = 0.0f;
                if constexpr ( LI == Linearisation::Newton ) {
                    newtonStencilY = fvCoeffs.Mom[Y].AU[X][e](i, j, k) * fields.U[X]( ig+1, jg  , kg  )
                                   + fvCoeffs.Mom[Y].AU[X][p](i, j, k) * fields.U[X]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Y].AU[X][w](i, j, k) * fields.U[X]( ig-1, jg  , kg  )
                    
                                   + fvCoeffs.Mom[Y].AU[Z][t](i, j, k) * fields.U[Z]( ig  , jg  , kg+1)
                                   + fvCoeffs.Mom[Y].AU[Z][p](i, j, k) * fields.U[Z]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Y].AU[Z][b](i, j, k) * fields.U[Z]( ig  , jg  , kg-1);      
                }
                residualV   += mask(i, j, k) 
                             * abs( 
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

                                        + newtonStencilY

                                        - fvCoeffs.Mom[Y].B(i, j, k)  );

                scalingFactorV += mask(i, j, k) * abs( fvCoeffs.Mom[Y].AU[Y][p](i, j, k) * fields.U[Y]( ig  , jg  , kg  ) );


                // W momentm
                floatType newtonStencilZ = 0.0f;
                if constexpr ( LI == Linearisation::Newton ) {
                    newtonStencilZ = fvCoeffs.Mom[Z].AU[X][e](i, j, k) * fields.U[X]( ig+1, jg  , kg  )
                                   + fvCoeffs.Mom[Z].AU[X][p](i, j, k) * fields.U[X]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Z].AU[X][w](i, j, k) * fields.U[X]( ig-1, jg  , kg  )
                    
                                   + fvCoeffs.Mom[Z].AU[Y][n](i, j, k) * fields.U[Y]( ig  , jg+1, kg  )
                                   + fvCoeffs.Mom[Z].AU[Y][p](i, j, k) * fields.U[Y]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Z].AU[Y][s](i, j, k) * fields.U[Y]( ig  , jg-1, kg  );    
                }
                residualW   += mask(i, j, k)
                             * abs( 
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

                                        + newtonStencilZ

                                        - fvCoeffs.Mom[Z].B(i, j, k)  );

                scalingFactorW +=  mask(i, j, k) * abs( fvCoeffs.Mom[Z].AU[Z][p](i, j, k) * fields.U[Z]( ig  , jg  , kg  ) );


                // Continuity 
                floatType pressureWideStencil = 0.0f;
                if constexpr ( MI == MomentumInterpolation::Implicit ) {
                    pressureWideStencil = fvCoeffs.Cont.AP[nn](i, j, k) * fields.P( ig  , jg+2, kg  ) 
                                        + fvCoeffs.Cont.AP[ee](i, j, k) * fields.P( ig+2, jg  , kg  ) 
                                        + fvCoeffs.Cont.AP[ss](i, j, k) * fields.P( ig  , jg-2, kg  ) 
                                        + fvCoeffs.Cont.AP[ww](i, j, k) * fields.P( ig-2, jg  , kg  ) 
                                        + fvCoeffs.Cont.AP[tt](i, j, k) * fields.P( ig  , jg  , kg+2) 
                                        + fvCoeffs.Cont.AP[bb](i, j, k) * fields.P( ig  , jg  , kg-2);
                }
                residualP   += mask(i, j, k) 
                             * abs( 
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

                                  + pressureWideStencil
                                
                                  - fvCoeffs.Cont.B(i, j, k) );

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

    // ForAllFieldData( [&] (intType f) {
    //     residuals[f] /= static_cast<floatType>( fvCoeffs.nCells[X] * fvCoeffs.nCells[Y] * fvCoeffs.nCells[Z] ); 
    // } );

    return residuals;
}



// Calculate the local discrete equation residual field
template< MomentumInterpolation MI,
          Linearisation LI > [[ maybe_unused ]]
inline FieldData<Tensor3D> StencilResidualsField( const FieldData<Tensor3D> &fields,
                                                  const FVCoefficients &fvCoeffs,
                                                  const Tensor3D &mask )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    FieldData<Tensor3D> residuals( Tensor3D( fvCoeffs.nCells[X], fvCoeffs.nCells[Y], fvCoeffs.nCells[Z] ) );
    ForAllFieldData( [&] (intType f) {
        residuals[f].setZero();
    } );

    for ( intType k = 0; k != fvCoeffs.nCells[Z]; k++ ) {
        for ( intType j = 0; j != fvCoeffs.nCells[Y]; j++ ) {
            for ( intType i = 0; i != fvCoeffs.nCells[X]; i++ ) {

                intType ig{ G(i) }, jg{ G(j) }, kg{ G(k) };

                // U momentum
                floatType newtonStencilX = 0.0f;
                if constexpr ( LI == Linearisation::Newton ) {
                    newtonStencilX = fvCoeffs.Mom[X].AU[Y][n](i, j, k) * fields.U[Y]( ig  , jg+1, kg  )
                                   + fvCoeffs.Mom[X].AU[Y][p](i, j, k) * fields.U[Y]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[X].AU[Y][s](i, j, k) * fields.U[Y]( ig  , jg-1, kg  )

                                   + fvCoeffs.Mom[X].AU[Z][t](i, j, k) * fields.U[Z]( ig  , jg  , kg+1)
                                   + fvCoeffs.Mom[X].AU[Z][p](i, j, k) * fields.U[Z]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[X].AU[Z][b](i, j, k) * fields.U[Z]( ig  , jg  , kg-1);
                }
                residuals.U[X](i, j, k) += mask(i, j, k)
                                        * abs( 
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

                                                    + newtonStencilX

                                                    - fvCoeffs.Mom[X].B(i, j, k)  );



                // V momentum
                floatType newtonStencilY = 0.0f;
                if constexpr ( LI == Linearisation::Newton ) {
                    newtonStencilY = fvCoeffs.Mom[Y].AU[X][e](i, j, k) * fields.U[X]( ig+1, jg  , kg  )
                                   + fvCoeffs.Mom[Y].AU[X][p](i, j, k) * fields.U[X]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Y].AU[X][w](i, j, k) * fields.U[X]( ig-1, jg  , kg  )
                    
                                   + fvCoeffs.Mom[Y].AU[Z][t](i, j, k) * fields.U[Z]( ig  , jg  , kg+1)
                                   + fvCoeffs.Mom[Y].AU[Z][p](i, j, k) * fields.U[Z]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Y].AU[Z][b](i, j, k) * fields.U[Z]( ig  , jg  , kg-1);      
                }
                residuals.U[Y](i, j, k) += mask(i, j, k) 
                                        * abs( 
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

                                                    + newtonStencilY

                                                    - fvCoeffs.Mom[Y].B(i, j, k)  );


                // W momentm
                floatType newtonStencilZ = 0.0f;
                if constexpr ( LI == Linearisation::Newton ) {
                    newtonStencilZ = fvCoeffs.Mom[Z].AU[X][e](i, j, k) * fields.U[X]( ig+1, jg  , kg  )
                                   + fvCoeffs.Mom[Z].AU[X][p](i, j, k) * fields.U[X]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Z].AU[X][w](i, j, k) * fields.U[X]( ig-1, jg  , kg  )
                    
                                   + fvCoeffs.Mom[Z].AU[Y][n](i, j, k) * fields.U[Y]( ig  , jg+1, kg  )
                                   + fvCoeffs.Mom[Z].AU[Y][p](i, j, k) * fields.U[Y]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Z].AU[Y][s](i, j, k) * fields.U[Y]( ig  , jg-1, kg  );    
                }
                residuals.U[Z](i, j, k) += mask(i, j, k)
                                        * abs( 
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

                                                    + newtonStencilZ

                                                    - fvCoeffs.Mom[Z].B(i, j, k)  );


                // Continuity 
                floatType pressureWideStencil = 0.0f;
                if constexpr ( MI == MomentumInterpolation::Implicit ) {
                    pressureWideStencil = fvCoeffs.Cont.AP[nn](i, j, k) * fields.P( ig  , jg+2, kg  ) 
                                        + fvCoeffs.Cont.AP[ee](i, j, k) * fields.P( ig+2, jg  , kg  ) 
                                        + fvCoeffs.Cont.AP[ss](i, j, k) * fields.P( ig  , jg-2, kg  ) 
                                        + fvCoeffs.Cont.AP[ww](i, j, k) * fields.P( ig-2, jg  , kg  ) 
                                        + fvCoeffs.Cont.AP[tt](i, j, k) * fields.P( ig  , jg  , kg+2) 
                                        + fvCoeffs.Cont.AP[bb](i, j, k) * fields.P( ig  , jg  , kg-2);
                }
                residuals.P(i, j, k) += mask(i, j, k) 
                                     * abs( 
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

                                         + pressureWideStencil
                                        
                                         - fvCoeffs.Cont.B(i, j, k) );

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



// Normalise the residual by the first iteration
[[ maybe_unused ]]
inline void NormaliseResiduals( FieldData<floatType> &residuals,
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