#include "Multigrid.h"

#include "../FiniteVolume/Mesh.h"
#include "../Tools/FVTools.h"

namespace CFD
{

namespace
{



} // end anonymous namespace

template< MomentumInterpolation MI, Linearisation LI >
void SetMGLevels( std::vector< GridLevelData<MI, LI> > &mgLevels, 
                  const InputData &inputData )
{

    const InputData::MultigridSettings &mgSettings = inputData.multigridSettings;

    // Reserve data. The linear solver holds references to the fields which can break if the vector is resized
    mgLevels.reserve( mgSettings.maxCoarseLevels + 1 );

    for ( intType level = 0; level != mgSettings.maxCoarseLevels + 1; level++ ) {

        if ( level == 0 ) {
            mgLevels.emplace_back();
            mgLevels[level].mesh = CreateMesh( inputData );
            mgLevels[level].isFinestLevel   = true;
            mgLevels[level].isCoarsestLevel = false;
        } else if ( MeshCanBeCoarsened( mgLevels[level-1].mesh ) ) {
            mgLevels.emplace_back();
            mgLevels[level].mesh = CoarsenMesh( mgLevels[level-1].mesh, inputData.schemes.faceInterpolationScheme );
            mgLevels[level].isCoarsestLevel = false;
            mgLevels[level].isFinestLevel   = false;
        } else {
            break;
        }
        GridLevelData<MI, LI> &mgl = mgLevels[level];
        mgl.level = level;


        // Boundary condition data
        mgl.bcData = SetBoundaryConditionData(inputData, mgl.mesh);

        // Allocate and initialise fields
        mgl.fieldsRestricted = InitialiseFields(mgl.mesh, inputData);
        mgl.fields           = InitialiseFields(mgl.mesh, inputData);
        mgl.fieldsOld        = mgl.fields;

        // Face fluxes and advected velocities
        mgl.faceFluxes = InitialiseFaceFluxes(mgl.mesh, mgl.fields.U, mgl.bcData);
        if ( inputData.schemes.linearisation == Linearisation::Newton ) 
            mgl.faceAdvectedVelocities = InitialiseAdvectedFaceVelocities( mgl.mesh, mgl.fields.U, mgl.faceFluxes, mgl.bcData );

        // Allocate and initialise residuals
        mgl.residuals = FieldData<Tensor3D>( Tensor3D( mgl.mesh.nCells(0) + 2*nGhost, 
                                                       mgl.mesh.nCells(1) + 2*nGhost, 
                                                       mgl.mesh.nCells(2) + 2*nGhost).setZero() );

        // Immersed boundary data
        mgl.ibData = CreateImmersedBoundaryData(inputData, mgl.mesh);

        // Mask out solid geometries
        MaskFields(mgl.fieldsRestricted, mgl.ibData.mask);
        MaskFields(mgl.fields          , mgl.ibData.mask);
        MaskFields(mgl.fieldsOld       , mgl.ibData.mask);

        // Finite volume coefficients
        mgl.fvCoeffs = InitialiseFVCoefficients(mgl.mesh, mgl.fields, mgl.faceAdvectedVelocities, mgl.faceFluxes, mgl.ibData, mgl.bcData, inputData);

        // Linear Solver
        mgl.linearSolver = std::make_unique< LinearSolver<MI, LI> >(mgl.fields, mgl.fieldsOld, mgl.ibData.mask, mgl.fvCoeffs, inputData.linearSolverSettings);

    }
    mgLevels.back().isCoarsestLevel = true;
    mgLevels.shrink_to_fit();
}
template void SetMGLevels( std::vector< GridLevelData<MomentumInterpolation::Implicit    , Linearisation::Picard> > &, const InputData & );
template void SetMGLevels( std::vector< GridLevelData<MomentumInterpolation::SemiExplicit, Linearisation::Picard> > &, const InputData & );
template void SetMGLevels( std::vector< GridLevelData<MomentumInterpolation::Implicit    , Linearisation::Newton> > &, const InputData & );
template void SetMGLevels( std::vector< GridLevelData<MomentumInterpolation::SemiExplicit, Linearisation::Newton> > &, const InputData & );





Tensor3D RestrictField( const Tensor3D &fineField, 
                        const Mesh &fineMesh,
                        const Mesh &coarseMesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    Tensor3D coarseField( coarseMesh.nCells(0) + 2*CFD::nGhost, 
                          coarseMesh.nCells(1) + 2*CFD::nGhost, 
                          coarseMesh.nCells(2) + 2*CFD::nGhost );


    // Iterate coarse grid
    for ( intType kC = 0, kF = 0, kFIncrement = 1; kC != coarseMesh.nCells(Z); kC++, kF += kFIncrement ) {

        intType kFp1;
        floatType lambdaZ;
        if ( coarseMesh.cellLengths[Z](kC) == fineMesh.cellLengths[Z](kF) ) { // cell not agglomerated
            kFp1 = kF;
            kFIncrement = 1;
            lambdaZ = 0.0f;
        } else {
            kFp1 = kF + 1;
            kFIncrement = 2;
            lambdaZ = ( coarseMesh.cellCenters[Z](kC) - fineMesh.cellCenters[Z](kF) ) / ( fineMesh.cellCenters[Z](kFp1) - fineMesh.cellCenters[Z](kF) );
        }

        for ( intType jC = 0, jF = 0, jFIncrement = 1; jC != coarseMesh.nCells(Y); jC++, jF += jFIncrement ) {

            intType jFp1;
            floatType lambdaY;
            if ( coarseMesh.cellLengths[Y](jC) == fineMesh.cellLengths[Y](jF) ) {
                jFp1 = jF;
                jFIncrement = 1;
                lambdaY = 0.0f;
            } else {
                jFp1 = jF + 1;
                jFIncrement = 2;
                lambdaY = ( coarseMesh.cellCenters[Y](jC) - fineMesh.cellCenters[Y](jF) ) / ( fineMesh.cellCenters[Y](jFp1) - fineMesh.cellCenters[Y](jF) );
            }

            for ( intType iC = 0, iF = 0, iFIncrement = 1; iC != coarseMesh.nCells(X); iC++, iF += iFIncrement ) {

                intType iFp1;
                floatType lambdaX;
                if ( coarseMesh.cellLengths[X](iC) == fineMesh.cellLengths[X](iF) ) { 
                    iFp1 = iF;
                    iFIncrement = 1;
                    lambdaX = 0.0f;
                } else {
                    iFp1 = iF + 1;
                    iFIncrement = 2;
                    lambdaX = ( coarseMesh.cellCenters[X](iC) - fineMesh.cellCenters[X](iF) ) / ( fineMesh.cellCenters[X](iFp1) - fineMesh.cellCenters[X](iF) );
                }

                // Points to interpolation from
                floatType c000 = fineField( G(iF  , jF  , kF  ) ),
                          c100 = fineField( G(iFp1, jF  , kF  ) ),
                          c010 = fineField( G(iF  , jFp1, kF  ) ),
                          c001 = fineField( G(iF  , jF  , kFp1) ),
                          c101 = fineField( G(iFp1, jF  , kFp1) ),
                          c011 = fineField( G(iF  , jFp1, kFp1) ),
                          c110 = fineField( G(iFp1, jFp1, kF  ) ),
                          c111 = fineField( G(iFp1, jFp1, kFp1) );

                // Linear interpolation in z direction
                floatType c00 = ( 1-lambdaZ ) * c000  +  lambdaZ * c001,
                          c10 = ( 1-lambdaZ ) * c100  +  lambdaZ * c101,
                          c01 = ( 1-lambdaZ ) * c010  +  lambdaZ * c011,
                          c11 = ( 1-lambdaZ ) * c110  +  lambdaZ * c111;

                
                // Linear interpolation in y direction
                floatType c0 = ( 1-lambdaY ) * c00  +  lambdaY * c01,
                          c1 = ( 1-lambdaY ) * c10  +  lambdaY * c11;


                // Linear interpolation in x direction
                floatType c = ( 1-lambdaX ) * c0  +  lambdaX * c1;

                coarseField( G(iC, jC, kC) ) = c;

            }
        }
    }

    return coarseField;
}



Tensor3D ProlongateField( const Tensor3D &coarseField,
                          const Mesh &coarseMesh,
                          const Mesh &fineMesh )
{

    using enum Axis::ENUMDATA;
    using FVT::G;

    Tensor3D fineField( fineMesh.nCells(0) + 2*CFD::nGhost, 
                        fineMesh.nCells(1) + 2*CFD::nGhost, 
                        fineMesh.nCells(2) + 2*CFD::nGhost );

    EnumVector<Axis, bool> firstCellNotAgglomerated;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        firstCellNotAgglomerated[axis] = coarseMesh.cellLengths[axis](0) == fineMesh.cellLengths[axis](0);
    } );

    // Iterate fine grid
    for ( intType kF = 0; kF != fineMesh.nCells(Z); kF++ ) {
        
        intType kC = firstCellNotAgglomerated[Z] ? static_cast<intType>( ceil( static_cast<floatType>(kF) / 2.0f ) )
                                                 : static_cast<intType>( floor( kF / 2 ) );

        for ( intType jF = 0; jF != fineMesh.nCells(Y); jF++ ) {
            
            intType jC = firstCellNotAgglomerated[Y] ? static_cast<intType>( ceil( static_cast<floatType>(jF) / 2.0f ) )
                                                     : static_cast<intType>( floor( jF / 2 ) );


            for ( intType iF = 0; iF != fineMesh.nCells(X); iF++ ) {

                intType iC = firstCellNotAgglomerated[X] ? static_cast<intType>( ceil( static_cast<floatType>(iF) / 2.0f ) )
                                                         : static_cast<intType>( floor( iF / 2 ) );

                // Injection
                fineField( G(iF, jF, kF) ) = coarseField( G(iC, jC, kC) );

            }
        }
    }

    return fineField;
}




FieldData<Tensor3D> ComputeFineGridCorrection( const FieldData<Tensor3D> &coarseGridSolution,
                                               const FieldData<Tensor3D> &restrictedFineGridApproximation,
                                               const Mesh &coarseMesh,
                                               const Mesh &fineMesh )
{

    FieldData<Tensor3D> coarseGridError;
    ForAllFieldData( [&] (intType f) {
        coarseGridError[f] = coarseGridSolution[f] - restrictedFineGridApproximation[f];
    } );

    FieldData<Tensor3D> fineGridCorrection;
    ForAllFieldData( [&] (intType f) {
        fineGridCorrection[f] = ProlongateField( coarseGridError[f], coarseMesh, fineMesh );
    } );

    return fineGridCorrection;

}


template< MomentumInterpolation MI, Linearisation LI >
void TransformToCoarseGridEquations( FVCoefficients &fvCoeffs, 
                                     const FieldData<Tensor3D> &fieldsRestricted,
                                     const FieldData<Tensor3D> &residuals,
                                     const Tensor3D &mask )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    for ( intType k = 0; k != fvCoeffs.nCells[Z]; k++ ) {
        for ( intType j = 0; j != fvCoeffs.nCells[Y]; j++ ) {
            for ( intType i = 0; i != fvCoeffs.nCells[X]; i++ ) {

                intType ig{ G(i) }, jg{ G(j) }, kg{ G(k) };

                // U momentum
                floatType newtonStencilX = 0.0f;
                if constexpr ( LI == Linearisation::Newton ) {
                    newtonStencilX = fvCoeffs.Mom[X].AU[Y][n](i, j, k) * fieldsRestricted.U[Y]( ig  , jg+1, kg  )
                                   + fvCoeffs.Mom[X].AU[Y][p](i, j, k) * fieldsRestricted.U[Y]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[X].AU[Y][s](i, j, k) * fieldsRestricted.U[Y]( ig  , jg-1, kg  )

                                   + fvCoeffs.Mom[X].AU[Z][t](i, j, k) * fieldsRestricted.U[Z]( ig  , jg  , kg+1)
                                   + fvCoeffs.Mom[X].AU[Z][p](i, j, k) * fieldsRestricted.U[Z]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[X].AU[Z][b](i, j, k) * fieldsRestricted.U[Z]( ig  , jg  , kg-1);
                }
                fvCoeffs.Mom[X].F(i, j, k) += mask(i, j, k)
                                              * (  residuals.U[X](ig, jg, kg)

                                                 + fvCoeffs.Mom[X].AU[X][p](i, j, k) * fieldsRestricted.U[X]( ig  , jg  , kg  ) 
                                                 + fvCoeffs.Mom[X].AU[X][n](i, j, k) * fieldsRestricted.U[X]( ig  , jg+1, kg  ) 
                                                 + fvCoeffs.Mom[X].AU[X][e](i, j, k) * fieldsRestricted.U[X]( ig+1, jg  , kg  ) 
                                                 + fvCoeffs.Mom[X].AU[X][s](i, j, k) * fieldsRestricted.U[X]( ig  , jg-1, kg  ) 
                                                 + fvCoeffs.Mom[X].AU[X][w](i, j, k) * fieldsRestricted.U[X]( ig-1, jg  , kg  ) 
                                                 + fvCoeffs.Mom[X].AU[X][t](i, j, k) * fieldsRestricted.U[X]( ig  , jg  , kg+1) 
                                                 + fvCoeffs.Mom[X].AU[X][b](i, j, k) * fieldsRestricted.U[X]( ig  , jg  , kg-1) 

                                                 + fvCoeffs.Mom[X].AP[e](i) * fieldsRestricted.P( ig+1, jg  , kg  )
                                                 + fvCoeffs.Mom[X].AP[p](i) * fieldsRestricted.P( ig  , jg  , kg  )
                                                 + fvCoeffs.Mom[X].AP[w](i) * fieldsRestricted.P( ig-1, jg  , kg  )

                                                 + newtonStencilX

                                                 + fvCoeffs.Mom[X].B(i, j, k) );



                // V momentum
                floatType newtonStencilY = 0.0f;
                if constexpr ( LI == Linearisation::Newton ) {
                    newtonStencilY = fvCoeffs.Mom[Y].AU[X][e](i, j, k) * fieldsRestricted.U[X]( ig+1, jg  , kg  )
                                   + fvCoeffs.Mom[Y].AU[X][p](i, j, k) * fieldsRestricted.U[X]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Y].AU[X][w](i, j, k) * fieldsRestricted.U[X]( ig-1, jg  , kg  )
                    
                                   + fvCoeffs.Mom[Y].AU[Z][t](i, j, k) * fieldsRestricted.U[Z]( ig  , jg  , kg+1)
                                   + fvCoeffs.Mom[Y].AU[Z][p](i, j, k) * fieldsRestricted.U[Z]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Y].AU[Z][b](i, j, k) * fieldsRestricted.U[Z]( ig  , jg  , kg-1);      
                }
                fvCoeffs.Mom[Y].F(i, j, k) += mask(i, j, k) 
                                              * (   residuals.U[Y](ig, jg, kg)

                                                  + fvCoeffs.Mom[Y].AU[Y][p](i, j, k) * fieldsRestricted.U[Y]( ig  , jg  , kg  ) 
                                                  + fvCoeffs.Mom[Y].AU[Y][n](i, j, k) * fieldsRestricted.U[Y]( ig  , jg+1, kg  ) 
                                                  + fvCoeffs.Mom[Y].AU[Y][e](i, j, k) * fieldsRestricted.U[Y]( ig+1, jg  , kg  ) 
                                                  + fvCoeffs.Mom[Y].AU[Y][s](i, j, k) * fieldsRestricted.U[Y]( ig  , jg-1, kg  ) 
                                                  + fvCoeffs.Mom[Y].AU[Y][w](i, j, k) * fieldsRestricted.U[Y]( ig-1, jg  , kg  ) 
                                                  + fvCoeffs.Mom[Y].AU[Y][t](i, j, k) * fieldsRestricted.U[Y]( ig  , jg  , kg+1) 
                                                  + fvCoeffs.Mom[Y].AU[Y][b](i, j, k) * fieldsRestricted.U[Y]( ig  , jg  , kg-1) 

                                                  + fvCoeffs.Mom[Y].AP[n](j) * fieldsRestricted.P( ig  , jg+1, kg  )
                                                  + fvCoeffs.Mom[Y].AP[p](j) * fieldsRestricted.P( ig  , jg  , kg  )
                                                  + fvCoeffs.Mom[Y].AP[s](j) * fieldsRestricted.P( ig  , jg-1, kg  )

                                                  + newtonStencilY

                                                  + fvCoeffs.Mom[Y].B(i, j, k) );



                // W momentm
                floatType newtonStencilZ = 0.0f;
                if constexpr ( LI == Linearisation::Newton ) {
                    newtonStencilZ = fvCoeffs.Mom[Z].AU[X][e](i, j, k) * fieldsRestricted.U[X]( ig+1, jg  , kg  )
                                   + fvCoeffs.Mom[Z].AU[X][p](i, j, k) * fieldsRestricted.U[X]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Z].AU[X][w](i, j, k) * fieldsRestricted.U[X]( ig-1, jg  , kg  )
                    
                                   + fvCoeffs.Mom[Z].AU[Y][n](i, j, k) * fieldsRestricted.U[Y]( ig  , jg+1, kg  )
                                   + fvCoeffs.Mom[Z].AU[Y][p](i, j, k) * fieldsRestricted.U[Y]( ig  , jg  , kg  )
                                   + fvCoeffs.Mom[Z].AU[Y][s](i, j, k) * fieldsRestricted.U[Y]( ig  , jg-1, kg  );    
                }
                fvCoeffs.Mom[Z].F(i, j, k) += mask(i, j, k)
                                              * (   residuals.U[Z](ig, jg, kg)

                                                  + fvCoeffs.Mom[Z].AU[Z][p](i, j, k) * fieldsRestricted.U[Z]( ig  , jg  , kg  ) 
                                                  + fvCoeffs.Mom[Z].AU[Z][n](i, j, k) * fieldsRestricted.U[Z]( ig  , jg+1, kg  ) 
                                                  + fvCoeffs.Mom[Z].AU[Z][e](i, j, k) * fieldsRestricted.U[Z]( ig+1, jg  , kg  ) 
                                                  + fvCoeffs.Mom[Z].AU[Z][s](i, j, k) * fieldsRestricted.U[Z]( ig  , jg-1, kg  ) 
                                                  + fvCoeffs.Mom[Z].AU[Z][w](i, j, k) * fieldsRestricted.U[Z]( ig-1, jg  , kg  ) 
                                                  + fvCoeffs.Mom[Z].AU[Z][t](i, j, k) * fieldsRestricted.U[Z]( ig  , jg  , kg+1) 
                                                  + fvCoeffs.Mom[Z].AU[Z][b](i, j, k) * fieldsRestricted.U[Z]( ig  , jg  , kg-1) 

                                                  + fvCoeffs.Mom[Z].AP[t](k) * fieldsRestricted.P( ig  , jg  , kg+1)
                                                  + fvCoeffs.Mom[Z].AP[p](k) * fieldsRestricted.P( ig  , jg  , kg  )
                                                  + fvCoeffs.Mom[Z].AP[b](k) * fieldsRestricted.P( ig  , jg  , kg-1)

                                                  + newtonStencilZ

                                                  + fvCoeffs.Mom[Z].B(i, j, k) );



                // Continuity 
                floatType pressureWideStencil = 0.0f;
                if constexpr ( MI == MomentumInterpolation::Implicit ) {
                    pressureWideStencil = fvCoeffs.Cont.AP[nn](i, j, k) * fieldsRestricted.P( ig  , jg+2, kg  ) 
                                        + fvCoeffs.Cont.AP[ee](i, j, k) * fieldsRestricted.P( ig+2, jg  , kg  ) 
                                        + fvCoeffs.Cont.AP[ss](i, j, k) * fieldsRestricted.P( ig  , jg-2, kg  ) 
                                        + fvCoeffs.Cont.AP[ww](i, j, k) * fieldsRestricted.P( ig-2, jg  , kg  ) 
                                        + fvCoeffs.Cont.AP[tt](i, j, k) * fieldsRestricted.P( ig  , jg  , kg+2) 
                                        + fvCoeffs.Cont.AP[bb](i, j, k) * fieldsRestricted.P( ig  , jg  , kg-2);
                }
                fvCoeffs.Cont.F(i, j, k) += mask(i, j, k) 
                                          * (   residuals.P(ig, jg, kg)

                                              + fvCoeffs.Cont.AU[X][e](i) * fieldsRestricted.U[X]( ig+1, jg  , kg  )
                                              + fvCoeffs.Cont.AU[X][p](i) * fieldsRestricted.U[X]( ig  , jg  , kg  )
                                              + fvCoeffs.Cont.AU[X][w](i) * fieldsRestricted.U[X]( ig-1, jg  , kg  )

                                              + fvCoeffs.Cont.AU[Y][n](j) * fieldsRestricted.U[Y]( ig  , jg+1, kg  )
                                              + fvCoeffs.Cont.AU[Y][p](j) * fieldsRestricted.U[Y]( ig  , jg  , kg  )
                                              + fvCoeffs.Cont.AU[Y][s](j) * fieldsRestricted.U[Y]( ig  , jg-1, kg  )

                                              + fvCoeffs.Cont.AU[Z][t](k) * fieldsRestricted.U[Z]( ig  , jg  , kg+1)
                                              + fvCoeffs.Cont.AU[Z][p](k) * fieldsRestricted.U[Z]( ig  , jg  , kg  )
                                              + fvCoeffs.Cont.AU[Z][b](k) * fieldsRestricted.U[Z]( ig  , jg  , kg-1)

                                              + fvCoeffs.Cont.AP[p](i, j, k) * fieldsRestricted.P( ig  , jg  , kg  )
                                              + fvCoeffs.Cont.AP[n](i, j, k) * fieldsRestricted.P( ig  , jg+1, kg  ) 
                                              + fvCoeffs.Cont.AP[e](i, j, k) * fieldsRestricted.P( ig+1, jg  , kg  ) 
                                              + fvCoeffs.Cont.AP[s](i, j, k) * fieldsRestricted.P( ig  , jg-1, kg  ) 
                                              + fvCoeffs.Cont.AP[w](i, j, k) * fieldsRestricted.P( ig-1, jg  , kg  ) 
                                              + fvCoeffs.Cont.AP[t](i, j, k) * fieldsRestricted.P( ig  , jg  , kg+1) 
                                              + fvCoeffs.Cont.AP[b](i, j, k) * fieldsRestricted.P( ig  , jg  , kg-1)

                                              + pressureWideStencil
                                        
                                              + fvCoeffs.Cont.B(i, j, k) );

            }
        }
    }
    
}
template void TransformToCoarseGridEquations<MomentumInterpolation::Implicit    , Linearisation::Picard>( FVCoefficients &, const FieldData<Tensor3D> &, const FieldData<Tensor3D> &, const Tensor3D & );
template void TransformToCoarseGridEquations<MomentumInterpolation::SemiExplicit, Linearisation::Picard>( FVCoefficients &, const FieldData<Tensor3D> &, const FieldData<Tensor3D> &, const Tensor3D & );
template void TransformToCoarseGridEquations<MomentumInterpolation::Implicit    , Linearisation::Newton>( FVCoefficients &, const FieldData<Tensor3D> &, const FieldData<Tensor3D> &, const Tensor3D & );
template void TransformToCoarseGridEquations<MomentumInterpolation::SemiExplicit, Linearisation::Newton>( FVCoefficients &, const FieldData<Tensor3D> &, const FieldData<Tensor3D> &, const Tensor3D & );


}   // end namespace CFD