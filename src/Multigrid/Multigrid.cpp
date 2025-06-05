#include "Multigrid.h"

#include "../FiniteVolume/Mesh.h"
#include "../Core/FVTools.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"

namespace CAMIRA
{

namespace
{



} // end anonymous namespace

template< MomentumInterpolation MI >
void SetMGLevels( std::vector< GridLevelData<MI > > &mgLevels, 
                  const InputData &inputData,
                  const AxisTransformationMap &axisTransofrmation )
{

    const InputData::MultigridSettings &mgSettings = inputData.multigridSettings;

    // Write the geometry to file
    if ( inputData.outputGeometry ) {
        WriteGeometryToFile(inputData, axisTransofrmation);
    }

    // Reserve data. The linear solver holds references to the fields which can break if the vector is resized
    mgLevels.reserve( mgSettings.maxCoarseLevels + 1 );

    for ( size_t level = 0; level != mgSettings.maxCoarseLevels + 1; level++ ) {

        if ( level == 0 ) {
            mgLevels.emplace_back();
            mgLevels[level].mesh = CreateMesh( inputData, axisTransofrmation );
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
        GridLevelData<MI> &mgl = mgLevels[level];
        mgl.level = level;


        // Boundary condition data
        mgl.bcData = SetBoundaryConditionData(inputData, mgl.mesh);


        // Immersed boundary data
        mgl.ibData = CreateImmersedBoundaryData(inputData, mgl.mesh);


        // Allocate and initialise fields
        mgl.fields = InitialiseFields(mgl.mesh, inputData, axisTransofrmation);
        if ( level == 0 ) {
            MaskFields(mgl.fields, mgl.ibData.mask);
        } else {
            RestrictFields( mgl.fields, mgLevels[level-1].fields, mgLevels[level-1].mesh, mgLevels[level].mesh, mgLevels[level].ibData.mask );
        }
        SetGhostCells(mgl.fields, mgl.mesh, mgl.bcData);
        mgl.fieldsRestricted = mgl.fields;
        
        switch ( inputData.schemes.timeScheme) {
            case TimeSchemes::BackwardsEuler:
                mgl.fieldsPrevTime        = mgl.fields;
                break;

            case TimeSchemes::BackwardsThreeLevel:
                mgl.fieldsPrevTime        = mgl.fields;
                mgl.fieldsPrevPrevTime    = mgl.fields;
                break;

            case TimeSchemes::Steady:
                /* NULL */
                break;
        }
        mgl.fieldsOld = mgl.fields;

        // Finite volume coefficients
        mgl.fvCoeffs = InitialiseFVCoefficients(mgl.mesh, inputData);

        // Face fluxes and advected velocities
        mgl.faceFluxes = InitialiseFaceFluxes(mgl.mesh, mgl.fields.U, mgl.bcData);

        // Set the coefficients that depend on linearisation
        UpdateFVCoefficients( mgl.fvCoeffs, mgl.mesh, mgl.fields, mgl.fieldsPrevTime, mgl.fieldsPrevPrevTime, mgl.faceFluxes, mgl.ibData );


        // Allocate and initialise residualsRestricted
        mgl.residuals = FieldData<Tensor3D>( Tensor3D( mgl.mesh.nCells(0) + 2*nGhost, 
                                                       mgl.mesh.nCells(1) + 2*nGhost, 
                                                       mgl.mesh.nCells(2) + 2*nGhost).setZero() );
        mgl.residualsRestricted = mgl.residuals;

        mgl.fineGridCorrection = FieldData<Tensor3D>( Tensor3D( mgl.mesh.nCells(0) + 2*nGhost, 
                                                                mgl.mesh.nCells(1) + 2*nGhost, 
                                                                mgl.mesh.nCells(2) + 2*nGhost).setZero() );

        mgl.coarseGridRightHandSide = FieldData<Tensor3D>( Tensor3D( mgl.mesh.nCells(0) + 2*nGhost, 
                                                                     mgl.mesh.nCells(1) + 2*nGhost, 
                                                                     mgl.mesh.nCells(2) + 2*nGhost).setZero() );

        // Linear Solver
        switch ( inputData.smootherSettings.type ) {
            case Smoothers::nestedLineSymmetricSerial:
                mgl.linearSolver = std::make_unique< nestedLineSymmetricSolverSerial<MI> >(mgl.fields, mgl.fieldsOld, mgl.ibData.mask, mgl.fvCoeffs, mgl.mesh, mgl.bcData, inputData.smootherSettings);
                break;
            
            case Smoothers::domainSymmetricSerial:
                mgl.linearSolver = std::make_unique< domainSymmetricSolverSerial<MI> >(mgl.fields, mgl.fieldsOld, mgl.ibData.mask, mgl.fvCoeffs, mgl.mesh, mgl.bcData, inputData.smootherSettings);
                break;

            case Smoothers::domainSymmetricParallel:
                mgl.linearSolver = std::make_unique< domainSymmetricSolverParallel<MI> >(mgl.fields, mgl.fieldsOld, mgl.ibData.mask, mgl.fvCoeffs, mgl.mesh, mgl.bcData, inputData.smootherSettings);
                break;
        }

        
        

    }
    mgLevels.back().isCoarsestLevel = true;
    mgLevels.shrink_to_fit();
}
template void SetMGLevels( std::vector< GridLevelData<MomentumInterpolation::Implicit     > > &, const InputData &, const AxisTransformationMap & );
template void SetMGLevels( std::vector< GridLevelData<MomentumInterpolation::SemiExplicit > > &, const InputData &, const AxisTransformationMap & );




void RestrictField( Tensor3D &coarseField,
                    const Tensor3D &fineField, 
                    const Mesh &fineMesh,
                    const Mesh &coarseMesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    enum class AgglomerationState { 
        AllAgglomerated, 
        LoNotAgglomerated, 
        HiNotAgglomerated 
    };

    EnumVector<Axis, AgglomerationState> agglomerationStates;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        if        ( fineMesh.nCells[axis] == 2 * coarseMesh.nCells[axis] ) {
            agglomerationStates[axis] = AgglomerationState::AllAgglomerated;
        } else if ( fineMesh.cellLengths[axis](0) == coarseMesh.cellLengths[axis](0) ) {
            agglomerationStates[axis] = AgglomerationState::LoNotAgglomerated;
        } else {
            agglomerationStates[axis] = AgglomerationState::HiNotAgglomerated;
        }
        
    } ); 


    auto SetFineIdicesAndInterpFactor = [&] ( floatType &lambda, intType &iF, intType &iFp1, const intType iC, const Axis::ENUMDATA axis ) {
        switch (agglomerationStates[axis]) {
            case AgglomerationState::AllAgglomerated:
                iF     = 2 * iC;
                iFp1   = iF + 1;
                lambda = ( coarseMesh.cellCenters[axis](iC) - fineMesh.cellCenters[axis](iF) ) / ( fineMesh.cellCenters[axis](iFp1) - fineMesh.cellCenters[axis](iF) );
                break;

            case AgglomerationState::LoNotAgglomerated:
                if ( iC == 0 ) {
                    iF     = 0;
                    iFp1   = 0;
                    lambda = 0.0f;
                } else {
                    iF     = 2 * iC - 1;
                    iFp1   = iF + 1;
                    lambda = ( coarseMesh.cellCenters[axis](iC) - fineMesh.cellCenters[axis](iF) ) / ( fineMesh.cellCenters[axis](iFp1) - fineMesh.cellCenters[axis](iF) );
                }
                break;

            case AgglomerationState::HiNotAgglomerated:
                iF = 2 * iC;
                if ( iC < coarseMesh.nCells[axis]-1 ) {
                    iFp1   = iF + 1;
                    lambda = ( coarseMesh.cellCenters[axis](iC) - fineMesh.cellCenters[axis](iF) ) / ( fineMesh.cellCenters[axis](iFp1) - fineMesh.cellCenters[axis](iF) );
                } else {
                    iFp1   = iF;
                    lambda = 0.0f;
                }
                break;
        }
    };


    // Iterate coarse grid
    #pragma omp parallel for
    for ( intType kC = 0; kC != coarseMesh.nCells(Z); kC++ ) {

        intType kF(0), kFp1(0);
        floatType lambdaZ(0.0f);
        SetFineIdicesAndInterpFactor(lambdaZ, kF, kFp1, kC, Z);

        for ( intType jC = 0; jC != coarseMesh.nCells(Y); jC++ ) {

            intType jF(0), jFp1(0);
            floatType lambdaY(0.0f);
            SetFineIdicesAndInterpFactor(lambdaY, jF, jFp1, jC, Y);

            for ( intType iC = 0; iC != coarseMesh.nCells(X); iC++ ) {

                intType iF(0), iFp1(0);
                floatType lambdaX(0.0f);
                SetFineIdicesAndInterpFactor(lambdaX, iF, iFp1, iC, X);

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

}



void ProlongateField( Tensor3D &fineField, 
                      const Tensor3D &coarseField,
                      const Mesh &coarseMesh,
                      const Mesh &fineMesh )
{

    using enum Axis::ENUMDATA;
    using FVT::G;

    EnumVector<Axis, bool> firstCellNotAgglomerated;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        firstCellNotAgglomerated[axis] = coarseMesh.cellLengths[axis](0) == fineMesh.cellLengths[axis](0);
    } );

    auto SetCoarseIdicesAndInterpFactor = [&] ( floatType &lambda, intType &iC, intType &iCp1, const intType iF, const Axis::ENUMDATA axis ) {
        if ( iF == 0 ) {    // Edge nodes get direction injection
            iC = 0;
            iCp1 = iC;
            lambda = 0.0f;
        } else if ( iF == fineMesh.nCells(axis) - 1 ) {
            iC = coarseMesh.nCells(axis) - 1;
            iCp1 = iC;
            lambda = 0.0f;
        } else {
            iC = firstCellNotAgglomerated[axis] ? static_cast<intType>( floor( static_cast<floatType>(iF) / 2.0f ) )
                                                : static_cast<intType>( floor( static_cast<floatType>(iF - 1) / 2.0f ) );
            iCp1 = iC + 1;
            lambda = ( fineMesh.cellCenters[axis](iF) - coarseMesh.cellCenters[axis](iC) ) 
                   / ( coarseMesh.cellCenters[axis](iCp1) - coarseMesh.cellCenters[axis](iC) );
        }
    };

    // Iterate fine grid
    #pragma omp parallel for
    for ( intType kF = 0; kF != fineMesh.nCells(Z); kF++ ) {
        
        intType kC, kCp1;
        floatType lambdaZ;
        SetCoarseIdicesAndInterpFactor( lambdaZ, kC, kCp1, kF, Z );

        for ( intType jF = 0; jF != fineMesh.nCells(Y); jF++ ) {
            
            intType jC, jCp1;
            floatType lambdaY;
            SetCoarseIdicesAndInterpFactor( lambdaY, jC, jCp1, jF, Y );

            for ( intType iF = 0; iF != fineMesh.nCells(X); iF++ ) {

                intType iC, iCp1;
                floatType lambdaX;
                SetCoarseIdicesAndInterpFactor( lambdaX, iC, iCp1, iF, X );

                // Points to interpolation from
                floatType c000 = coarseField( G(iC  , jC  , kC  ) ),
                          c100 = coarseField( G(iCp1, jC  , kC  ) ),
                          c010 = coarseField( G(iC  , jCp1, kC  ) ),
                          c001 = coarseField( G(iC  , jC  , kCp1) ),
                          c101 = coarseField( G(iCp1, jC  , kCp1) ),
                          c011 = coarseField( G(iC  , jCp1, kCp1) ),
                          c110 = coarseField( G(iCp1, jCp1, kC  ) ),
                          c111 = coarseField( G(iCp1, jCp1, kCp1) );

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

                fineField( G(iF, jF, kF) ) = c;

            }
        }
    }
}



void RestrictFields( FieldData<Tensor3D> &fieldsRestricted,
                     const FieldData<Tensor3D> &fields,
                     const Mesh &fineMesh,
                     const Mesh &coarseMesh,
                     const Tensor3D &mask )
{
    ForAllFieldData( [&] (intType f) {
        RestrictField( fieldsRestricted[f],
                       fields[f],
                       fineMesh, 
                       coarseMesh );
    } );
    MaskFields( fieldsRestricted, mask );
}



void ComputeFineGridCorrection( FieldData<Tensor3D> &fineGridCorrection,
                                const FieldData<Tensor3D> &coarseGridSolution,
                                const FieldData<Tensor3D> &restrictedFineGridApproximation,
                                const Mesh &coarseMesh,
                                const Mesh &fineMesh,
                                const Tensor3D &mask )
{
    FieldData<Tensor3D> coarseGridError;
    ForAllFieldData( [&] (intType f) {
        coarseGridError[f]    = coarseGridSolution[f] - restrictedFineGridApproximation[f];
        ProlongateField( fineGridCorrection[f], coarseGridError[f], coarseMesh, fineMesh );
    } );
    MaskFields( fineGridCorrection, mask );
}



template< MomentumInterpolation MI >
void CalculateCoarseGridRightHandSide( FieldData<Tensor3D> &coarseGridRightHandSide,
                                       const FVCoefficients &fvCoeffs, 
                                       const FieldData<Tensor3D> &fieldsRestricted,
                                       const FieldData<Tensor3D> &residualsRestricted,
                                       const Tensor3D &mask )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    #pragma omp parallel for collapse(3)
    for ( intType k = 0; k != fvCoeffs.nCells[Z]; k++ ) {
        for ( intType j = 0; j != fvCoeffs.nCells[Y]; j++ ) {
            for ( intType i = 0; i != fvCoeffs.nCells[X]; i++ ) {

                intType ig{ G(i) }, jg{ G(j) }, kg{ G(k) };

                // U momentum
                auto & xMomCoeffs = fvCoeffs.Mom[X];
                coarseGridRightHandSide.U[X](ig, jg, kg) = mask(ig, jg, kg)
                                                         * (  residualsRestricted.U[X](ig, jg, kg)

                                                            + xMomCoeffs.AU[p](ig, jg, kg) * fieldsRestricted.U[X]( ig  , jg  , kg  ) 
                                                            + xMomCoeffs.AU[n](ig, jg, kg) * fieldsRestricted.U[X]( ig  , jg+1, kg  ) 
                                                            + xMomCoeffs.AU[e](ig, jg, kg) * fieldsRestricted.U[X]( ig+1, jg  , kg  ) 
                                                            + xMomCoeffs.AU[s](ig, jg, kg) * fieldsRestricted.U[X]( ig  , jg-1, kg  ) 
                                                            + xMomCoeffs.AU[w](ig, jg, kg) * fieldsRestricted.U[X]( ig-1, jg  , kg  ) 
                                                            + xMomCoeffs.AU[t](ig, jg, kg) * fieldsRestricted.U[X]( ig  , jg  , kg+1) 
                                                            + xMomCoeffs.AU[b](ig, jg, kg) * fieldsRestricted.U[X]( ig  , jg  , kg-1) 

                                                            + xMomCoeffs.AP[e](ig) * fieldsRestricted.P( ig+1, jg  , kg  )
                                                            + xMomCoeffs.AP[p](ig) * fieldsRestricted.P( ig  , jg  , kg  )
                                                            + xMomCoeffs.AP[w](ig) * fieldsRestricted.P( ig-1, jg  , kg  )

                                                            + xMomCoeffs.B(ig, jg, kg) );



                // V momentum
                auto & yMomCoeffs = fvCoeffs.Mom[Y];
                coarseGridRightHandSide.U[Y](ig, jg, kg) = mask(ig, jg, kg) 
                                                         * (  residualsRestricted.U[Y](ig, jg, kg)

                                                            + yMomCoeffs.AU[p](ig, jg, kg) * fieldsRestricted.U[Y]( ig  , jg  , kg  ) 
                                                            + yMomCoeffs.AU[n](ig, jg, kg) * fieldsRestricted.U[Y]( ig  , jg+1, kg  ) 
                                                            + yMomCoeffs.AU[e](ig, jg, kg) * fieldsRestricted.U[Y]( ig+1, jg  , kg  ) 
                                                            + yMomCoeffs.AU[s](ig, jg, kg) * fieldsRestricted.U[Y]( ig  , jg-1, kg  ) 
                                                            + yMomCoeffs.AU[w](ig, jg, kg) * fieldsRestricted.U[Y]( ig-1, jg  , kg  ) 
                                                            + yMomCoeffs.AU[t](ig, jg, kg) * fieldsRestricted.U[Y]( ig  , jg  , kg+1) 
                                                            + yMomCoeffs.AU[b](ig, jg, kg) * fieldsRestricted.U[Y]( ig  , jg  , kg-1) 

                                                            + yMomCoeffs.AP[n](jg) * fieldsRestricted.P( ig  , jg+1, kg  )
                                                            + yMomCoeffs.AP[p](jg) * fieldsRestricted.P( ig  , jg  , kg  )
                                                            + yMomCoeffs.AP[s](jg) * fieldsRestricted.P( ig  , jg-1, kg  )

                                                            + yMomCoeffs.B(ig, jg, kg) );


                // W momentum
                auto & zMomCoeffs = fvCoeffs.Mom[Z];
                coarseGridRightHandSide.U[Z](ig, jg, kg) = mask(ig, jg, kg)
                                                         * (  residualsRestricted.U[Z](ig, jg, kg)

                                                            + zMomCoeffs.AU[p](ig, jg, kg) * fieldsRestricted.U[Z]( ig  , jg  , kg  ) 
                                                            + zMomCoeffs.AU[n](ig, jg, kg) * fieldsRestricted.U[Z]( ig  , jg+1, kg  ) 
                                                            + zMomCoeffs.AU[e](ig, jg, kg) * fieldsRestricted.U[Z]( ig+1, jg  , kg  ) 
                                                            + zMomCoeffs.AU[s](ig, jg, kg) * fieldsRestricted.U[Z]( ig  , jg-1, kg  ) 
                                                            + zMomCoeffs.AU[w](ig, jg, kg) * fieldsRestricted.U[Z]( ig-1, jg  , kg  ) 
                                                            + zMomCoeffs.AU[t](ig, jg, kg) * fieldsRestricted.U[Z]( ig  , jg  , kg+1) 
                                                            + zMomCoeffs.AU[b](ig, jg, kg) * fieldsRestricted.U[Z]( ig  , jg  , kg-1) 

                                                            + zMomCoeffs.AP[t](kg) * fieldsRestricted.P( ig  , jg  , kg+1)
                                                            + zMomCoeffs.AP[p](kg) * fieldsRestricted.P( ig  , jg  , kg  )
                                                            + zMomCoeffs.AP[b](kg) * fieldsRestricted.P( ig  , jg  , kg-1)

                                                            + zMomCoeffs.B(ig, jg, kg) );



                // Continuity 
                auto &contCoeffs = fvCoeffs.Cont;
                floatType pressureWideStencil = 0.0f;
                if constexpr ( MI == MomentumInterpolation::Implicit ) {
                    pressureWideStencil = contCoeffs.AP[nn](ig, jg, kg) * fieldsRestricted.P( ig  , jg+2, kg  ) 
                                        + contCoeffs.AP[ee](ig, jg, kg) * fieldsRestricted.P( ig+2, jg  , kg  ) 
                                        + contCoeffs.AP[ss](ig, jg, kg) * fieldsRestricted.P( ig  , jg-2, kg  ) 
                                        + contCoeffs.AP[ww](ig, jg, kg) * fieldsRestricted.P( ig-2, jg  , kg  ) 
                                        + contCoeffs.AP[tt](ig, jg, kg) * fieldsRestricted.P( ig  , jg  , kg+2) 
                                        + contCoeffs.AP[bb](ig, jg, kg) * fieldsRestricted.P( ig  , jg  , kg-2);
                }
               coarseGridRightHandSide.P(ig, jg, kg) = mask(ig, jg, kg) 
                                                     * (  residualsRestricted.P(ig, jg, kg)

                                                        + contCoeffs.AU[X][e](ig) * fieldsRestricted.U[X]( ig+1, jg  , kg  )
                                                        + contCoeffs.AU[X][p](ig) * fieldsRestricted.U[X]( ig  , jg  , kg  )
                                                        + contCoeffs.AU[X][w](ig) * fieldsRestricted.U[X]( ig-1, jg  , kg  )

                                                        + contCoeffs.AU[Y][n](jg) * fieldsRestricted.U[Y]( ig  , jg+1, kg  )
                                                        + contCoeffs.AU[Y][p](jg) * fieldsRestricted.U[Y]( ig  , jg  , kg  )
                                                        + contCoeffs.AU[Y][s](jg) * fieldsRestricted.U[Y]( ig  , jg-1, kg  )

                                                        + contCoeffs.AU[Z][t](kg) * fieldsRestricted.U[Z]( ig  , jg  , kg+1)
                                                        + contCoeffs.AU[Z][p](kg) * fieldsRestricted.U[Z]( ig  , jg  , kg  )
                                                        + contCoeffs.AU[Z][b](kg) * fieldsRestricted.U[Z]( ig  , jg  , kg-1)

                                                        + contCoeffs.AP[p](ig, jg, kg) * fieldsRestricted.P( ig  , jg  , kg  )
                                                        + contCoeffs.AP[n](ig, jg, kg) * fieldsRestricted.P( ig  , jg+1, kg  ) 
                                                        + contCoeffs.AP[e](ig, jg, kg) * fieldsRestricted.P( ig+1, jg  , kg  ) 
                                                        + contCoeffs.AP[s](ig, jg, kg) * fieldsRestricted.P( ig  , jg-1, kg  ) 
                                                        + contCoeffs.AP[w](ig, jg, kg) * fieldsRestricted.P( ig-1, jg  , kg  ) 
                                                        + contCoeffs.AP[t](ig, jg, kg) * fieldsRestricted.P( ig  , jg  , kg+1) 
                                                        + contCoeffs.AP[b](ig, jg, kg) * fieldsRestricted.P( ig  , jg  , kg-1)

                                                        + pressureWideStencil
                                        
                                                        + contCoeffs.B(ig, jg, kg) );

            }
        }
    }
}
template void CalculateCoarseGridRightHandSide<MomentumInterpolation::Implicit     >( FieldData<Tensor3D> &, const FVCoefficients &, const FieldData<Tensor3D> &, const FieldData<Tensor3D> &, const Tensor3D & );
template void CalculateCoarseGridRightHandSide<MomentumInterpolation::SemiExplicit >( FieldData<Tensor3D> &, const FVCoefficients &, const FieldData<Tensor3D> &, const FieldData<Tensor3D> &, const Tensor3D & );




}   // end namespace CAMIRA