#include "Multigrid.h"

#include "../FiniteVolume/Mesh.h"
#include "../Core/FVTools.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"

namespace CFD
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
        if ( level == 0 ) {
            mgl.fields = InitialiseFields(mgl.mesh, inputData, axisTransofrmation);
            MaskFields(mgl.fields, mgl.ibData.mask);
        } else {
            mgl.fields = RestrictFields( mgLevels[level-1].fields, mgLevels[level-1].mesh, mgLevels[level].mesh, mgLevels[level].ibData.mask );
        }
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
        mgl.residualsRestricted = FieldData<Tensor3D>( Tensor3D( mgl.mesh.nCells(0) + 2*nGhost, 
                                                                 mgl.mesh.nCells(1) + 2*nGhost, 
                                                                 mgl.mesh.nCells(2) + 2*nGhost).setZero() );


        // Linear Solver
        switch ( inputData.smootherSettings.type ) {
            case Smoothers::nestedLineSymmetric:
                mgl.linearSolver = std::make_unique< nestedLineSymmetricSolverSerial<MI> >(mgl.fields, mgl.fieldsOld, mgl.ibData.mask, mgl.fvCoeffs, mgl.mesh, mgl.bcData, inputData.smootherSettings);
                break;
            
            case Smoothers::domainSymmetric:
                mgl.linearSolver = std::make_unique< domainSymmetricSolverSerial<MI> >(mgl.fields, mgl.fieldsOld, mgl.ibData.mask, mgl.fvCoeffs, mgl.mesh, mgl.bcData, inputData.smootherSettings);
                break;
        }

        // mgl.linearSolver = std::make_unique< domainSymmetricSolverParallel<MI> >(mgl.fields, mgl.fieldsOld, mgl.ibData.mask, mgl.fvCoeffs, mgl.mesh, mgl.bcData, inputData.smootherSettings);
        

    }
    mgLevels.back().isCoarsestLevel = true;
    mgLevels.shrink_to_fit();
}
template void SetMGLevels( std::vector< GridLevelData<MomentumInterpolation::Implicit     > > &, const InputData &, const AxisTransformationMap & );
template void SetMGLevels( std::vector< GridLevelData<MomentumInterpolation::SemiExplicit > > &, const InputData &, const AxisTransformationMap & );




Tensor3D RestrictField( const Tensor3D &fineField, 
                        const Mesh &fineMesh,
                        const Mesh &coarseMesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    Tensor3D coarseField( coarseMesh.nCells(0) + 2*CFD::nGhost, 
                          coarseMesh.nCells(1) + 2*CFD::nGhost, 
                          coarseMesh.nCells(2) + 2*CFD::nGhost );
    coarseField.setZero();

    auto SetFineIdicesAndInterpFactor = [&] ( floatType &lambda, intType &iFp1, intType &iFIncrement, const intType &iF, const intType iC, const Axis::ENUMDATA axis ) {
        if ( coarseMesh.cellLengths[axis](iC) == fineMesh.cellLengths[axis](iF) ) { 
            iFp1 = iF;
            iFIncrement = 1;
            lambda = 0.0f;
        } else {
            iFp1 = iF + 1;
            iFIncrement = 2;
            lambda = ( coarseMesh.cellCenters[axis](iC) - fineMesh.cellCenters[axis](iF) ) / ( fineMesh.cellCenters[axis](iFp1) - fineMesh.cellCenters[axis](iF) );
        }
    };


    // Iterate coarse grid
    for ( intType kC = 0, kF = 0, kFIncrement = 2; kC != coarseMesh.nCells(Z); kC++, kF += kFIncrement ) {

        intType kFp1;
        floatType lambdaZ;
        SetFineIdicesAndInterpFactor(lambdaZ, kFp1, kFIncrement, kF, kC, Z);

        for ( intType jC = 0, jF = 0, jFIncrement = 2; jC != coarseMesh.nCells(Y); jC++, jF += jFIncrement ) {

            intType jFp1;
            floatType lambdaY;
            SetFineIdicesAndInterpFactor(lambdaY, jFp1, jFIncrement, jF, jC, Y);

            for ( intType iC = 0, iF = 0, iFIncrement = 2; iC != coarseMesh.nCells(X); iC++, iF += iFIncrement ) {

                intType iFp1;
                floatType lambdaX;
                SetFineIdicesAndInterpFactor(lambdaX, iFp1, iFIncrement, iF, iC, X);

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
    fineField.setZero();

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

    return fineField;
}



// Tensor3D ProlongateField( const Tensor3D &coarseField,
//                           const Mesh &coarseMesh,
//                           const Mesh &fineMesh )
// {

//     using enum Axis::ENUMDATA;
//     using FVT::G;

//     Tensor3D fineField( fineMesh.nCells(0) + 2*CFD::nGhost, 
//                         fineMesh.nCells(1) + 2*CFD::nGhost, 
//                         fineMesh.nCells(2) + 2*CFD::nGhost );
//     fineField.setZero();

//     EnumVector<Axis, bool> firstCellNotAgglomerated;
//     EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
//         firstCellNotAgglomerated[axis] = coarseMesh.cellLengths[axis](0) == fineMesh.cellLengths[axis](0);
//     } );

//     auto GetCoarseIndex = [&] ( const intType iF, const Axis::ENUMDATA axis ) -> intType {
//         return firstCellNotAgglomerated[axis] ? static_cast<intType>( ceil( static_cast<floatType>(iF) / 2.0f ) )
//                                               : static_cast<intType>( floor( iF / 2 ) );
//     };

//     // Iterate fine grid
//     for ( intType kF = 0; kF != fineMesh.nCells(Z); kF++ ) {
        
//         intType kC = GetCoarseIndex(kF, Z);

//         for ( intType jF = 0; jF != fineMesh.nCells(Y); jF++ ) {
            
//             intType jC = GetCoarseIndex(jF, Y);

//             for ( intType iF = 0; iF != fineMesh.nCells(X); iF++ ) {

//                 intType iC = GetCoarseIndex(iF, X);

//                 // Injection
//                 fineField( G(iF, jF, kF) ) = coarseField( G(iC, jC, kC) );

//             }
//         }
//     }

//     return fineField;
// }


FieldData<Tensor3D> RestrictFields( const FieldData<Tensor3D> fields,
                                    const Mesh &fineMesh,
                                    const Mesh &coarseMesh,
                                    const Tensor3D &mask )
{
    FieldData<Tensor3D> fieldsRestricted;
    ForAllFieldData( [&] (intType f) {
        fieldsRestricted[f] = RestrictField( fields[f],
                                             fineMesh, 
                                             coarseMesh );
    } );
    MaskFields( fieldsRestricted, mask );
    return fieldsRestricted;
}



FieldData<Tensor3D> ComputeFineGridCorrection( const FieldData<Tensor3D> &coarseGridSolution,
                                               const FieldData<Tensor3D> &restrictedFineGridApproximation,
                                               const Mesh &coarseMesh,
                                               const Mesh &fineMesh,
                                               const Tensor3D &mask )
{
    FieldData<Tensor3D> coarseGridError, fineGridCorrection;
    ForAllFieldData( [&] (intType f) {
        coarseGridError[f]    = coarseGridSolution[f] - restrictedFineGridApproximation[f];
        fineGridCorrection[f] = ProlongateField( coarseGridError[f], coarseMesh, fineMesh );
    } );
    MaskFields( fineGridCorrection, mask );
    
    return fineGridCorrection;
}


template< MomentumInterpolation MI >
void TransformToCoarseGridEquations( FVCoefficients &fvCoeffs, 
                                     const FieldData<Tensor3D> &fieldsRestricted,
                                     const FieldData<Tensor3D> &residualsRestricted,
                                     const Tensor3D &mask )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    for ( intType k = 0; k != fvCoeffs.nCells[Z]; k++ ) {
        for ( intType j = 0; j != fvCoeffs.nCells[Y]; j++ ) {
            for ( intType i = 0; i != fvCoeffs.nCells[X]; i++ ) {

                intType ig{ G(i) }, jg{ G(j) }, kg{ G(k) };

                // U momentum
                auto & xMomCoeffs = fvCoeffs.Mom[X];
                xMomCoeffs.F(ig, jg, kg) += mask(ig, jg, kg)
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
                yMomCoeffs.F(ig, jg, kg) += mask(ig, jg, kg) 
                                          * (   residualsRestricted.U[Y](ig, jg, kg)

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



                // W momentm
                auto & zMomCoeffs = fvCoeffs.Mom[Z];
                zMomCoeffs.F(ig, jg, kg) += mask(ig, jg, kg)
                                          * (   residualsRestricted.U[Z](ig, jg, kg)

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
                auto & contCoeffs = fvCoeffs.Cont;
                floatType pressureWideStencil = 0.0f;
                if constexpr ( MI == MomentumInterpolation::Implicit ) {
                    pressureWideStencil = contCoeffs.AP[nn](ig, jg, kg) * fieldsRestricted.P( ig  , jg+2, kg  ) 
                                        + contCoeffs.AP[ee](ig, jg, kg) * fieldsRestricted.P( ig+2, jg  , kg  ) 
                                        + contCoeffs.AP[ss](ig, jg, kg) * fieldsRestricted.P( ig  , jg-2, kg  ) 
                                        + contCoeffs.AP[ww](ig, jg, kg) * fieldsRestricted.P( ig-2, jg  , kg  ) 
                                        + contCoeffs.AP[tt](ig, jg, kg) * fieldsRestricted.P( ig  , jg  , kg+2) 
                                        + contCoeffs.AP[bb](ig, jg, kg) * fieldsRestricted.P( ig  , jg  , kg-2);
                }
                contCoeffs.F(ig, jg, kg) += mask(ig, jg, kg) 
                                            * (   residualsRestricted.P(ig, jg, kg)

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
template void TransformToCoarseGridEquations<MomentumInterpolation::Implicit     >( FVCoefficients &, const FieldData<Tensor3D> &, const FieldData<Tensor3D> &, const Tensor3D & );
template void TransformToCoarseGridEquations<MomentumInterpolation::SemiExplicit >( FVCoefficients &, const FieldData<Tensor3D> &, const FieldData<Tensor3D> &, const Tensor3D & );


template< MomentumInterpolation MI >
FieldData<Tensor3D> CalculateCoarseGridRightHandSide( FVCoefficients &fvCoeffs, 
                                                      const FieldData<Tensor3D> &fieldsRestricted,
                                                      const FieldData<Tensor3D> &residualsRestricted,
                                                      const Tensor3D &mask )
{
    using enum Axis::ENUMDATA;
    using enum TransportCoefficients::ENUMDATA;

    FieldData<Tensor3D> coarseGridRightHandSide( Tensor3D( fvCoeffs.nCells[X] + 2*nGhost,
                                                           fvCoeffs.nCells[Y] + 2*nGhost,
                                                           fvCoeffs.nCells[Z] + 2*nGhost ).setZero() );    // Does not really need to setZero

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
                coarseGridRightHandSide.U[Y](ig, jg, kg) += mask(ig, jg, kg) 
                                              * (   residualsRestricted.U[Y](ig, jg, kg)

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
                coarseGridRightHandSide.U[Z](ig, jg, kg) += mask(ig, jg, kg)
                                              * (   residualsRestricted.U[Z](ig, jg, kg)

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
               coarseGridRightHandSide.P(ig, jg, kg) += mask(ig, jg, kg) 
                                          * (   residualsRestricted.P(ig, jg, kg)

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

    return coarseGridRightHandSide;
}
template FieldData<Tensor3D> CalculateCoarseGridRightHandSide<MomentumInterpolation::Implicit     >( FVCoefficients &, const FieldData<Tensor3D> &, const FieldData<Tensor3D> &, const Tensor3D & );
template FieldData<Tensor3D> CalculateCoarseGridRightHandSide<MomentumInterpolation::SemiExplicit >( FVCoefficients &, const FieldData<Tensor3D> &, const FieldData<Tensor3D> &, const Tensor3D & );




}   // end namespace CFD