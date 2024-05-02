#include "Multigrid.h"

#include "../FiniteVolume/Mesh.h"
#include "../Tools/FVTools.h"

namespace CFD
{

namespace
{


Tensor3D MomentumEquationOperator( const MomentumEquation &momEqn,
                                   const FieldData<Tensor3D> &fields )
{

}



Tensor3D ContinuityEquationOperator( const ContinuityEquation &conEqn,
                                     const FieldData<Tensor3D> &fields )
{


}


} // end anonymous namespace


template< MomentumInterpolation MI, Linearisation LI >
std::vector< GridLevelData<MI, LI> > CreateMGLevels( const InputData &inputData )
{

    const InputData::MultigridSettings &mgSettings = inputData.multigridSettings;

    std::vector<GridLevelData<MI, LI>> mgLevels;

    mgLevels.emplace_back();

    for ( intType level = 0; level != mgSettings.maxCoarseLevels; level++ ) {

        if ( !MeshCanBeCoarsened( mgLevels[level-1].mesh ) ) {
            break;
        }
        mgLevels.emplace_back();
        auto &mgl = mgLevels[level];

        mgl.isFinestLevel = false;
        mgl.isFinestLevel = false;

        // Mesh
        if ( level == 0 ) {
            mgl.mesh = CreateMesh( inputData );
            mgl.isFinestLevel = true;
        } else {
            mgl.mesh = CoarsenMesh( mgLevels[level-1].mesh, inputData.schemes.faceInterpolationScheme );
        }
        

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
        mgl.residuals = FieldData<Tensor3D>( Tensor3D( mgl.mesh.nCells(0), 
                                                       mgl.mesh.nCells(1), 
                                                       mgl.mesh.nCells(2) ).setZero() );

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

    return mgLevels;
}
template std::vector< GridLevelData<MomentumInterpolation::Implicit    , Linearisation::Picard> > CreateMGLevels( const InputData & );
template std::vector< GridLevelData<MomentumInterpolation::SemiExplicit, Linearisation::Picard> > CreateMGLevels( const InputData & );
template std::vector< GridLevelData<MomentumInterpolation::Implicit    , Linearisation::Newton> > CreateMGLevels( const InputData & );
template std::vector< GridLevelData<MomentumInterpolation::SemiExplicit, Linearisation::Newton> > CreateMGLevels( const InputData & );




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
            lambdaZ = ( coarseMesh.cellCenters[Z](kC) - fineMesh.cellCenters[Z](kF) ) / ( fineMesh.cellCenters[Z](kFp1) - fineMesh.cellCenters[Z](kFp1) );
        }

        for ( intType jC = 0, jF = 0, jFIncrement = 1; jC != coarseMesh.nCells(Y); jC++, jF += jFIncrement ) {

            intType jFp1;
            floatType lambdaY;
            if ( coarseMesh.cellLengths[Y](jC) == fineMesh.cellLengths[Y](jF) ) {
                jFp1 = jF;
                jFIncrement = 1;
                lambdaY = 0.0f;
            } else {
                jFp1 = kF + 1;
                jFIncrement = 2;
                lambdaY = ( coarseMesh.cellCenters[Y](jC) - fineMesh.cellCenters[Y](jF) ) / ( fineMesh.cellCenters[Y](jFp1) - fineMesh.cellCenters[Y](jFp1) );
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
                    lambdaX = ( coarseMesh.cellCenters[X](iC) - fineMesh.cellCenters[X](iF) ) / ( fineMesh.cellCenters[X](iFp1) - fineMesh.cellCenters[X](iFp1) );
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
        firstCellNotAgglomerated[axis] = coarseMesh.cellLengths[Z](0) == fineMesh.cellLengths[Z](0);
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



void TransformToCoarseGridEquations( FVCoefficients &fvCoeffs, 
                                     const FieldData<Tensor3D> &fieldsRestricted,
                                     const FieldData<Tensor3D> &residuals )
{

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        fvCoeffs.Mom[axis].F += residuals.U[axis]
                              + MomentumEquationOperator( fvCoeffs.Mom[axis], fieldsRestricted );

    } );

    fvCoeffs.Cont.F += residuals.P
                     + ContinuityEquationOperator( fvCoeffs.Cont, fieldsRestricted );



}



}   // end namespace CFD