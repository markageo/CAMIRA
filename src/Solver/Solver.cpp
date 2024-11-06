#include "Solver.h"
#include "ConvergenceLogging.h"

#include "../Types.h"
#include "../Macros.h"
#include "../IO/InputProcessing.h"
#include "../Tools/SweepTransformations.h"
#include "../Tools/FieldProbe.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../Multigrid/Multigrid.h"

#include "../IO/ArrayIO.h"
#include "../IO/VTKWriter.h"

#include "ResidualFunctions.h"

#include <iostream>

namespace CFD
{

namespace
{

template< MomentumInterpolation MI, Linearisation LI >
void SetFineGridEquations( GridLevelData<MI, LI> &gridLevelData )
{
    auto &gld = gridLevelData; 
    UpdateIBData( gld.ibData, gld.fields );
    UpdateFaceFluxes( gld.faceFluxes, gld.mesh, gld.fields.U, gld.bcData);
    // UpdateFaceFluxesWithMWI( gld.faceFluxes, gld.mesh, gld.fields, gld.fvCoeffs, gld.bcData);
    SetIBFaceFluxes( gld.faceFluxes, gld.ibData );
    if constexpr ( LI == Linearisation::Newton ) {
        UpdateFaceAdvectedVelocities( gld.faceAdvectedVelocities, gld.mesh, gld.fvCoeffs, gld.fields.U, gld.faceFluxes, gld.bcData);
        switch ( gld.fvCoeffs.Mom[Axis::X].advectionScheme ) {
            case AdvectionSchemes::Upwind:
                SetIBFaceAdvectedVelocities<AdvectionSchemes::Upwind >( gld.faceAdvectedVelocities, gld.faceFluxes, gld.fields, gld.fvCoeffs, gld.mesh, gld.ibData );
                break;

            case AdvectionSchemes::Central:
                SetIBFaceAdvectedVelocities<AdvectionSchemes::Central>( gld.faceAdvectedVelocities, gld.faceFluxes, gld.fields, gld.fvCoeffs, gld.mesh, gld.ibData );
                break;

            case AdvectionSchemes::SOU:
                SetIBFaceAdvectedVelocities<AdvectionSchemes::SOU    >( gld.faceAdvectedVelocities, gld.faceFluxes, gld.fields, gld.fvCoeffs, gld.mesh, gld.ibData );
                break;

            case AdvectionSchemes::QUICK:
                SetIBFaceAdvectedVelocities<AdvectionSchemes::QUICK  >( gld.faceAdvectedVelocities, gld.faceFluxes, gld.fields, gld.fvCoeffs, gld.mesh, gld.ibData );
                break;
        }
    }
    SetGhostCells( gridLevelData.fields, gridLevelData.mesh, gridLevelData.bcData );
    UpdateFVCoefficients( gld.fvCoeffs, gld.mesh, gld.fields, gld.fieldsOld, gld.fieldsOldOld, gld.faceAdvectedVelocities, gld.faceFluxes, gld.ibData, gld.bcData);
}



template< MomentumInterpolation MI, Linearisation LI >
void SetCoarseGridRightHandSideInStencil( GridLevelData<MI, LI> &gridLevelData )
{
    auto &gld = gridLevelData;

    // Set fvCoeffs based on the restricted fine grid approximation for the RHS
    UpdateIBData( gld.ibData, gld.fieldsRestricted );
    UpdateFaceFluxes( gld.faceFluxes, gld.mesh, gld.fieldsRestricted.U, gld.bcData);
    // UpdateFaceFluxesWithMWI( gld.faceFluxes, gld.mesh, gld.fields, gld.fvCoeffs, gld.bcData);
    if constexpr ( LI == Linearisation::Newton ) {
        UpdateFaceAdvectedVelocities( gld.faceAdvectedVelocities, gld.mesh, gld.fvCoeffs, gld.fieldsRestricted.U, gld.faceFluxes, gld.bcData);
        switch ( gld.fvCoeffs.Mom[Axis::X].advectionScheme ) {
            case AdvectionSchemes::Upwind:
                SetIBFaceAdvectedVelocities<AdvectionSchemes::Upwind >( gld.faceAdvectedVelocities, gld.faceFluxes, gld.fieldsRestricted, gld.fvCoeffs, gld.mesh, gld.ibData );
                break;

            case AdvectionSchemes::Central:
                SetIBFaceAdvectedVelocities<AdvectionSchemes::Central>( gld.faceAdvectedVelocities, gld.faceFluxes, gld.fieldsRestricted, gld.fvCoeffs, gld.mesh, gld.ibData );
                break;

            case AdvectionSchemes::SOU:
                SetIBFaceAdvectedVelocities<AdvectionSchemes::SOU    >( gld.faceAdvectedVelocities, gld.faceFluxes, gld.fieldsRestricted, gld.fvCoeffs, gld.mesh, gld.ibData );
                break;

            case AdvectionSchemes::QUICK:
                SetIBFaceAdvectedVelocities<AdvectionSchemes::QUICK  >( gld.faceAdvectedVelocities, gld.faceFluxes, gld.fieldsRestricted, gld.fvCoeffs, gld.mesh, gld.ibData );
                break;
        }
    }
    SetIBFaceFluxes( gld.faceFluxes, gld.ibData );
    SetGhostCells( gridLevelData.fieldsRestricted, gridLevelData.mesh, gridLevelData.bcData );
    UpdateFVCoefficients( gld.fvCoeffs, gld.mesh, gld.fieldsRestricted, gld.fieldsOld, gld.fieldsOldOld, gld.faceAdvectedVelocities, gld.faceFluxes, gld.ibData, gld.bcData);
}



template< MomentumInterpolation MI, Linearisation LI >
void SetCoarseGridEquations( GridLevelData<MI, LI> &gridLevelData,
                             const FieldData<Tensor3D> &coarseGridRightHandSide )
{
    auto &gld = gridLevelData;

    // Set fvCoeffs based on the latest solution to the coarse grid problem
    SetFineGridEquations<MI, LI>( gld );

    // Add the terms that appear on the RHS of the coarse grid equation
    EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {
        gld.fvCoeffs.Mom[axis].F += coarseGridRightHandSide.U[axis];
    } );
    gld.fvCoeffs.Cont.F += coarseGridRightHandSide.P;
}



template< MomentumInterpolation MI, Linearisation LI >
void Smooth( GridLevelData<MI, LI> &gridLevelData,
             const FieldData<floatType> &maxResiduals,
             const intType maxIterations,
             const MultigridEquation mgEquationType )
{
    FieldData<floatType> residualsScaleFactor;
    FieldData<Tensor3D> coarseGridRightHandSide;

    if ( mgEquationType == MultigridEquation::NoTauCorrection ) {
        SetFineGridEquations<MI, LI>( gridLevelData );
    } else {
        SetCoarseGridRightHandSideInStencil( gridLevelData );
        coarseGridRightHandSide = CalculateCoarseGridRightHandSide<MI, LI>( gridLevelData.fvCoeffs,
                                                                            gridLevelData.fieldsRestricted,
                                                                            gridLevelData.residualsRestricted,
                                                                            gridLevelData.ibData.mask );
        SetCoarseGridEquations<MI, LI>( gridLevelData, coarseGridRightHandSide );
    }

    FieldData<floatType > residuals = ScaledL1NormResiduals<MI, LI>( gridLevelData.fields, 
                                                                     gridLevelData.fvCoeffs, 
                                                                     gridLevelData.ibData.mask);
    SetResidualsNormalisationFactor( residualsScaleFactor, residuals );
    NormaliseResiduals( residuals, residualsScaleFactor );

    intType nIterations = 1;
    for ( /* NULL */; nIterations <= maxIterations; nIterations++ ) {

        gridLevelData.linearSolver->UpdateState();
        gridLevelData.linearSolver->Solve();

        if ( mgEquationType == MultigridEquation::NoTauCorrection ) {
            SetFineGridEquations<MI, LI>( gridLevelData );
        } else {
            SetCoarseGridEquations<MI, LI>( gridLevelData, coarseGridRightHandSide );
        }

        residuals = ScaledL1NormResiduals<MI, LI>( gridLevelData.fields,
                                                   gridLevelData.fvCoeffs,
                                                   gridLevelData.ibData.mask );
        NormaliseResiduals( residuals, residualsScaleFactor );

        if ( ResidualsDiverged(residuals) ) {
            break;
        }

        if ( nIterations + 1 > maxIterations ) {
            break;
        }

        if ( MetResidualTolerence(residuals, maxResiduals) ) {
            break;
        }

    }
    std::cout << std::string(gridLevelData.level, ' ') << " Level " << gridLevelData.level << ", relative residuals: " << residuals.U[0] << ", " 
                                                                                                                       << residuals.U[1] << ", " 
                                                                                                                       << residuals.U[2] << ", " 
                                                                                                                       << residuals.P    << ", "
                                                                                                                       << " after " << nIterations << " iterations.\n";
}



template< MomentumInterpolation MI, Linearisation LI >
void SmoothWithFixedIterations( GridLevelData<MI, LI> &gridLevelData,
                                const intType maxIterations,
                                const MultigridEquation mgEquationType )
{
    FieldData<floatType> residualsScaleFactor;
    FieldData<Tensor3D> coarseGridRightHandSide;

    if ( mgEquationType == MultigridEquation::NoTauCorrection ) {
        SetFineGridEquations<MI, LI>( gridLevelData );
    } else {
        SetCoarseGridRightHandSideInStencil( gridLevelData );
        coarseGridRightHandSide = CalculateCoarseGridRightHandSide<MI, LI>( gridLevelData.fvCoeffs,
                                                                            gridLevelData.fieldsRestricted,
                                                                            gridLevelData.residualsRestricted,
                                                                            gridLevelData.ibData.mask );
        SetCoarseGridEquations<MI, LI>( gridLevelData, coarseGridRightHandSide );
    }

    for ( intType nIterations = 1; nIterations <= maxIterations; nIterations++ ) {

        gridLevelData.linearSolver->UpdateState();
        gridLevelData.linearSolver->Solve();

        if ( mgEquationType == MultigridEquation::NoTauCorrection ) {
            SetFineGridEquations<MI, LI>( gridLevelData );
        } else {
            SetCoarseGridEquations<MI, LI>( gridLevelData, coarseGridRightHandSide );
        }
    }
    std::cout << std::string(gridLevelData.level, ' ') << " Level " << gridLevelData.level << ", performed " << maxIterations << " iterations" << "\n";
}



template< MultigridCycleType MGCycle, MomentumInterpolation MI, Linearisation LI >
void Cycle( std::vector< GridLevelData<MI, LI> > &mgLevels,
            const intType level,
            const MultigridEquation mgEquationType,
            const InputData::MultigridSettings &mgSettings )
{
    // Presmoothing 
    SmoothWithFixedIterations<MI, LI>( mgLevels[level], mgSettings.maxPreSmoothingIterations, mgEquationType );

    // Calculate residual
    FieldData<Tensor3D> residuals = ResidualsField<MI, LI>( mgLevels[level].fields, 
                                                            mgLevels[level].fvCoeffs, 
                                                            mgLevels[level].ibData.mask );

    // Restrict residuals and solutions
    mgLevels[level+1].residualsRestricted = RestrictFields( residuals                   , mgLevels[level].mesh, mgLevels[level+1].mesh, mgLevels[level+1].ibData.mask );
    mgLevels[level+1].fieldsRestricted    = RestrictFields( mgLevels[level].fields      , mgLevels[level].mesh, mgLevels[level+1].mesh, mgLevels[level+1].ibData.mask );
    if ( mgLevels[level+1].fieldsOld.P.size() != 0 ) {
        mgLevels[level+1].fieldsOld    = RestrictFields( mgLevels[level].fieldsOld   , mgLevels[level].mesh, mgLevels[level+1].mesh, mgLevels[level+1].ibData.mask );
    }
    if ( mgLevels[level+1].fieldsOldOld.P.size() != 0 ) {
        mgLevels[level+1].fieldsOldOld = RestrictFields( mgLevels[level].fieldsOldOld, mgLevels[level].mesh, mgLevels[level+1].mesh, mgLevels[level+1].ibData.mask );
    }
    mgLevels[level+1].fields = mgLevels[level+1].fieldsRestricted;

    // Solve coarse grid problem
    if ( mgLevels[level+1].isCoarsestLevel ) {  // This is bad if there is 0 coarse levels
        Smooth<MI, LI>( mgLevels[level+1], mgSettings.maxCoarseGridResiduals, mgSettings.maxCoarseGridIterations, MultigridEquation::TauCorrection );
    } else {
        Cycle<MGCycle, MI, LI>( mgLevels, level+1, MultigridEquation::TauCorrection, mgSettings );
    }

    // Compute fine grid correction 
    FieldData<Tensor3D> fineGridCorrection = ComputeFineGridCorrection( mgLevels[level+1].fields, 
                                                                        mgLevels[level+1].fieldsRestricted, 
                                                                        mgLevels[level+1].mesh, 
                                                                        mgLevels[level].mesh,
                                                                        mgLevels[level].ibData.mask );

    // Correct fine grid approximation
    ForAllFieldData( [&] (intType f) {
        mgLevels[level].fields[f] += fineGridCorrection[f];
    } );

    if constexpr ( MGCycle == MultigridCycleType::F ||
                   MGCycle == MultigridCycleType::W  )
    {

        // Re-smoothing
        if ( mgLevels[level].isFinestLevel ) {
            SmoothWithFixedIterations<MI, LI>( mgLevels[level], mgSettings.maxFineGridIterations, mgEquationType );
        } else {
            SmoothWithFixedIterations<MI, LI>( mgLevels[level], mgSettings.maxPostSmoothingIterations, mgEquationType );
        }

        // Calculate residual
        residuals = ResidualsField<MI, LI>( mgLevels[level].fields, 
                                            mgLevels[level].fvCoeffs, 
                                            mgLevels[level].ibData.mask );

        // Restrict residuals and solution
        mgLevels[level+1].residualsRestricted = RestrictFields( residuals                   , mgLevels[level].mesh, mgLevels[level+1].mesh, mgLevels[level+1].ibData.mask );
        mgLevels[level+1].fieldsRestricted    = RestrictFields( mgLevels[level].fields      , mgLevels[level].mesh, mgLevels[level+1].mesh, mgLevels[level+1].ibData.mask );
        if ( mgLevels[level+1].fieldsOld.P.size() != 0 ) {
            mgLevels[level+1].fieldsOld    = RestrictFields( mgLevels[level].fieldsOld   , mgLevels[level].mesh, mgLevels[level+1].mesh, mgLevels[level+1].ibData.mask );
        }
        if ( mgLevels[level+1].fieldsOldOld.P.size() != 0 ) {
            mgLevels[level+1].fieldsOldOld = RestrictFields( mgLevels[level].fieldsOldOld, mgLevels[level].mesh, mgLevels[level+1].mesh, mgLevels[level+1].ibData.mask );
        }
        mgLevels[level+1].fields = mgLevels[level+1].fieldsRestricted;

        // Solve coarse grid problem
        if ( mgLevels[level+1].isCoarsestLevel ) {  // This is bad if there is 0 coarse levels
            Smooth<MI, LI>( mgLevels[level+1], mgSettings.maxCoarseGridResiduals, mgSettings.maxCoarseGridIterations, MultigridEquation::TauCorrection );
        } else {
            if        constexpr ( MGCycle == MultigridCycleType::F ) {
                Cycle<MultigridCycleType::V, MI, LI>( mgLevels, level+1, MultigridEquation::TauCorrection, mgSettings );
            } else if constexpr ( MGCycle == MultigridCycleType::W ) {
                Cycle<MultigridCycleType::W, MI, LI>( mgLevels, level+1, MultigridEquation::TauCorrection, mgSettings );
            }
        }

        // Compute fine grid correction 
        fineGridCorrection = ComputeFineGridCorrection( mgLevels[level+1].fields, 
                                                        mgLevels[level+1].fieldsRestricted, 
                                                        mgLevels[level+1].mesh, 
                                                        mgLevels[level].mesh,
                                                        mgLevels[level].ibData.mask );

        // Correct fine grid approximation
        ForAllFieldData( [&] (intType f) {
            mgLevels[level].fields[f] += fineGridCorrection[f];
        } );

    }

    // Post-smoothing
    if ( mgLevels[level].isFinestLevel ) {
        SmoothWithFixedIterations<MI, LI>( mgLevels[level], mgSettings.maxFineGridIterations, mgEquationType );
    } else {
        SmoothWithFixedIterations<MI, LI>( mgLevels[level], mgSettings.maxPostSmoothingIterations, mgEquationType );
    }

}



// FMG
template< MomentumInterpolation MI, Linearisation LI >
void MultigridCycle( std::vector< GridLevelData<MI, LI> > &mgLevels,
                     const InputData::MultigridSettings &mgSettings,
                     const intType iteration )
{
    if ( iteration == 1 ) {
        
        intType coarsestLevel = mgLevels.size() - 1;

        // Solve coarsest level
        Smooth<MI, LI>( mgLevels[coarsestLevel], mgSettings.maxCoarseGridResiduals, mgSettings.maxCoarseGridIterations, MultigridEquation::NoTauCorrection );

        // Prolongate
        ForAllFieldData( [&] (intType f) {
            mgLevels[coarsestLevel-1].fields[f] = ProlongateField( mgLevels[coarsestLevel].fields[f], mgLevels[coarsestLevel].mesh, mgLevels[coarsestLevel-1].mesh );
        } );
        MaskFields( mgLevels[coarsestLevel-1].fields, mgLevels[coarsestLevel-1].ibData.mask );

        for ( intType level = coarsestLevel-1; level != 0; level-- ) {

            // VCycle
            Cycle<MultigridCycleType::V, MI, LI>( mgLevels, level, MultigridEquation::NoTauCorrection, mgSettings );

            // Prolongate
            ForAllFieldData( [&] (intType f) {
                mgLevels[level-1].fields[f] = ProlongateField( mgLevels[level].fields[f], mgLevels[level].mesh, mgLevels[level-1].mesh );
            } );
            MaskFields( mgLevels[level-1].fields, mgLevels[level-1].ibData.mask );

        }

    } else  {
     
        switch ( mgSettings.cycle ) {
            case MultigridCycleType::V:
                Cycle<MultigridCycleType::V, MI, LI>( mgLevels, 0, MultigridEquation::NoTauCorrection, mgSettings );
                break;
            case MultigridCycleType::F:
                Cycle<MultigridCycleType::F, MI, LI>( mgLevels, 0, MultigridEquation::NoTauCorrection, mgSettings );
                break;
            case MultigridCycleType::W:
                Cycle<MultigridCycleType::W, MI, LI>( mgLevels, 0, MultigridEquation::NoTauCorrection, mgSettings );
                break;
        }

    }
}


}   // end anonymous namespace



template< MomentumInterpolation MI, Linearisation LI >
void SolveSteady( const InputData &inputData,
                  const AxisTransformationMap &axisTransformation )
{
    using enum Axis::ENUMDATA;
    TIC("Pre processing")
    // Extract from input data
    const intType maxOuterIterations = inputData.schemes.maxOuterIterations;
    const FieldData<floatType> maxOuterResiduals = inputData.schemes.maxOuterResiduals;

    // Multigrid level data
    std::vector< GridLevelData<MI, LI> > mgLevels; 
    SetMGLevels( mgLevels, inputData, axisTransformation );

    // References to finest grid
    auto& fields = mgLevels[0].fields;
    auto& mesh   = mgLevels[0].mesh;
    auto& bcData = mgLevels[0].bcData;

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
    
    FieldWriter fieldWriter( fields, mesh, bcData, axisTransformation, inputData.fieldOutputFilename );
    ResidualLogFile residualsLogFile( inputData.residualHistoryFilename, axisTransformation );
    ConsoleLog consoleLog( axisTransformation );

    bool writeFields = ( inputData.fieldWriteInterval > 0 );

    // Calculate initial residual
    residualsOuter   = ScaledL1NormResiduals<MI, LI>( mgLevels[0].fields, 
                                                      mgLevels[0].fvCoeffs, 
                                                      mgLevels[0].ibData.mask);
    SetResidualsNormalisationFactor( residualsScaleFactor, residualsOuter );
    NormaliseResiduals( residualsOuter, residualsScaleFactor );
    massFluxResidual = BoundaryMassFluxResidual( mgLevels[0].faceFluxes, 
                                                 mgLevels[0].mesh);
    TOC()

    TIC("Solver Loop")
    consoleLog.WriteHeader();
    consoleLog.WriteResiduals( residualsOuter, massFluxResidual, 0 );
    residualsLogFile.WriteData( residualsOuter, massFluxResidual, 0 );
    fieldWriter.WriteDataIteration( 0 );
    for ( intType nOuterIterations = 1; nOuterIterations <= maxOuterIterations; nOuterIterations++ )
    {
        TIC("Multigrid Cycling")
        if ( mgLevels.size() == 1 ) {
            SmoothWithFixedIterations<MI, LI>( mgLevels[0], 1, MultigridEquation::NoTauCorrection );
        } else {
            MultigridCycle( mgLevels, inputData.multigridSettings, nOuterIterations );
        }
        TOC()
        
        TIC("Residuals and Logging")
        residualsOuter   = ScaledL1NormResiduals<MI, LI>( mgLevels[0].fields, 
                                                          mgLevels[0].fvCoeffs, 
                                                          mgLevels[0].ibData.mask); 
        NormaliseResiduals( residualsOuter, residualsScaleFactor );
        massFluxResidual = BoundaryMassFluxResidual( mgLevels[0].faceFluxes, 
                                                     mgLevels[0].mesh);

        probeValues      = SetFieldProbeValues( mgLevels[0].fields, fieldProbes); 
        
        consoleLog.WriteResiduals( residualsOuter, massFluxResidual, nOuterIterations );
        residualsLogFile.WriteData( residualsOuter, massFluxResidual, nOuterIterations );
        for ( size_t p = 0; p != fieldProbes.size(); p++ ) {
            probeLogFiles[p].WriteData( probeValues[p], nOuterIterations );
        }
        TOC() 
        
        if ( ResidualsDiverged(residualsOuter) ) {
            fieldWriter.WriteDataIteration( nOuterIterations );
            std::cout << "*** SOLUTION DIVERGED ***" << "\n\n";
            break;
        }

        if ( nOuterIterations + 1 > maxOuterIterations ) {
            fieldWriter.WriteDataIteration( nOuterIterations );
            std::cout << "*** REACHED ITERATION LIMIT ***" << "\n\n";
            break;
        }

        if ( MetResidualTolerence(residualsOuter, maxOuterResiduals) ) {
            fieldWriter.WriteDataIteration( nOuterIterations );
            std::cout << "*** SOLUTION CONVERGED ***" << "\n\n";
            break;
        }

        if ( writeFields && (nOuterIterations % inputData.fieldWriteInterval) == 0 ) {
            fieldWriter.WriteDataIteration( nOuterIterations );
        }  

    }
    TOC()

}
template void SolveSteady<MomentumInterpolation::Implicit    , Linearisation::Picard>( const InputData &, const AxisTransformationMap &);
template void SolveSteady<MomentumInterpolation::SemiExplicit, Linearisation::Picard>( const InputData &, const AxisTransformationMap &);
template void SolveSteady<MomentumInterpolation::Implicit    , Linearisation::Newton>( const InputData &, const AxisTransformationMap &);
template void SolveSteady<MomentumInterpolation::SemiExplicit, Linearisation::Newton>( const InputData &, const AxisTransformationMap &);



template< MomentumInterpolation MI, Linearisation LI >
void SolveTransient( const InputData &inputData,
                     const AxisTransformationMap &axisTransformation )
{
    using enum Axis::ENUMDATA;
    TIC("Pre processing")
    // Extract from input data
    const intType maxOuterIterations = inputData.schemes.maxOuterIterations;
    const FieldData<floatType> maxOuterResiduals = inputData.schemes.maxOuterResiduals;

    // Multigrid level data
    std::vector< GridLevelData<MI, LI> > mgLevels; 
    SetMGLevels( mgLevels, inputData, axisTransformation );

    // References to finest grid
    auto& fields = mgLevels[0].fields;
    auto& mesh   = mgLevels[0].mesh;
    auto& bcData = mgLevels[0].bcData;

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
    
    FieldWriter fieldWriter( fields, mesh, bcData, axisTransformation, inputData.fieldOutputFilename );
    ResidualLogFile residualsLogFile( inputData.residualHistoryFilename, axisTransformation );
    ConsoleLog consoleLog( axisTransformation );

    bool writeFields = ( inputData.fieldWriteInterval > 0 );

    // Calculate initial residual
    residualsOuter   = ScaledL1NormResiduals<MI, LI>( mgLevels[0].fields, 
                                                      mgLevels[0].fvCoeffs, 
                                                      mgLevels[0].ibData.mask);
    SetResidualsNormalisationFactor( residualsScaleFactor, residualsOuter );
    NormaliseResiduals( residualsOuter, residualsScaleFactor );
    massFluxResidual = BoundaryMassFluxResidual( mgLevels[0].faceFluxes, 
                                                 mgLevels[0].mesh);
    TOC()

    TIC("Solver Loop")
    consoleLog.WriteHeader();
    residualsLogFile.WriteData( residualsOuter, massFluxResidual, 0 );
    consoleLog.WriteResiduals( residualsOuter, massFluxResidual, 0 );
    // fieldWriter.WriteDataTime( 0.0f );
    fieldWriter.WriteDataIteration( 0 );
    bool abortTimestepping = false;

    for ( intType timeStepNumber = 1; timeStepNumber < inputData.schemes.numberOfTimesteps; timeStepNumber++ ) {
        
        floatType currentTime = static_cast<floatType>(timeStepNumber) * inputData.schemes.timeStep;
        std::cout << "\n" << "Time = " << currentTime << "\n";

        for ( intType nOuterIterations = 1; nOuterIterations <= maxOuterIterations; nOuterIterations++ )
        {
            TIC("Multigrid Cycling")
            if ( mgLevels.size() == 1 ) {
                SmoothWithFixedIterations<MI, LI>( mgLevels[0], 1, MultigridEquation::NoTauCorrection );
            } else {
                MultigridCycle( mgLevels, inputData.multigridSettings, nOuterIterations );
            }
            TOC()
            
            TIC("Residuals and Logging")
            residualsOuter   = ScaledL1NormResiduals<MI, LI>( mgLevels[0].fields, 
                                                              mgLevels[0].fvCoeffs, 
                                                              mgLevels[0].ibData.mask); 
            NormaliseResiduals( residualsOuter, residualsScaleFactor );
            massFluxResidual = BoundaryMassFluxResidual( mgLevels[0].faceFluxes, 
                                                        mgLevels[0].mesh);
            consoleLog.WriteResiduals( residualsOuter, massFluxResidual, nOuterIterations );
            residualsLogFile.WriteData( residualsOuter, massFluxResidual, 0 );
            TOC() 
            
            if ( ResidualsDiverged(residualsOuter) ) {
                abortTimestepping = true;
                std::cout << "Solution Diverged for timestep" << "\n";
                break;
            }

            if ( nOuterIterations + 1 > maxOuterIterations ) {
                std::cout << "Reached iteration limit for timestep." << "\n";
                break;
            }

            if ( MetResidualTolerence(residualsOuter, maxOuterResiduals) ) {
                std::cout << "Residuals converged for timestep." << "\n";
                break;
            }

        }

        probeValues = SetFieldProbeValues( mgLevels[0].fields, fieldProbes); 
        for ( size_t p = 0; p != fieldProbes.size(); p++ ) {
            probeLogFiles[p].WriteData( probeValues[p], timeStepNumber );
        }
        residualsLogFile.WriteData( residualsOuter, massFluxResidual, timeStepNumber );

        if ( writeFields && (timeStepNumber % inputData.fieldWriteInterval) == 0 ) {
            // fieldWriter.WriteDataTime( currentTime );
            fieldWriter.WriteDataIteration( timeStepNumber );
        }  

        if ( abortTimestepping ) {
            std::cout << "*** TIMESTEPPING ABORTED ***" << "\n\n";
            break;
        }

        mgLevels[0].fieldsOldOld = mgLevels[0].fieldsOld;
        mgLevels[0].fieldsOld    = mgLevels[0].fields;

    }


    TOC()

}
template void SolveTransient<MomentumInterpolation::Implicit    , Linearisation::Picard>( const InputData &, const AxisTransformationMap &);
template void SolveTransient<MomentumInterpolation::SemiExplicit, Linearisation::Picard>( const InputData &, const AxisTransformationMap &);
template void SolveTransient<MomentumInterpolation::Implicit    , Linearisation::Newton>( const InputData &, const AxisTransformationMap &);
template void SolveTransient<MomentumInterpolation::SemiExplicit, Linearisation::Newton>( const InputData &, const AxisTransformationMap &);



} // end namespace CFD