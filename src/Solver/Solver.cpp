#include "Solver.h"
#include "ConvergenceLogging.h"

#include "../Core/Types.h"
#include "../Core/Macros.h"
#include "../IO/InputProcessing.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"
#include "../DerivedQuantities/DerivedQuantities.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../FiniteVolume/GhostCells.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../Multigrid/Multigrid.h"

#include "../IO/ArrayIO.h"
#include "../IO/VTKWriter.h"

#include "ResidualFunctions.h"

#include <iostream>
#include <memory>

namespace CFD
{

namespace
{

// To be returned from smoother functions
enum class OperationStatus {
    Sucess, Failure
};


template< MomentumInterpolation MI>
void SetStencil( GridLevelData<MI> &gridLevelData,
                 FieldData<Tensor3D> &fields )
{
    auto &gld = gridLevelData; 

    TIC("Updating Coefficients")
    SetGhostCells( fields, gridLevelData.mesh, gridLevelData.bcData );
    UpdateIBData( gld.ibData, fields );
    UpdateFaceFluxes( gld.faceFluxes, gld.mesh, fields.U, gld.bcData);
    SetIBFaceFluxes( gld.faceFluxes, gld.ibData );
    UpdateFVCoefficients( gld.fvCoeffs, gld.mesh, fields, gld.fieldsPrevTime, gld.fieldsPrevPrevTime, gld.faceFluxes, gld.ibData);
    TOC()
}


template<MomentumInterpolation MI>
void SetCoarseGridRightHandSide( GridLevelData<MI> &gridLevelData,
                                 const MultigridEquation mgEquationType )
{
    if ( mgEquationType == MultigridEquation::NoTauCorrection )
        return;

    // Set the coarse grid right hand side in the stencil coefficients
    SetStencil<MI>(gridLevelData, gridLevelData.fieldsRestricted);

    CalculateCoarseGridRightHandSide<MI >( gridLevelData.coarseGridRightHandSide,
                                           gridLevelData.fvCoeffs,
                                           gridLevelData.fieldsRestricted,
                                           gridLevelData.residualsRestricted,
                                           gridLevelData.ibData.mask );
}


template< MomentumInterpolation MI >
void SetFineGridEquations( GridLevelData<MI> &gridLevelData )
{
    SetStencil<MI>(gridLevelData, gridLevelData.fields);
}


template< MomentumInterpolation MI >
void SetCoarseGridEquations( GridLevelData<MI> &gridLevelData )
{
    auto &gld = gridLevelData;

    // Set fvCoeffs based on the latest solution to the coarse grid problem
    SetStencil<MI>(gridLevelData, gridLevelData.fields);

    // Add the terms that appear on the RHS of the coarse grid equation
    EnumFor<Axis>( [&] ( Axis::ENUMDATA axis ) {
        gld.fvCoeffs.Mom[axis].F += gridLevelData.coarseGridRightHandSide.U[axis];
    } );
    gld.fvCoeffs.Cont.F += gridLevelData.coarseGridRightHandSide.P; 
}



template< MomentumInterpolation MI >
OperationStatus Smooth( GridLevelData<MI > &gridLevelData,
                        const FieldData<floatType> &maxResiduals,
                        const intType maxIterations,
                        const MultigridEquation mgEquationType )
{
    OperationStatus returnFlag = OperationStatus::Sucess;
    FieldData<floatType> residualsScaleFactor;

    FieldData<floatType > residuals = ScaledL1NormResiduals<MI >( gridLevelData.fields, 
                                                                     gridLevelData.fvCoeffs, 
                                                                     gridLevelData.ibData.mask);
    SetResidualsNormalisationFactor( residualsScaleFactor, residuals );
    NormaliseResiduals( residuals, residualsScaleFactor );

    intType nIterations = 0;
    while ( nIterations < maxIterations ) {
        nIterations++;

        gridLevelData.fieldsOld = gridLevelData.fields;
        if ( mgEquationType == MultigridEquation::NoTauCorrection ) {
            SetFineGridEquations<MI>( gridLevelData );
        } else {
            SetCoarseGridEquations<MI>( gridLevelData );
        }

        residuals = ScaledL1NormResiduals<MI>( gridLevelData.fields,
                                               gridLevelData.fvCoeffs,
                                               gridLevelData.ibData.mask );
        NormaliseResiduals( residuals, residualsScaleFactor );

        if ( ResidualsDiverged(residuals) ) {
            returnFlag = OperationStatus::Failure;
            break;
        }

        if ( MetResidualTolerence(residuals, maxResiduals) ) {
            returnFlag = OperationStatus::Sucess;
            break;
        }
        
        gridLevelData.linearSolver->UpdateState();
        gridLevelData.linearSolver->Solve();
    }
    std::cout << std::string(gridLevelData.level, ' ') << " Level " << gridLevelData.level << ", relative residuals: " << residuals.U[0] << ", " 
                                                                                                                       << residuals.U[1] << ", " 
                                                                                                                       << residuals.U[2] << ", " 
                                                                                                                       << residuals.P    << ", "
                                                                                                                       << " after " << nIterations << " iterations.";
    if ( returnFlag == OperationStatus::Failure ) {
        std::cout << "  (Solution discarded)\n";
    } else {
        std::cout << "\n";
    }
    return returnFlag;
}



template< MomentumInterpolation MI >
OperationStatus SmoothWithFixedIterations( GridLevelData<MI> &gridLevelData,
                                           const intType maxIterations,
                                           const MultigridEquation mgEquationType )
{
    for ( intType nIterations = 1; nIterations <= maxIterations; nIterations++ ) {

        gridLevelData.fieldsOld = gridLevelData.fields;
        if ( mgEquationType == MultigridEquation::NoTauCorrection ) {
            SetFineGridEquations<MI>( gridLevelData );
        } else {
            SetCoarseGridEquations<MI>( gridLevelData );
        }

        gridLevelData.linearSolver->UpdateState();
        gridLevelData.linearSolver->Solve();

    }
    std::cout << std::string(gridLevelData.level, ' ') << " Level " << gridLevelData.level << ", performed " << maxIterations << " iterations" << "\n";

    return OperationStatus::Sucess;
}



template< MomentumInterpolation MI >
void RestrictLevel( std::vector< GridLevelData<MI> > &mgLevels,
                    const size_t fineLevel,
                    const size_t coarseLevel )
{
    RestrictFields( mgLevels[coarseLevel].residualsRestricted, mgLevels[fineLevel].residuals, mgLevels[fineLevel].mesh, mgLevels[coarseLevel].mesh, mgLevels[coarseLevel].ibData.mask );
    RestrictFields( mgLevels[coarseLevel].fieldsRestricted   , mgLevels[fineLevel].fields   , mgLevels[fineLevel].mesh, mgLevels[coarseLevel].mesh, mgLevels[coarseLevel].ibData.mask );
    if ( mgLevels[coarseLevel].fieldsPrevTime.P.size() != 0 ) {
        RestrictFields( mgLevels[coarseLevel].fieldsPrevTime, mgLevels[fineLevel].fieldsPrevTime     , mgLevels[fineLevel].mesh, mgLevels[coarseLevel].mesh, mgLevels[coarseLevel].ibData.mask );
    }
    if ( mgLevels[coarseLevel].fieldsPrevPrevTime.P.size() != 0 ) {
        RestrictFields( mgLevels[coarseLevel].fieldsPrevPrevTime, mgLevels[fineLevel].fieldsPrevPrevTime, mgLevels[fineLevel].mesh, mgLevels[coarseLevel].mesh, mgLevels[coarseLevel].ibData.mask );
    }
    mgLevels[coarseLevel].fields = mgLevels[coarseLevel].fieldsRestricted;
}




template< MultigridCycleType MGCycle, MomentumInterpolation MI >
void Cycle( std::vector< GridLevelData<MI> > &mgLevels,
            const size_t level,
            const MultigridEquation mgEquationType,
            const InputData::MultigridSettings &mgSettings )
{
    OperationStatus smootherStatus = OperationStatus::Sucess;

    // Coarse grid right hand side is constant for this level
    SetCoarseGridRightHandSide( mgLevels[level], mgEquationType );

    // Presmoothing 
    SmoothWithFixedIterations<MI>( mgLevels[level], mgSettings.preSmoothingIterations, mgEquationType );

    // Calculate residual
    ResidualsField<MI>( mgLevels[level].residuals,
                        mgLevels[level].fields, 
                        mgLevels[level].fvCoeffs, 
                        mgLevels[level].ibData.mask );

    // Restrict residuals and solutions
    RestrictLevel(mgLevels, level, level+1);


    // Solve coarse grid problem
    if ( mgLevels[level+1].isCoarsestLevel ) {  // This is bad if there is 0 coarse levels
        SetCoarseGridRightHandSide( mgLevels[level+1], MultigridEquation::TauCorrection );
        smootherStatus = Smooth<MI>( mgLevels[level+1], mgSettings.maxCoarseGridResiduals, mgSettings.maxCoarseGridIterations, MultigridEquation::TauCorrection );
    } else {
        Cycle<MGCycle, MI>( mgLevels, level+1, MultigridEquation::TauCorrection, mgSettings );
    }


    if ( smootherStatus == OperationStatus::Sucess) {
        ComputeFineGridCorrection( mgLevels[level].fineGridCorrection,
                                   mgLevels[level+1].fields, 
                                   mgLevels[level+1].fieldsRestricted, 
                                   mgLevels[level+1].mesh, 
                                   mgLevels[level].mesh,
                                   mgLevels[level].ibData.mask );

        // Correct fine grid approximation
        ForAllFieldData( [&] (intType f) {
            mgLevels[level].fields[f] += mgLevels[level].fineGridCorrection[f];
        } );
    }

    if constexpr ( MGCycle == MultigridCycleType::F ||
                   MGCycle == MultigridCycleType::W  )
    {

        // Re-smoothing
        if ( mgLevels[level].isFinestLevel ) {
            SmoothWithFixedIterations<MI>( mgLevels[level], mgSettings.fineGridIterations, mgEquationType );
        } else {
            SmoothWithFixedIterations<MI>( mgLevels[level], mgSettings.postSmoothingIterations, mgEquationType );
        }

        // Calculate residual
        ResidualsField<MI>( mgLevels[level].residuals,
                            mgLevels[level].fields, 
                            mgLevels[level].fvCoeffs, 
                            mgLevels[level].ibData.mask );

        // Restrict residuals and solution
        RestrictLevel(mgLevels, level, level+1);

        // Solve coarse grid problem
        smootherStatus = OperationStatus::Sucess;
        if ( mgLevels[level+1].isCoarsestLevel ) {  // This is bad if there is 0 coarse levels
            SetCoarseGridRightHandSide( mgLevels[level+1], MultigridEquation::TauCorrection );
            smootherStatus = Smooth<MI>( mgLevels[level+1], mgSettings.maxCoarseGridResiduals, mgSettings.maxCoarseGridIterations, MultigridEquation::TauCorrection );
        } else {
            if        constexpr ( MGCycle == MultigridCycleType::F ) {
                Cycle<MultigridCycleType::V, MI>( mgLevels, level+1, MultigridEquation::TauCorrection, mgSettings );
            } else if constexpr ( MGCycle == MultigridCycleType::W ) {
                Cycle<MultigridCycleType::W, MI>( mgLevels, level+1, MultigridEquation::TauCorrection, mgSettings );
            }
        }

        // Compute fine grid correction 
        if ( smootherStatus == OperationStatus::Sucess) {
            ComputeFineGridCorrection( mgLevels[level].fineGridCorrection,
                                       mgLevels[level+1].fields, 
                                       mgLevels[level+1].fieldsRestricted, 
                                       mgLevels[level+1].mesh, 
                                       mgLevels[level].mesh,
                                       mgLevels[level].ibData.mask );

            // Correct fine grid approximation
            ForAllFieldData( [&] (intType f) {
                mgLevels[level].fields[f] += mgLevels[level].fineGridCorrection[f];
            } );
        }

    }

    // Post-smoothing
    if ( mgLevels[level].isFinestLevel ) {
        SmoothWithFixedIterations<MI>( mgLevels[level], mgSettings.fineGridIterations, mgEquationType );
    } else {
        SmoothWithFixedIterations<MI>( mgLevels[level], mgSettings.postSmoothingIterations, mgEquationType );
    }

}



template< MomentumInterpolation MI >
void FMGInitialise( std::vector< GridLevelData<MI> > &mgLevels,
                    const InputData::MultigridSettings &mgSettings)
{
    size_t coarsestLevel = mgLevels.size() - 1;

    // Solve coarsest level
    OperationStatus smootherStatus = Smooth<MI>( mgLevels[coarsestLevel], mgSettings.maxCoarseGridResiduals, mgSettings.maxCoarseGridIterations, MultigridEquation::NoTauCorrection );

    // Prolongate
    if ( smootherStatus == OperationStatus::Sucess ) {
        ForAllFieldData( [&] (intType f) {
            ProlongateField(mgLevels[coarsestLevel-1].fields[f], mgLevels[coarsestLevel].fields[f], mgLevels[coarsestLevel].mesh, mgLevels[coarsestLevel-1].mesh );
        } );
        MaskFields( mgLevels[coarsestLevel-1].fields, mgLevels[coarsestLevel-1].ibData.mask );
    }
    

    for ( size_t level = coarsestLevel-1; level != 0; level-- ) {

        // VCycle
        Cycle<MultigridCycleType::V, MI>( mgLevels, level, MultigridEquation::NoTauCorrection, mgSettings );

        // Prolongate
        ForAllFieldData( [&] (intType f) {
            ProlongateField(mgLevels[level-1].fields[f], mgLevels[level].fields[f], mgLevels[level].mesh, mgLevels[level-1].mesh );
        } );
        MaskFields( mgLevels[level-1].fields, mgLevels[level-1].ibData.mask );

    }
}



template< MomentumInterpolation MI >
void MultigridCycle( std::vector< GridLevelData<MI> > &mgLevels,
                     const InputData::MultigridSettings &mgSettings)
{
    switch ( mgSettings.cycle ) {
        case MultigridCycleType::V:
            Cycle<MultigridCycleType::V, MI>( mgLevels, 0, MultigridEquation::NoTauCorrection, mgSettings );
            break;
        case MultigridCycleType::F:
            Cycle<MultigridCycleType::F, MI>( mgLevels, 0, MultigridEquation::NoTauCorrection, mgSettings );
            break;
        case MultigridCycleType::W:
            Cycle<MultigridCycleType::W, MI>( mgLevels, 0, MultigridEquation::NoTauCorrection, mgSettings );
            break;
    }
}


}   // end anonymous namespace



template< MomentumInterpolation MI >
void SolveSteady( const InputData &inputData,
                  const AxisTransformationMap &axisTransformation )
{
    using enum Axis::ENUMDATA;
    TIC("Pre processing")
    // Extract from input data
    const intType maxOuterIterations = inputData.schemes.maxOuterIterations;
    const FieldData<floatType> maxOuterResiduals = inputData.schemes.maxOuterResiduals;

    // Multigrid level data
    std::vector< GridLevelData<MI> > mgLevels; 
    SetMGLevels( mgLevels, inputData, axisTransformation );

    // References to finest grid
    auto& fields = mgLevels[0].fields;
    auto& mesh   = mgLevels[0].mesh;
    auto& bcData = mgLevels[0].bcData;
    auto& ibData = mgLevels[0].ibData;
    auto& mask   = mgLevels[0].ibData.mask;

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
    
    ForceCalculator forceCalculator( ibData, mesh, fields, inputData.rho, inputData.nu );
    std::unique_ptr<ForceLogFile> forceLogFilePtr;
    if ( inputData.calculateForces )
        forceLogFilePtr = std::make_unique<ForceLogFile>( inputData.forceCalculatorFilename, axisTransformation );

    FieldWriter fieldWriter( fields, mask, mesh, bcData, axisTransformation, inputData.fieldOutputFilename );
    ResidualLogFile residualsLogFile( inputData.residualHistoryFilename, axisTransformation );
    ConsoleLog consoleLog( axisTransformation );

    bool writeFields = ( inputData.fieldWriteInterval > 0 );

    // Calculate initial residual
    residualsOuter   = ScaledL1NormResiduals<MI>( mgLevels[0].fields, 
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

    intType nOuterIterations = 0;
    while ( nOuterIterations < maxOuterIterations ) {
        nOuterIterations++;

        // Solve
        TIC("Multigrid Cycling")
        if ( mgLevels.size() == 1 ) {
            SmoothWithFixedIterations<MI>( mgLevels[0], 1, MultigridEquation::NoTauCorrection );
        } else {
            if ( nOuterIterations == 1 )
                FMGInitialise(mgLevels, inputData.multigridSettings );
            MultigridCycle( mgLevels, inputData.multigridSettings );
        }
        TOC()
        
        TIC("Residuals and Logging")
        // Residuals
        SetStencil<MI>(mgLevels[0], mgLevels[0].fields);
        
        residualsOuter   = ScaledL1NormResiduals<MI>( mgLevels[0].fields, 
                                                      mgLevels[0].fvCoeffs, 
                                                      mgLevels[0].ibData.mask); 
        NormaliseResiduals( residualsOuter, residualsScaleFactor );
        massFluxResidual = BoundaryMassFluxResidual( mgLevels[0].faceFluxes, 
                                                     mgLevels[0].mesh);
        consoleLog.WriteResiduals( residualsOuter, massFluxResidual, nOuterIterations );
        residualsLogFile.WriteData( residualsOuter, massFluxResidual, nOuterIterations );

        // Probes
        for ( size_t p = 0; p != fieldProbes.size(); p++ ) {
            probeLogFiles[p].WriteData( ProbeAllFieldValues( mgLevels[0].fields, fieldProbes[p] ), 
                                        nOuterIterations );
        }
        
        // Forces
        if ( forceLogFilePtr ) {
            forceLogFilePtr->WriteData( forceCalculator.GetForce(), nOuterIterations );
        }
        TOC() 
        
        // Stopping criteria
        if ( ResidualsDiverged(residualsOuter) ) {
            fieldWriter.WriteDataIteration( nOuterIterations );
            std::cout << "*** SOLUTION DIVERGED ***" << "\n\n";
            TOC()
            return;
        }

        if ( MetResidualTolerence(residualsOuter, maxOuterResiduals) ) {
            fieldWriter.WriteDataIteration( nOuterIterations );
            std::cout << "*** SOLUTION CONVERGED ***" << "\n\n";
            TOC()
            return;
        }

        if ( writeFields && (nOuterIterations % inputData.fieldWriteInterval) == 0 ) {
            fieldWriter.WriteDataIteration( nOuterIterations );
        }  

    }

    fieldWriter.WriteDataIteration( nOuterIterations );
    std::cout << "*** REACHED ITERATION LIMIT ***" << "\n\n";
    TOC()

}
template void SolveSteady<MomentumInterpolation::Implicit     >( const InputData &, const AxisTransformationMap &);
template void SolveSteady<MomentumInterpolation::SemiExplicit >( const InputData &, const AxisTransformationMap &);





template< MomentumInterpolation MI >
void SolveTransient( const InputData &inputData,
                     const AxisTransformationMap &axisTransformation )
{
    using enum Axis::ENUMDATA;
    TIC("Pre processing")
    // Extract from input data
    const intType maxOuterIterations = inputData.schemes.maxOuterIterations;
    const FieldData<floatType> maxOuterResiduals = inputData.schemes.maxOuterResiduals;

    // Multigrid level data
    std::vector< GridLevelData<MI> > mgLevels; 
    SetMGLevels( mgLevels, inputData, axisTransformation );

    // References to finest grid
    auto& fields = mgLevels[0].fields;
    auto& mesh   = mgLevels[0].mesh;
    auto& bcData = mgLevels[0].bcData;
    auto& ibData = mgLevels[0].ibData;
    auto& mask   = mgLevels[0].ibData.mask;

    // Initialise residuals
    FieldData<floatType> residualsOuter(0.0f), residualsScaleFactor(1.0f);
    floatType massFluxResidual(0.0f);

    // Logging objects
    std::vector< FieldProbe > fieldProbes;
    std::vector< ProbeLogFile > probeLogFiles;
    for ( const auto &probeData : inputData.probes ) {
        fieldProbes.emplace_back( mesh, probeData.location );
        probeLogFiles.emplace_back( probeData.filename, axisTransformation, fieldProbes.back() );
    }
    
    ForceCalculator forceCalculator( ibData, mesh, fields, inputData.rho, inputData.nu );
    std::unique_ptr<ForceLogFile> forceLogFilePtr;
    if ( inputData.calculateForces )
        forceLogFilePtr = std::make_unique<ForceLogFile>( inputData.forceCalculatorFilename, axisTransformation );


    FieldWriter fieldWriter( fields, mask, mesh, bcData, axisTransformation, inputData.fieldOutputFilename );
    ResidualLogFile residualsLogFile( inputData.residualHistoryFilename, axisTransformation );
    ConsoleLog consoleLog( axisTransformation );

    bool writeFields = ( inputData.fieldWriteInterval > 0 );
    TOC()

    TIC("Solver Loop")
    consoleLog.WriteHeader();
    fieldWriter.WriteDataIteration( 0 );
    bool abortTimestepping = false;

    for ( intType timeStepNumber = 1; timeStepNumber < inputData.schemes.numberOfTimesteps; timeStepNumber++ ) {
        
        floatType currentTime = static_cast<floatType>(timeStepNumber) * inputData.schemes.timeStep;
        std::cout << "\n" << "Time = " << currentTime << ", timestep no. " << timeStepNumber << "\n";

        for ( intType nOuterIterations = 1; nOuterIterations <= maxOuterIterations; nOuterIterations++ )
        {
            TIC("Multigrid Cycling")
            if ( mgLevels.size() == 1 ) {
                SmoothWithFixedIterations<MI>( mgLevels[0], 1, MultigridEquation::NoTauCorrection );
            } else {
                MultigridCycle( mgLevels, inputData.multigridSettings );
            }
            TOC()
            
            SetStencil<MI>(mgLevels[0], mgLevels[0].fields);
            residualsOuter   = ScaledL1NormResiduals<MI>( mgLevels[0].fields, 
                                                          mgLevels[0].fvCoeffs, 
                                                          mgLevels[0].ibData.mask); 
            // if ( nOuterIterations == 1 )
            //     SetResidualsNormalisationFactor( residualsScaleFactor, residualsOuter );
            NormaliseResiduals( residualsOuter, residualsScaleFactor );
            massFluxResidual = BoundaryMassFluxResidual( mgLevels[0].faceFluxes, 
                                                         mgLevels[0].mesh);
            consoleLog.WriteResiduals( residualsOuter, massFluxResidual, nOuterIterations );
            
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

        for ( size_t p = 0; p != fieldProbes.size(); p++ ) {
            probeLogFiles[p].WriteData( ProbeAllFieldValues( mgLevels[0].fields, fieldProbes[p] ), 
                                        timeStepNumber );
        }

        if ( forceLogFilePtr ) {
            forceLogFilePtr->WriteData( forceCalculator.GetForce(), timeStepNumber );
        }

        residualsLogFile.WriteData( residualsOuter, massFluxResidual, timeStepNumber );

        if ( writeFields && (timeStepNumber % inputData.fieldWriteInterval) == 0 ) {
            fieldWriter.WriteDataIteration( timeStepNumber );
        }  

        if ( abortTimestepping ) {
            std::cout << "*** TIMESTEPPING ABORTED ***" << "\n\n";
            break;
        }

        mgLevels[0].fieldsPrevPrevTime = mgLevels[0].fieldsPrevTime;
        mgLevels[0].fieldsPrevTime    = mgLevels[0].fields;

    }


    TOC()

}
template void SolveTransient<MomentumInterpolation::Implicit    >( const InputData &, const AxisTransformationMap &);
template void SolveTransient<MomentumInterpolation::SemiExplicit>( const InputData &, const AxisTransformationMap &);



} // end namespace CFD