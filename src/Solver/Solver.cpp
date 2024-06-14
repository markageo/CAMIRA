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

#include "TriadSolver.h"
#include "LineSolver.h"
#include "PlaneSolver.h"
#include "LinearSolver.h"
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
    if constexpr ( LI == Linearisation::Newton ) {
        UpdateFaceAdvectedVelocities( gld.faceAdvectedVelocities, gld.mesh, gld.fields.U, gld.faceFluxes, gld.bcData);
    }
    SetIBFaceFluxes( gld.faceFluxes, gld.ibData );
    UpdateFVCoefficients( gld.fvCoeffs, gld.mesh, gld.fields, gld.faceAdvectedVelocities, gld.faceFluxes, gld.ibData, gld.bcData);
}


template< MomentumInterpolation MI, Linearisation LI >
void SetCoarseGridEquations( GridLevelData<MI, LI> &gridLevelData )
{
    auto &gld = gridLevelData;

    // First set fvCoeffs based on the restricted fine grid approximation for the RHS
    UpdateIBData( gld.ibData, gld.fieldsRestricted );
    UpdateFaceFluxes( gld.faceFluxes, gld.mesh, gld.fieldsRestricted.U, gld.bcData);
    if constexpr ( LI == Linearisation::Newton ) {
        UpdateFaceAdvectedVelocities( gld.faceAdvectedVelocities, gld.mesh, gld.fieldsRestricted.U, gld.faceFluxes, gld.bcData);
    }
    SetIBFaceFluxes( gld.faceFluxes, gld.ibData );
    UpdateFVCoefficients( gld.fvCoeffs, gld.mesh, gld.fieldsRestricted, gld.faceAdvectedVelocities, gld.faceFluxes, gld.ibData, gld.bcData);


    // Calculate the extra terms that appear on the RHS of the coarse grid equation
    FieldData<Tensor3D> coarseGridRightHandSide = CalculateCoarseGridRightHandSide<MI, LI>( gld.fvCoeffs,
                                                                                            gld.fieldsRestricted,
                                                                                            gld.residualsRestricted,
                                                                                            gld.ibData.mask );


    // Now set fvCoeffs based on the latest solution to the coarse grid problem
    UpdateIBData( gld.ibData, gld.fields );
    UpdateFaceFluxes( gld.faceFluxes, gld.mesh, gld.fields.U, gld.bcData);
    if constexpr ( LI == Linearisation::Newton ) {
        UpdateFaceAdvectedVelocities( gld.faceAdvectedVelocities, gld.mesh, gld.fields.U, gld.faceFluxes, gld.bcData);
    }
    SetIBFaceFluxes( gld.faceFluxes, gld.ibData );
    UpdateFVCoefficients( gld.fvCoeffs, gld.mesh, gld.fields, gld.faceAdvectedVelocities, gld.faceFluxes, gld.ibData, gld.bcData);


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

    if ( mgEquationType == MultigridEquation::NoTauCorrection ) {
        SetFineGridEquations<MI, LI>( gridLevelData );
    } else {
        SetCoarseGridEquations<MI, LI>( gridLevelData );
    }

    FieldData<floatType > residuals = ScaledL1NormResiduals<MI, LI>( gridLevelData.fields, 
                                                                     gridLevelData.fvCoeffs, 
                                                                     gridLevelData.ibData.mask);
    SetResidualsNormalisationFactor( residualsScaleFactor, residuals );
    NormaliseResiduals( residuals, residualsScaleFactor );

    for ( intType nIterations = 1; nIterations <= maxIterations; nIterations++ ) {

        gridLevelData.linearSolver->UpdateState();
        gridLevelData.linearSolver->Solve();
        // SetGhostCells( gridLevelData.fields, gridLevelData.mesh, gridLevelData.bcData ); // Deferred correction needs this

        if ( mgEquationType == MultigridEquation::NoTauCorrection ) {
            SetFineGridEquations<MI, LI>( gridLevelData );
        } else {
            SetCoarseGridEquations<MI, LI>( gridLevelData );
        }

        residuals = ScaledL1NormResiduals<MI, LI>( gridLevelData.fields,
                                                   gridLevelData.fvCoeffs,
                                                   gridLevelData.ibData.mask );
        NormaliseResiduals( residuals, residualsScaleFactor );

        gridLevelData.fieldsOld = gridLevelData.fields;

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
                                                                                                                       << residuals.P << "\n";

}



template< MomentumInterpolation MI, Linearisation LI >
void VCycle( std::vector< GridLevelData<MI, LI> > &mgLevels,
             const intType startingLevel,
             const intType level, 
             const InputData::MultigridSettings &mgSettings,
             const intType iteration,
             const AxisTransformationMap &axisTransformation )
{
    MultigridEquation mgEquationType = ( level == startingLevel ) ? MultigridEquation::NoTauCorrection : MultigridEquation::TauCorrection;

    // Presmoothing 
    TIC("Presmoothing")
    Smooth<MI, LI>( mgLevels[level], mgSettings.maxPreSmoothingResiduals, mgSettings.maxPreSmoothingIterations, mgEquationType );
    TOC()

    // Calculate residual
    TIC("Residual calculation")
    FieldData<Tensor3D> residuals = ResidualsField<MI, LI>( mgLevels[level].fields, 
                                                            mgLevels[level].fvCoeffs, 
                                                            mgLevels[level].ibData.mask );
    TOC()

    // Restrict residual
    TIC("Residual restriction")
    ForAllFieldData( [&] (intType f) {
        mgLevels[level+1].residualsRestricted[f] = RestrictField( residuals[f],
                                                                  mgLevels[level].mesh, 
                                                                  mgLevels[level+1].mesh );
    } );
    TOC()

    // Restrict solution
    TIC("Solution restriction")
    ForAllFieldData( [&] (intType f) {
        mgLevels[level+1].fieldsRestricted[f] = RestrictField( mgLevels[level].fields[f], 
                                                               mgLevels[level].mesh, 
                                                               mgLevels[level+1].mesh );
        mgLevels[level+1].fields[f] = mgLevels[level+1].fieldsRestricted[f];
    } );
    TOC()

    // Solve coarse grid problem
    if ( mgLevels[level+1].isCoarsestLevel ) {  // This is bad if there is 0 coarse levels
        TIC("Coarse grid smoothing")
        Smooth<MI, LI>( mgLevels[level+1], mgSettings.maxCoarseGridResiduals, mgSettings.maxCoarseGridIterations, MultigridEquation::TauCorrection );
        TOC()
    } else {
        VCycle<MI, LI>( mgLevels, startingLevel, level+1, mgSettings, iteration, axisTransformation );
    }

    // Compute fine grid correction 
    TIC("Fine grid correction computation")
    FieldData<Tensor3D> fineGridCorrection = ComputeFineGridCorrection( mgLevels[level+1].fields, 
                                                                        mgLevels[level+1].fieldsRestricted, 
                                                                        mgLevels[level+1].mesh, 
                                                                        mgLevels[level].mesh );
    TOC()

    // Correct fine grid approximation
    TIC("Fine grid correction application")
    ForAllFieldData( [&] (intType f) {
        mgLevels[level].fields[f] += fineGridCorrection[f];
    } );
    TOC()

    // Postsmoothing
    TIC("Post smoothing")
    if ( mgLevels[level].isFinestLevel ) {
        Smooth<MI, LI>( mgLevels[level], mgSettings.maxFineGridResiduals, mgSettings.maxFineGridIterations, mgEquationType );
    } else {
        Smooth<MI, LI>( mgLevels[level], mgSettings.maxPostSmoothingResiduals, mgSettings.maxPostSmoothingIterations, mgEquationType );
    }
    TOC()


}

// VCycle
template< MomentumInterpolation MI, Linearisation LI >
void MultigridCycle( std::vector< GridLevelData<MI, LI> > &mgLevels,
                     const InputData::MultigridSettings &mgSettings,
                     const intType iteration,
                     const AxisTransformationMap &axisTransformation )
{
    VCycle<MI, LI>( mgLevels, 0, 0, mgSettings, iteration, axisTransformation );
}


// // FMG
// template< MomentumInterpolation MI, Linearisation LI >
// void MultigridCycle( std::vector< GridLevelData<MI, LI> > &mgLevels,
//                      const InputData::MultigridSettings &mgSettings,
//                      const intType iteration,
//                      const AxisTransformationMap &axisTransformation )
// {
//     if ( iteration == 1 ) {
        
//         intType coarsestLevel = mgLevels.size() - 1;

//         // Solve coarsest level
//         Smooth<MI, LI>( mgLevels[coarsestLevel], mgSettings.maxCoarseGridResiduals, mgSettings.maxCoarseGridIterations, MultigridEquation::NoTauCorrection );

//         // Prolongate
//         ForAllFieldData( [&] (intType f) {
//             mgLevels[coarsestLevel-1].fields[f] = ProlongateField( mgLevels[coarsestLevel].fields[f], mgLevels[coarsestLevel].mesh, mgLevels[coarsestLevel-1].mesh );
//         } );

//         for ( intType level = coarsestLevel-1; level != 0; level-- ) {

//             // VCycle
//             VCycle<MI, LI>( mgLevels, level, level, mgSettings, iteration, axisTransformation );

//             // Prolongate
//             ForAllFieldData( [&] (intType f) {
//                 mgLevels[level-1].fields[f] = ProlongateField( mgLevels[level].fields[f], mgLevels[level].mesh, mgLevels[level-1].mesh );
//             } );

//         }

//     } else  {
     
//         VCycle<MI, LI>( mgLevels, 0, 0, mgSettings, iteration, axisTransformation );

//     }
// }


}   // end anonymous namespace



template< MomentumInterpolation MI, Linearisation LI >
void SweepSolve( const InputData &inputData,
                 const AxisTransformationMap &axisTransformation )
{
    using enum Axis::ENUMDATA;
    TIC("Pre processing")
    // Extract from input data
    const intType maxOuterIterations = inputData.schemes.maxOuterIterations;
    const FieldData<floatType> maxOuterResiduals = inputData.schemes.maxOuterResiduals;

    // Multigrid level data
    std::vector< GridLevelData<MI, LI> > mgLevels;  // FIND A BETTER WAY TO DO THIS
    SetMGLevels( mgLevels, inputData );

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

    // Outer iterations
    bool writeFields = ( inputData.fieldWriteInterval > 0 );
    if ( writeFields ) {
        fieldWriter.WriteData( 0 );
    }

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
    for ( intType nOuterIterations = 1; nOuterIterations <= maxOuterIterations; nOuterIterations++ )
    {
        if ( mgLevels.size() == 1 ) {
            Smooth<MI, LI>( mgLevels[0], 0, 1, MultigridEquation::NoTauCorrection );
        } else {
            MultigridCycle( mgLevels, inputData.multigridSettings, nOuterIterations, axisTransformation );
        }
        

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
        
        if ( ResidualsDiverged(residualsOuter) ) {
            fieldWriter.WriteData( nOuterIterations );
            std::cout << "*** SOLUTION DIVERGED ***" << "\n\n";
            break;
        }

        if ( nOuterIterations + 1 > maxOuterIterations ) {
            fieldWriter.WriteData( nOuterIterations );
            std::cout << "*** REACHED ITERATION LIMIT ***" << "\n\n";
            break;
        }

        if ( MetResidualTolerence(residualsOuter, maxOuterResiduals) ) {
            fieldWriter.WriteData( nOuterIterations );
            std::cout << "*** SOLUTION CONVERGED ***" << "\n\n";
            break;
        }

        if ( writeFields && (nOuterIterations % inputData.fieldWriteInterval) == 0 ) {
            fieldWriter.WriteData( nOuterIterations );
        }   
        
    }
    TOC()


}
template void SweepSolve<MomentumInterpolation::Implicit    , Linearisation::Picard>( const InputData &, const AxisTransformationMap &);
template void SweepSolve<MomentumInterpolation::SemiExplicit, Linearisation::Picard>( const InputData &, const AxisTransformationMap &);
template void SweepSolve<MomentumInterpolation::Implicit    , Linearisation::Newton>( const InputData &, const AxisTransformationMap &);
template void SweepSolve<MomentumInterpolation::SemiExplicit, Linearisation::Newton>( const InputData &, const AxisTransformationMap &);


} // end namespace CFD