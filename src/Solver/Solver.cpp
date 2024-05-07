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

// Update face fluxes, immersed boundary data, and finite volume equation coefficients
template< Linearisation LI >
void UpdateFVEquations( FVCoefficients &fvCoeffs,
                        IBData &ibData,
                        EnumVector<Axis, Tensor3D> &faceFluxes,
                        EnumVector< Axis, EnumVector< Axis, Tensor3D> > &faceAdvectedVelocities,
                        const FieldData<Tensor3D> &fields,
                        const Mesh &mesh,
                        const BoundaryConditionData &bcData )
{
    UpdateIBData( ibData, fields );
    UpdateFaceFluxes(faceFluxes, mesh, fields.U, bcData);
    if constexpr ( LI == Linearisation::Newton ) {
        UpdateFaceAdvectedVelocities(faceAdvectedVelocities, mesh, fields.U, faceFluxes, bcData);
    }
    SetIBFaceFluxes( faceFluxes, ibData );
    UpdateFVCoefficients(fvCoeffs, mesh, fields, faceAdvectedVelocities, faceFluxes, ibData, bcData);
}



template< MomentumInterpolation MI, Linearisation LI >
void Smooth( GridLevelData<MI, LI> &gridLevelData,
             const FieldData<floatType> &maxResiduals,
             const intType maxIterations )
{
    FieldData<floatType> residualsScaleFactor;

    for ( intType nIterations = 1; nIterations <= maxIterations; nIterations++ ) {

        gridLevelData.linearSolver->UpdateState();
        gridLevelData.linearSolver->Solve();
        // SetGhostCells( gridLevelData.fields, gridLevelData.mesh, gridLevelData.bcData ); // Deferred correction needs this
        UpdateFVEquations<LI>( gridLevelData.fvCoeffs, 
                               gridLevelData.ibData, 
                               gridLevelData.faceFluxes, 
                               gridLevelData.faceAdvectedVelocities, 
                               gridLevelData.fields, 
                               gridLevelData.mesh, 
                               gridLevelData.bcData );

        if ( !gridLevelData.isFinestLevel ) {
            TransformToCoarseGridEquations<MI, LI>( gridLevelData.fvCoeffs, 
                                                    gridLevelData.fieldsRestricted,
                                                    gridLevelData.residuals,
                                                    gridLevelData.ibData.mask );
        }

        FieldData<floatType> residuals = ScaledL1NormResiduals<MI, LI>( gridLevelData.fields,
                                                                        gridLevelData.fvCoeffs,
                                                                        gridLevelData.ibData.mask );
        NormaliseResiduals( residuals, residualsScaleFactor, nIterations );

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

}



template< MomentumInterpolation MI, Linearisation LI >
void VCycle( std::vector< GridLevelData<MI, LI> > &mgLevels,
             const intType level,
             const InputData::MultigridSettings &mgSettings )
{

    if ( mgLevels[level].isCoarsestLevel ) {

        // Solve coarsest grid
        Smooth<MI, LI>( mgLevels[level], mgSettings.maxCoarseGridResiduals, mgSettings.maxCoarseGridIterations );

    } else {

        // Presmoothing 
        Smooth<MI, LI>( mgLevels[level], mgSettings.maxPreSmoothingResiduals, mgSettings.maxPreSmoothingIterations );

        // Calculate residual on finest grid
        if ( mgLevels[level].isFinestLevel ) {
            mgLevels[level+1].residuals= ResidualsField<MI, LI>( mgLevels[level].fields, 
                                                                 mgLevels[level].fvCoeffs, 
                                                                 mgLevels[level].ibData.mask );
        }
        
        // Restrict residual
        ForAllFieldData( [&] (intType f) {
            mgLevels[level+1].residuals[f] = RestrictField( mgLevels[level].residuals[f], 
                                                            mgLevels[level].mesh, 
                                                            mgLevels[level+1].mesh );
        } );

        // Restrict solution
        ForAllFieldData( [&] (intType f) {
            mgLevels[level+1].fieldsRestricted[f] = RestrictField( mgLevels[level].fieldsRestricted[f], 
                                                                   mgLevels[level].mesh, 
                                                                   mgLevels[level+1].mesh );
        } );

        // VCycle recursive call
        VCycle<MI, LI>(mgLevels, level+1, mgSettings);

        // Compute fine grid correction 
        FieldData<Tensor3D> fineGridCorrection = ComputeFineGridCorrection( mgLevels[level+1].fields, 
                                                                            mgLevels[level+1].fieldsRestricted, 
                                                                            mgLevels[level+1].mesh, 
                                                                            mgLevels[level].mesh );

        // Correct fine grid approximation
        ForAllFieldData( [&] (intType f) {
            mgLevels[level].fields[f] += fineGridCorrection[f];
        } );

        // Postsmoothing
        Smooth<MI, LI>( mgLevels[level], mgSettings.maxPostSmoothingResiduals, mgSettings.maxPostSmoothingIterations );

    }

}

template< MomentumInterpolation MI, Linearisation LI >
void MultigridCycle( std::vector< GridLevelData<MI, LI> > &mgLevels,
                     const InputData::MultigridSettings &mgSettings )
{
    VCycle<MI, LI>( mgLevels, 0, mgSettings );
}


}   // end anonymous namespace



template< MomentumInterpolation MI, Linearisation LI >
void SweepSolve( const InputData &inputData,
                 const AxisTransformationMap &axisTransformation )
{
    using enum Axis::ENUMDATA;
    
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

    TIC("Solver Loop")
    for ( intType nOuterIterations = 1; nOuterIterations <= maxOuterIterations; nOuterIterations++ )
    {
        MultigridCycle( mgLevels, inputData.multigridSettings );

        residualsOuter   = ScaledL1NormResiduals<MI, LI>( mgLevels[0].fields, 
                                                     mgLevels[0].fvCoeffs, 
                                                     mgLevels[0].ibData.mask); 

        NormaliseResiduals( residualsOuter, residualsScaleFactor, nOuterIterations );

        massFluxResidual = BoundaryMassFluxResidual( mgLevels[0].faceFluxes, 
                                                     mgLevels[0].mesh);

        probeValues      = SetFieldProbeValues( mgLevels[0].fields, fieldProbes); 
        
        for ( auto &mgLevel : mgLevels ) {
            mgLevel.fieldsOld = mgLevel.fields;
        }

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