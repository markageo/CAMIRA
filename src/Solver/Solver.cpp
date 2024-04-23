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
template< bool isNewtonLinearisation >
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
    if constexpr ( isNewtonLinearisation ) {
        UpdateFaceAdvectedVelocities(faceAdvectedVelocities, mesh, fields.U, faceFluxes, bcData);
    }
    SetIBFaceFluxes( faceFluxes, ibData );
    UpdateFVCoefficients(fvCoeffs, mesh, fields, faceAdvectedVelocities, faceFluxes, ibData, bcData);
}




template< MomentumInterpolation MI, Linearisation LI >
void VCycle( std::vector< GridLevelData<MI, LI> > &mgLevels,
             intType level )
{

    if ( mgLevels[level].isCoarsestGrid ) {

        // Solve coarsest grid

    } else {

        // Presmoothing 

        // Restrict residual

        // Restrict solution

        // VCycle recursive call
        VCycle();

        // Compute coarse grid approximation to the error

        // Interpolate error to the fine grid

        // Correct fine grid approximation

        // Postsmoothing

    }

}

template< MomentumInterpolation MI, Linearisation LI >
void MultigridCycle()
{

    // Perform a V Cycle
    VCycle<MI, LI>();

}


}   // end anonymous namespace



template< MomentumInterpolation MI, Linearisation LI >
void SweepSolve( FieldData<Tensor3D> &fields,               // TODO: Some of these arguments can be removed since they are now created with MG data. Maybe pass MG data into function?
                 const Mesh &mesh,
                 const BoundaryConditionData &bcData,
                 const InputData &inputData,
                 const AxisTransformationMap &axisTransformation)
{
    using enum Axis::ENUMDATA;
    
    // Extract from input data
    const intType maxOuterIterations = inputData.schemes.maxOuterIterations;
    const FieldData<floatType> maxOuterResiduals = inputData.schemes.maxOuterResiduals;

    // Multigrid level data
    std::vector< GridLevelData<MI, LI> > mgLevels = CreateMGLevels( inputData );

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
        MultigridCycle();

        residualsOuter   = StencilResiduals<MI, LI>(fields, fvCoeffs, ibData.mask); 
        NormaliseResiduals( residualsOuter, residualsScaleFactor, nOuterIterations );

        massFluxResidual = BoundaryMassFluxResidual(faceFluxes, mesh);

        probeValues      = SetFieldProbeValues(fields, fieldProbes); 
        
        fieldsOld = fields;

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
template void SweepSolve<MomentumInterpolation::Implicit    , Linearisation::Picard>( FieldData<Tensor3D> &, const Mesh &, const BoundaryConditionData &, const InputData &, const AxisTransformationMap &);
template void SweepSolve<MomentumInterpolation::SemiExplicit, Linearisation::Picard>( FieldData<Tensor3D> &, const Mesh &, const BoundaryConditionData &, const InputData &, const AxisTransformationMap &);
template void SweepSolve<MomentumInterpolation::Implicit    , Linearisation::Newton>( FieldData<Tensor3D> &, const Mesh &, const BoundaryConditionData &, const InputData &, const AxisTransformationMap &);
template void SweepSolve<MomentumInterpolation::SemiExplicit, Linearisation::Newton>( FieldData<Tensor3D> &, const Mesh &, const BoundaryConditionData &, const InputData &, const AxisTransformationMap &);


} // end namespace CFD