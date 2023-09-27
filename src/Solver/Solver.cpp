#include "Solver.h"
#include "ConvergenceLogging.h"
#include "FieldProbe.h"

#include "../Types.h"
#include "../Macros.h"
#include "../IO/InputProcessing.h"
#include "../Tools/SweepTransformations.h"
#include "../FiniteVolume/FiniteVolume.h"

#include "TriadSolver.h"
#include "LineSolver.h"
#include "PlaneSolver.h"
#include "LinearSolver.h"
#include "ResidualFunctions.h"

#include <iostream>

namespace CFD
{

template< MomentumInterpolation MI, Linearisation LI >
void SweepSolve( FieldData<array3D> &fields,
                 const Mesh &mesh,
                 const FieldData< BoundaryConditionData > &bcData,
                 const InputData &inputData,
                 const AxisTransformationMap &axisTransformation)
{
    using enum Axis::ENUMDATA;

    constexpr bool isNewtonLinearisation = ( LI == Linearisation::Newton );

    // Extract from input data
    const InputData::LinearSolverSettings linearSolverSettings = inputData.linearSolverSettings;
    const intType maxOuterIterations = inputData.schemes.maxOuterIterations;
    const FieldData<floatType> maxOuterResiduals = inputData.schemes.maxOuterResiduals;

    // Initialise
    EnumVector<Axis, array3D> faceFluxes = InitialiseFaceFluxes(mesh, fields.U, bcData);
    EnumVector< Axis, EnumVector< Axis, array3D> > faceAdvectedVelocities;
    if constexpr ( isNewtonLinearisation ) 
        faceAdvectedVelocities = InitialiseAdvectedFaceVelocities( mesh, fields.U, faceFluxes, bcData );

    FieldData<array3D> fieldsOld = fields;
    FVCoefficients fvCoeffs = InitialiseFVCoefficients(mesh, fields, faceAdvectedVelocities, faceFluxes, bcData, inputData);

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

    // Instantiate linear solver, this holds references to the fields
    LinearSolver<MI, LI> linearSolver(fields, fieldsOld, fvCoeffs, linearSolverSettings);


    // Outer iterations
    bool writeFields = ( inputData.fieldWriteInterval > 0 );
    if ( writeFields ) {
        fieldWriter.WriteData( 0 );
    }
        
    TIC("Solver Loop")
    for ( intType nOuterIterations = 1; nOuterIterations <= maxOuterIterations; nOuterIterations++ )
    {
        TIC("Solver Update State")
        linearSolver.UpdateState();
        TOC()

        TIC("Linear Solver")
        linearSolver.Solve();
        TOC()

        TIC("Update face values")
        UpdateFaceFluxes(faceFluxes, mesh, fields.U, bcData);
        if constexpr ( isNewtonLinearisation ) {
            UpdateFaceAdvectedVelocities(faceAdvectedVelocities, mesh, fields.U, faceFluxes, bcData);
        }
        TOC()

        TIC("Update Coefficients")
        UpdateFVCoefficients(fvCoeffs, mesh, fields, faceAdvectedVelocities, faceFluxes, bcData);
        TOC()

        TIC("Residuals and logging")
        residualsOuter   = StencilResiduals<MI, LI>(fields, fvCoeffs); 
        NormaliseResiduals( residualsOuter, residualsScaleFactor, nOuterIterations );

        massFluxResidual = BoundaryMassFluxResidual(faceFluxes, mesh);

        probeValues      = SetFieldProbeValues(fields, fieldProbes); 
        
        fieldsOld = fields;

        consoleLog.WriteResiduals( residualsOuter, massFluxResidual, nOuterIterations );
        residualsLogFile.WriteData( residualsOuter, massFluxResidual, nOuterIterations );
        for ( size_t p = 0; p != fieldProbes.size(); p++ ) {
            probeLogFiles[p].WriteData( probeValues[p], nOuterIterations );
        }

        TOC()
        
        if ( ResidualsDiverged(residualsOuter) ) {
            fieldWriter.WriteData( nOuterIterations );
            std::cout << "*** SOLUTION DIVERGED ***" << "\n\n";
            break;
        }

        if ( MetResidualTolerence(residualsOuter, maxOuterResiduals) ) {
            fieldWriter.WriteData( nOuterIterations );
            std::cout << "*** SOLUTION CONVERGED ***" << "\n\n";
            break;
        }

        TIC("Writing Fields")
        if ( writeFields && (nOuterIterations % inputData.fieldWriteInterval) == 0 ) {
            fieldWriter.WriteData( nOuterIterations );
        }
        TOC()
        

    }
    TOC()


}
template void SweepSolve<MomentumInterpolation::Implicit    , Linearisation::Picard>( FieldData<array3D> &, const Mesh &, const FieldData< BoundaryConditionData > &, const InputData &, const AxisTransformationMap &);
template void SweepSolve<MomentumInterpolation::SemiExplicit, Linearisation::Picard>( FieldData<array3D> &, const Mesh &, const FieldData< BoundaryConditionData > &, const InputData &, const AxisTransformationMap &);
template void SweepSolve<MomentumInterpolation::Implicit    , Linearisation::Newton>( FieldData<array3D> &, const Mesh &, const FieldData< BoundaryConditionData > &, const InputData &, const AxisTransformationMap &);
template void SweepSolve<MomentumInterpolation::SemiExplicit, Linearisation::Newton>( FieldData<array3D> &, const Mesh &, const FieldData< BoundaryConditionData > &, const InputData &, const AxisTransformationMap &);


} // end namespace CFD