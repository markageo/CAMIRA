/*---------------------------------------------------------------------------*\
   3D coupled Navier-Stokes solver on rectilinear grids

   Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"
#include "InputProcessing.h"
#include "SweepTransformations.h"
#include "FiniteVolume.h"
#include "VTKWriter.h"
#include "Solver.h"
#include "Utils.h"
#include <iostream>
#include <fmt/core.h>

int main(int argc, char const *argv[])
{

    /*-------------------------------------------------------------------------------------*\
                                         Input Processing
    \*-------------------------------------------------------------------------------------*/

    CFD::InputData inputData = CFD::InputDataFromCommandLine(argc, argv);

    CFD::AxisTransformationMap axisTransformation = CFD::TransformUserInputData( inputData );


    /*-------------------------------------------------------------------------------------*\
                                              Solve
    \*-------------------------------------------------------------------------------------*/

    TIC("Meshing");
    CFD::Mesh mesh( inputData );
    TOC();

    TIC("Field Allocation");
    CFD::FieldData<CFD::array3D> fields = CFD::InitialiseFields( mesh, inputData );
    TOC();

    TIC("Solver");
    CFD::SweepSolve(fields, mesh, inputData, axisTransformation);
    TOC();


    /*-------------------------------------------------------------------------------------*\
                                         Post-Processing
    \*-------------------------------------------------------------------------------------*/

    // Remove ghost cells from the fields
    CFD::ForAllFieldData( [&] (CFD::intType f) {
        CFD::RemoveGhostCells(fields[f], CFD::nGhost);
    } );

    CFD::FieldData<CFD::array3D> vertexFields = GetVertexFields( fields, mesh, inputData );

    // Undo the boundary condition transformation
    CFD::TransformToUserCoordinates(mesh, fields, vertexFields, axisTransformation);


    /*-------------------------------------------------------------------------------------*\
                                             Output
    \*-------------------------------------------------------------------------------------*/

    using enum CFD::Axis::ENUMDATA;

    // Data to pass to writer
    // VTK::VTKWriterConfig config( mesh.nCells[X], mesh.nCells[Y], mesh.nCells[Z] );
    //     config.SetWriteMode( VTK::WriteModes::ASCII );
    //     config.SetGridType( VTK::GridTypes::CELL_DATA );
    // VTK::gridVectorType<CFD::floatType> gridVector = {mesh.cellFaces[X].data(), mesh.cellFaces[Y].data(), mesh.cellFaces[Z].data()};
    // VTK::scalarMapType<CFD::floatType> scalarMap = { {"Pressure", fields.P.data()} };
    // VTK::vectorMapType<CFD::floatType> vectorMap = { {"Velocity", {fields.U[X].data(), fields.U[Y].data(), fields.U[Z].data()}} };


     VTK::VTKWriterConfig config( mesh.nCells[X]+1, mesh.nCells[Y]+1, mesh.nCells[Z]+1 );
        config.SetWriteMode( VTK::WriteModes::ASCII );
        config.SetGridType( VTK::GridTypes::POINT_DATA );
    VTK::gridVectorType<CFD::floatType> gridVector = {mesh.cellFaces[X].data(), mesh.cellFaces[Y].data(), mesh.cellFaces[Z].data()};
    VTK::scalarMapType<CFD::floatType> scalarMap = { {"Pressure", vertexFields.P.data()} };
    VTK::vectorMapType<CFD::floatType> vectorMap = { {"Velocity", {vertexFields.U[X].data(), vertexFields.U[Y].data(), vertexFields.U[Z].data()}} };

    // Write output
    VTK::VTKWriter writer(gridVector, scalarMap, vectorMap, config);
    TIC("Writer");
    writer.WriteData(inputData.fieldOutputFilename, "CFD simulation output");
    TOC();

// Display profiling information
#ifdef PROFILING
    std::cout << PROF::prof;
#endif

    return 0;
}