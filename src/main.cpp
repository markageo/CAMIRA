/*---------------------------------------------------------------------------*\
   3D coupled Navier-Stokes solver on rectilinear grids

   Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"
#include "Macros.h"

#include "IO/InputProcessing.h"
#include "IO/VTKWriter.h"
#include "Tools/SweepTransformations.h"
#include "Tools/FVTools.h"
#include "FiniteVolume/FiniteVolume.h"
#include "Solver/Solver.h"

#include <iostream>
#include <fmt/core.h>

int main(int argc, char const *argv[])
{

    /*-------------------------------------------------------------------------------------*\
                                         Input Processing
    \*-------------------------------------------------------------------------------------*/

    TIC("Input Processing")
    CFD::InputData inputData = CFD::InputDataFromCommandLine(argc, argv);

    CFD::AxisTransformationMap axisTransformation = CFD::TransformUserInputData(inputData);
    TOC()

    /*-------------------------------------------------------------------------------------*\
                                              Solve
    \*-------------------------------------------------------------------------------------*/

    TIC("Meshing");
    CFD::Mesh mesh(inputData);
    TOC();

    TIC("Field Allocation");
    CFD::FieldData<CFD::array3D> fields = CFD::InitialiseFields(mesh, inputData);
    TOC();

    TIC("Boundary Condition Processing")
    CFD::FieldData<CFD::BoundaryConditionData> bcData = SetBoundaryConditionData(inputData, mesh);
    TOC()

    TIC("Solver");
    switch (inputData.schemes.momentumInterpolation)
    {
        using MI = CFD::MomentumInterpolation;
    case (MI::Implicit):
        CFD::SweepSolve<MI::Implicit>(fields, mesh, bcData, inputData, axisTransformation);
        break;

    case (MI::SemiExplicit):
        CFD::SweepSolve<MI::SemiExplicit>(fields, mesh, bcData, inputData, axisTransformation);
        break;
    }
    TOC();

    /*-------------------------------------------------------------------------------------*\
                                         Post-Processing
    \*-------------------------------------------------------------------------------------*/

    TIC("Post Processing");
    // Remove ghost cells from the fields
    CFD::ForAllFieldData([&](CFD::intType f)
                         { CFD::FVT::RemoveGhostCells(fields[f], CFD::nGhost); });

    CFD::FieldData<CFD::array3D> vertexFields = GetVertexFields(fields, mesh, bcData);

    // Undo the boundary condition transformation
    CFD::TransformToUserCoordinates(mesh, fields, vertexFields, axisTransformation);
    TOC();

    /*-------------------------------------------------------------------------------------*\
                                             Output
    \*-------------------------------------------------------------------------------------*/

    using enum CFD::Axis::ENUMDATA;

    VTK::VTKWriterConfig config(mesh.nFacesNormal[X](X), mesh.nFacesNormal[Y](Y), mesh.nFacesNormal[Z](Z));
    config.SetWriteMode(VTK::WriteModes::BINARY);
    VTK::gridVectorType<CFD::floatType> gridVector = {mesh.cellFaces[X].data(), mesh.cellFaces[Y].data(), mesh.cellFaces[Z].data()};

    VTK::scalarMapType<CFD::floatType> scalarMap = {{"Pressure", VTK::GridTypes::CELL_DATA, fields.P.data()},
                                                    {"Pressure", VTK::GridTypes::POINT_DATA, vertexFields.P.data()}};

    VTK::vectorMapType<CFD::floatType> vectorMap = {{"Velocity", VTK::GridTypes::POINT_DATA, {vertexFields.U[X].data(), vertexFields.U[Y].data(), vertexFields.U[Z].data()}},
                                                    {"Velocity", VTK::GridTypes::CELL_DATA, {fields.U[X].data(), fields.U[Y].data(), fields.U[Z].data()}}};
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