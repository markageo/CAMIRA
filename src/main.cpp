 /*---------------------------------------------------------------------------*\
    3D coupled Navier-Stokes solver on rectilinear grids

    Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"
#include "Utils.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"
#include "VTKWriter.h"
#include "Solver.h"

#include <iostream>
#include <fmt/core.h>

#ifdef PROFILING
#include "profiler/profiler.h"
    namespace PROF {
        profiler<perf_counter::clock<time_units::SECONDS>> prof;
    }
#endif

int main(int argc, char const *argv[])
{

    /*-------------------------------------------------------------------------------------*\
                                         Input Processing
    \*-------------------------------------------------------------------------------------*/

    CFD::InputData inputData = CFD::InputDataFromCommandLine(argc, argv);

    /*-------------------------------------------------------------------------------------*\
                                           Solve
    \*-------------------------------------------------------------------------------------*/

    using F = CFD::Fields::ENUMDATA;

    TIC("Meshing");
    CFD::Mesh mesh(inputData);
    TOC();

    TIC("Field Allocation");
    CFD::ArrayAllocator<CFD::Fields, CFD::array3D> fields({F::U, F::V, F::W, F::P}, mesh.nCells + 2*CFD::nGhost);
    TOC();

    TIC("Solver");
    CFD::SweepSolve(fields, mesh, inputData);
    TOC();


    /*-------------------------------------------------------------------------------------*\
                                           Output
    \*-------------------------------------------------------------------------------------*/

    using AX = CFD::Axis::ENUMDATA;

    // Undo the boundary condition transformation
    CFD::TransformToUserCoordinates(mesh, fields, inputData.axisTransformation);

    // Remove ghost cells from the fields
    CFD::RemoveGhostCells(fields, CFD::nGhost);


    // Data to pass to writer
    VTK::VTKWriterConfig config(mesh.nCells[AX::X], mesh.nCells[AX::Y], mesh.nCells[AX::Z]);
        config.SetWriteMode(VTK::ASCII);
    VTK::gridVectorType<CFD::floatType> gridVector = {mesh.cellCenters[AX::X].data(), mesh.cellCenters[AX::Y].data(), mesh.cellCenters[AX::Z].data()};
    VTK::scalarMapType<CFD::floatType> scalarMap = { {"Pressure", fields[F::P].data()} };
    VTK::vectorMapType<CFD::floatType> vectorMap = { {"Velocity", { fields[F::U].data(), fields[F::V].data(), fields[F::W].data() } } };

    // Write output
    VTK::VTKWriter writer(gridVector, scalarMap, vectorMap, config);
    TIC("Writer");
    writer.WriteData("fields.vtk", "CFD simulation");
    TOC();

    // Display profiling information
    #ifdef PROFILING
        std::cout << PROF::prof;
    #endif

    return 0;
}
