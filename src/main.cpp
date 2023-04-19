/*---------------------------------------------------------------------------*\
    3D coupled Navier-Stokes solver on rectilinear grids

    Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"
#include "utils.h"
#include "InputProcessing.h"
#include "FiniteVolumeStructures.h"
#include "FiniteVolumeFunctions.h"
#include "VTKWriter.h"
#include "Solver.h"

#include <iostream>
#include <fmt/core.h>

#ifdef PROFILING
#include "profiler/profiler.h"
namespace PROF
{
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

    using AX = CFD::Axis::ENUMDATA;
    using F = CFD::Fields::ENUMDATA;

    TIC("Meshing");
    const CFD::Mesh mesh(inputData);
    TOC();

    TIC("Field Allocation");
    CFD::ArrayAllocator<CFD::Fields, CFD::array3D> fields({F::U, F::V, F::W, F::P}, mesh.nCells);

    // Faces are staggered in the negative direction:
    //   cellFaceVelocity_x(i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceVelocity_y(i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceVelocity_z(i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    CFD::ArrayAllocator<CFD::Fields, CFD::array3D> faceVelocities( {{F::U, {mesh.nCells(0) + 1, mesh.nCells(1)    , mesh.nCells(2)    }},
                                                                    {F::V, {mesh.nCells(0)    , mesh.nCells(1) + 1, mesh.nCells(2)    }},
                                                                    {F::W, {mesh.nCells(0)    , mesh.nCells(1)    , mesh.nCells(2) + 1}}} );
    TOC();



    /*-------------------------------------------------------------------------------------*\
                                           Output
    \*-------------------------------------------------------------------------------------*/

    // Data to pass to writer
    VTK::dataType VTKDataType = VTK::DOUBLE;
    if (std::is_same<CFD::floatType, float>::value)
    {
        VTKDataType = VTK::FLOAT;
    }
    VTK::VTKWriterConfig config(mesh.nCells[AX::X], mesh.nCells[AX::Y], mesh.nCells[AX::Z], VTKDataType);
        config.SetWriteMode("binary");
    VTK::gridVectorType gridVector = {mesh.cellCenters[AX::X].data(), mesh.cellCenters[AX::Y].data(), mesh.cellCenters[AX::Z].data()};
    VTK::scalarMapType scalarMap = {};
    VTK::vectorMapType vectorMap = {};

    // Write output
    VTK::VTKWriter writer(gridVector, scalarMap, vectorMap, config);
    TIC("Writer");
    writer.WriteData("mesh.vtk", "3D Rectilinear Grid");
    TOC();

// Display profiling information
#ifdef PROFILING
    std::cout << PROF::prof;
#endif

    return 0;
}
