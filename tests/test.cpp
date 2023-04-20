/*---------------------------------------------------------------------------*\
    3D coupled Navier-Stokes solver on rectilinear grids
    TESTS

    Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"
#include "utils.h"
#include "InputProcessing.h"
#include "FiniteVolumeStructures.h"
#include "FiniteVolumeFunctions.h"
#include "VTKWriter.h"
#include "Solver.h"

#include "TestFunctions.h"

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

    // Read input file from command line
    CFD::InputData inputData = CFD::InputDataFromCommandLine(argc, argv);

    // Useful enums
    using AX = CFD::Axis::ENUMDATA;
    using F = CFD::Fields::ENUMDATA;

    /*-------------------------------------------------------------------------------------*\
                                            Mesh Testing
    \*-------------------------------------------------------------------------------------*/

    // Generate mesh
    TIC("Meshing");
    const CFD::Mesh mesh(inputData);
    TOC();

    // Write mesh data to file for each axis
    TEST::WriteMesh(mesh, "tests/mesh/output/");


    /*-------------------------------------------------------------------------------------*\
                                      Face Velocities Testing
    \*-------------------------------------------------------------------------------------*/

    // Allocate fields and face velocities
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


    // Set values for the velocity fields
    TIC("Set Values");
    fields[F::U].setRandom();
    fields[F::V].setRandom();
    fields[F::W].setRandom();
    TOC();

    

    // Update the face velocities
    TIC("Face Velocity Update")
    CFD::UpdateFaceVelocities(faceVelocities, mesh, fields, inputData.boundaryConditions);
    TOC();

    UTIL::WriteArray("tests/test_velocities.dat", fields[F::U]);
    CFD::array3D array = UTIL::ReadArray<CFD::array3D>("tests/test_velocities.dat");
    UTIL::WriteArray("tests/test_velocities_rewritten.dat", array);


    /*-------------------------------------------------------------------------------------*\
                                Finite Volume Coeffiicents Testing
    \*-------------------------------------------------------------------------------------*/

    
    // Create initial finite volume coefficients
    TIC("Allocate and Initialise FV Coefficients")
    CFD::FVCoefficients fvCoeffs = CFD::InitialiseFVCoefficients(mesh, faceVelocities, inputData);
    TOC()

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
    VTK::scalarMapType scalarMap = {{"U", fields[F::U].data()}};
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
