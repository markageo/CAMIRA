/*---------------------------------------------------------------------------*\
    3D coupled Navier-Stokes solver on rectilinear grids
    TESTS

    Mark George
\*---------------------------------------------------------------------------*/

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <fmt/core.h>

#include "Types.h"
#include "utils.h"
#include "InputProcessing.h"
#include "FiniteVolumeStructures.h"
#include "FiniteVolumeFunctions.h"
#include "VTKWriter.h"
#include "Solver.h"
#include "TestFunctions.h"



#ifdef PROFILING
#include "profiler/profiler.h"
    namespace PROF {
        profiler<perf_counter::clock<time_units::SECONDS>> prof;
    }
#endif

int main(int argc, char const *argv[])
{

    /*-------------------------------------------------------------------------------------*\
                                        Config/Input Files
    \*-------------------------------------------------------------------------------------*/

    // Conifg file
    std::string configFilename;
    if (argc == 1) {
        std::cout << "Please enter config filename: " << "\n";
        std::cin >> configFilename;
        std::cout << "\n";
    }
    else if (argc == 2) {
        configFilename = argv[1];
    }
    else {
        throw std::invalid_argument("Invalid command line options.");
    }
    TEST::TestConfig testConfig = TEST::ReadConfig(configFilename);


    /*-------------------------------------------------------------------------------------*\
                                            Mesh Testing
    \*-------------------------------------------------------------------------------------*/

    // Read input file for mesh testing
    CFD::InputData inputData = CFD::ReadInputData(testConfig.testInputDirectory + testConfig.meshTestInputFilename);

    // Generate mesh
    TIC("Meshing");
    const CFD::Mesh mesh(inputData);
    TOC();

    if (testConfig.meshTest == TEST::write) {

        // Write mesh data to file for each axis
        TEST::WriteMesh(mesh, testConfig.meshTestOutputDirectory);

    } else if (testConfig.meshTest == TEST::test) {

        // Compare with the generated mesh data
    

        // Display the result

    } 
    

    /*-------------------------------------------------------------------------------------*\
                                          Fields Testing
    \*-------------------------------------------------------------------------------------*/

    using F = CFD::Fields::ENUMDATA;

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


    // Set the seed so the random arrays are reproducable
    srand(590);

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

    // Write face velocity
    UTIL::WriteArray("tests/velocities.dat", fields[F::U]);

    /*-------------------------------------------------------------------------------------*\
                                Finite Volume Coeffiicents Testing
    \*-------------------------------------------------------------------------------------*/

    
    // Create initial finite volume coefficients
    TIC("Allocate and Initialise FV Coefficients")
    CFD::FVCoefficients fvCoeffs = CFD::InitialiseFVCoefficients(mesh, faceVelocities, inputData);
    TOC()





    // Display profiling information
    #ifdef PROFILING
        std::cout << PROF::prof;
    #endif

    return 0;
}
