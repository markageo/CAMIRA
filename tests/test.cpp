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
    const CFD::Mesh mesh(inputData);

    // Write mesh to the output directory
    if (testConfig.meshTest != TEST::none) {
        TEST::WriteMesh(mesh, testConfig.meshTestOutputDirectory);
    } 
    
    // Compare with the mesh in the reference directory
    if (testConfig.meshTest == TEST::test) {
        if ( TEST::CompareMesh(testConfig.meshTestOutputDirectory, testConfig.meshTestReferenceDirectory) ){
            std::cout << "Meshes same" << "\n";
        } else {
            std::cout << "Meshes different!" << "\n";
        }
    } 
    

    /*-------------------------------------------------------------------------------------*\
                                          Fields Testing
    \*-------------------------------------------------------------------------------------*/

    using F = CFD::Fields::ENUMDATA;

    // Allocate fields and face velocities
    CFD::ArrayAllocator<CFD::Fields, CFD::array3D> fields({F::U, F::V, F::W, F::P}, mesh.nCells);


    // Faces are staggered in the negative direction:
    //   cellFaceVelocity_x(i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceVelocity_y(i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceVelocity_z(i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    CFD::ArrayAllocator<CFD::Fields, CFD::array3D> faceVelocities( {{F::U, {mesh.nCells(0) + 1, mesh.nCells(1)    , mesh.nCells(2)    }},
                                                                    {F::V, {mesh.nCells(0)    , mesh.nCells(1) + 1, mesh.nCells(2)    }},
                                                                    {F::W, {mesh.nCells(0)    , mesh.nCells(1)    , mesh.nCells(2) + 1}}} );


    // Set the seed so the random arrays are reproducable
    srand(590);

    // Set values for the velocity fields
    fields[F::U].setRandom();
    fields[F::V].setRandom();
    fields[F::W].setRandom();

    // Update the face velocities
    CFD::UpdateFaceVelocities(faceVelocities, mesh, fields, inputData.boundaryConditions);

    /*-------------------------------------------------------------------------------------*\
                                Finite Volume Coeffiicents Testing
    \*-------------------------------------------------------------------------------------*/
    
    // Create initial finite volume coefficients
    CFD::FVCoefficients fvCoeffs = CFD::InitialiseFVCoefficients(mesh, faceVelocities, inputData);


    return 0;
}
