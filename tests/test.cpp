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
#include "Utils.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"
#include "Solver.h"
#include "VTKWriter.h"
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
    CFD::Mesh mesh(inputData);

    // Write mesh to the output directory
    if (testConfig.meshTest != TEST::none) {
        TEST::WriteMesh(mesh, testConfig.meshTestOutputDirectory);
        std::cout << "Mesh written to output directory." << "\n";
    } 
    
    // Compare with the mesh in the reference directory
    if (testConfig.meshTest == TEST::test) {
        if ( TEST::CompareMesh(testConfig.meshTestOutputDirectory, testConfig.meshTestReferenceDirectory) ){
            std::cout << "Mesh test: PASSED!" << "\n";
        } else {
            std::cout << "Mesh test: FAILED! Mesh does not match reference!" << "\n";
        }
    } 
    std::cout << "\n";
    

    /*-------------------------------------------------------------------------------------*\
                         Face velocities and Finite Volume Coefficients
    \*-------------------------------------------------------------------------------------*/

    using F = CFD::Fields::ENUMDATA;
    using BC = CFD::BoundaryConditions::ENUMDATA;

    // Read and store boundary conditions for each case
    CFD::EnumVector<CFD::BoundaryConditions, CFD::InputData> velTestInputData;
    velTestInputData[BC::uniform] = CFD::ReadInputData(testConfig.testInputDirectory + testConfig.velUniformTestInputFilename);
    velTestInputData[BC::zeroGradient] = CFD::ReadInputData(testConfig.testInputDirectory + testConfig.velZeroGradientTestInputFilename);
    velTestInputData[BC::extrapolated] = CFD::ReadInputData(testConfig.testInputDirectory + testConfig.velExtrapolatedTestInputFilename);

    // Generate new mesh
    mesh = CFD::Mesh(velTestInputData[BC::uniform]);

    // Allocate fields
    CFD::ArrayAllocator<CFD::Fields, CFD::array3D> fields({F::U, F::V, F::W, F::P}, mesh.nCells + 2*CFD::nGhost);

    // Allocate face velocities for each bouundary condition type to be tested
    // Faces are staggered in the negative direction:
    //   cellFaceVelocity_x(i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceVelocity_y(i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceVelocity_z(i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    CFD::EnumVector< CFD::BoundaryConditions, CFD::ArrayAllocator<CFD::Fields, CFD::array3D> > 
        testFaceVelocities( CFD::ArrayAllocator<CFD::Fields, CFD::array3D>({{F::U, {mesh.nCells(0) + 1, mesh.nCells(1)    , mesh.nCells(2)    }},
                                                                            {F::V, {mesh.nCells(0)    , mesh.nCells(1) + 1, mesh.nCells(2)    }},             
                                                                            {F::W, {mesh.nCells(0)    , mesh.nCells(1)    , mesh.nCells(2) + 1}}}) );

    CFD::EnumVector< CFD::BoundaryConditions, CFD::FVCoefficients > testFVCoeffs( CFD::FVCoefficients(mesh.nCells) );

    // Set values for the fields to be random and write them out to the test directory
    // Make all the fields the same, they should be rotated to give the same results
    srand(590);
    Eigen::array< std::pair<int, int>, 3 > paddings;
    paddings[0] = std::make_pair(CFD::nGhost, CFD::nGhost);
    paddings[1] = std::make_pair(CFD::nGhost, CFD::nGhost);
    paddings[2] = std::make_pair(CFD::nGhost, CFD::nGhost);
    CFD::array3D randTemp(mesh.nCells(0), mesh.nCells(1), mesh.nCells(2));
    randTemp = randTemp.setRandom() * randTemp.constant( 2 ) + randTemp.constant( -1 ); 
    fields[F::U] = randTemp.pad(paddings);
    fields[F::V] = fields[F::U].shuffle( Eigen::array<int, 3>{1, 0, 2} );
    fields[F::W] = fields[F::U].shuffle( Eigen::array<int, 3>{2, 1, 0} );
    fields[F::P] = fields[F::U];
    TEST::WriteFields(fields, testConfig.faceVelTestOutputDirectory);


    // Update the face velocities and set finite volume coefficients
    CFD::BoundaryConditions::ENUMDATA boundaryCondition;
    for (int i = 0; i != CFD::BoundaryConditions::count; i++) {
        boundaryCondition = static_cast<CFD::BoundaryConditions::ENUMDATA>(i);
        CFD::UpdateFaceVelocities(testFaceVelocities[boundaryCondition], mesh, fields, velTestInputData[boundaryCondition]);
        testFVCoeffs[boundaryCondition] = CFD::InitialiseFVCoefficients( mesh, fields, testFaceVelocities[boundaryCondition], velTestInputData[boundaryCondition] );
    }

    // Face velocities
    if (testConfig.faceVelTest != TEST::none) {
        TEST::WriteFaceVels(testFaceVelocities, testConfig.faceVelTestOutputDirectory);
        TEST::WriteMesh(mesh, testConfig.faceVelTestOutputDirectory);
        std::cout << "Face velocities and associated mesh written to output directory" << "\n";
    } 

    if (testConfig.faceVelTest == TEST::test) {
        if ( TEST::CompareFaceVels(testConfig.faceVelTestOutputDirectory, testConfig.faceVelTestReferenceDirectory) ){
            std::cout << "Face velocity test: PASSED!" << "\n";
        } else {
            std::cout << "Face velocity test: FAILED! Face velocities do not match reference!" << "\n";
        }
    } 
    std::cout << "\n";


    // Finite volume coefficients
    if (testConfig.fvCoeffTest != TEST::none) {
        TEST::WriteFVCoeffs(testFVCoeffs, testConfig.fvCoeffTestOutputDirectory);
        std::cout << "Finite volume coefficients written to output directory" << "\n";
    } 

    if (testConfig.fvCoeffTest == TEST::test) {
        if ( TEST::CompareFVCoeffs(testConfig.fvCoeffTestOutputDirectory, testConfig.fvCoeffTestReferenceDirectory) ){
            std::cout << "Finite volume coefficients test: PASSED!" << "\n";
        } else {
            std::cout << "Finite volume coefficients test: FAILED! Face velocities do not match reference!" << "\n";
        }
    } 
    std::cout << "\n";





    return 0;
}
