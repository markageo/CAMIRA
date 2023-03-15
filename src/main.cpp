/*---------------------------------------------------------------------------*\
    3D coupled Navier-Stokes solver on rectilinear grids

    Mark George
\*---------------------------------------------------------------------------*/

#include "SimulationParameters.h"
#include "InputProcessing.h"
#include "FiniteVolumeStructures.h"
#include "VTKWriter.h"
#include "Solver.h"

#include <iostream>
#include <type_traits>

int main(int argc, char const *argv[]) 
{

    /*-------------------------------------------------------------------------------------*\
                                         Input Processing
    \*-------------------------------------------------------------------------------------*/

    // Input file from first command line argument
    std::string inputFilename;
    if (argc == 1) {
        std::cout << "Please enter input file name: " << "\n";
        std::cin >> inputFilename;
        std::cout << "\n";
        std::cin.ignore();
    } else if (argc == 2) {
        inputFilename = argv[1];
    } else {
        throw std::invalid_argument( "Invalid command line options." );
    }

    // User input data
    std::string inputFileRetryChoice;
    CFD::InputData inputData;
    while (true){
        try {

            std::cout << "Reading input file: '" + inputFilename + "' ..." << "\n\n";
            inputData = CFD::ReadInputData(inputFilename);
            std::cout << "Success!" << "\n\n";
            break;

        // } catch (std::runtime_error &e) {
        } catch (int &e) {

            // std::cout << "Failure reading input file! \n" 
            //           << e.what() 
            //           << "\n\n";

            std::cout << "Would you like to try again? (y/n)" 
                      << "\n";
            std::cin  >> inputFileRetryChoice;
            if (inputFileRetryChoice != "y") 
                exit(-1);
            std::cout << "\n";

            std::cout << "Please enter input file name: " << "\n";
            std::cin  >> inputFilename;
            std::cout << "\n";
            std::cin.ignore();

        }
    }
    std::cout << "Press enter to begin.";
    std::cin.ignore();
    std::cout << std::endl;


    /*-------------------------------------------------------------------------------------*\
                                           Solve
    \*-------------------------------------------------------------------------------------*/


    CFD::Mesh mesh(inputData);

    using FE = CFD::Fields::ENUMDATA;
    CFD::ArrayAllocator<FE> fields(mesh.nCells_x, mesh.nCells_y, mesh.nCells_z, {FE::U, FE::V, FE::W, FE::P});

    CFD::SweepSolve(fields, mesh, inputData);


    /*-------------------------------------------------------------------------------------*\
                                           Output
    \*-------------------------------------------------------------------------------------*/

    // Data to pass to writer
    VTK::dataType VTKDataType = VTK::DOUBLE;
    if (std::is_same<CFD::floatType, float>::value) {
        VTKDataType = VTK::FLOAT;
    } 
    VTK::VTKWriterConfig config( mesh.cellCenters_x.size(), mesh.cellCenters_y.size(), mesh.cellCenters_z.size(), VTKDataType);
        config.SetWriteMode("ascii");
        config.SetASCIIPrecision(8);
    VTK::gridVectorType gridVector = {mesh.cellCenters_x.data(), mesh.cellCenters_y.data(), mesh.cellCenters_z.data()};
    VTK::scalarMapType scalarMap = { };
    VTK::vectorMapType vectorMap = { };

    // Write output
    VTK::VTKWriter writer(gridVector, scalarMap, vectorMap, config);
    writer.WriteData("mesh.vtk", "3D Rectilinear Grid");


    /*-------------------------------------------------------------------------------------*\
                                           Testing
    \*-------------------------------------------------------------------------------------*/

    std::vector<CFD::Fields::ENUMDATA> fieldEnums = {CFD::Fields::ENUMDATA::U, CFD::Fields::ENUMDATA::V, CFD::Fields::ENUMDATA::P};
    CFD::ArrayAllocator<CFD::Fields::ENUMDATA>(mesh.nCells_x, mesh.nCells_y, mesh.nCells_z, fieldEnums);


    return 0;
}