/*---------------------------------------------------------------------------*\
    3D coupled Navier-Stokes solver on rectilinear grids

    Mark George
\*---------------------------------------------------------------------------*/

#include "SimulationParameters.h"
#include "InputProcessing.h"
#include "MeshStructure.h"
#include "VTKWriter.h"


#include <iostream>
#include <optional>
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
    std::optional<InputData> inputData_optional;
    while (true){
        try {
            std::cout << "Reading input file: '" + inputFilename + "'" << std::endl;
            inputData_optional = ReadInputData(inputFilename);
            if (!inputData_optional) 
                throw -1;
            std::cout << "Success!" << "\n\n";
            break;
        } catch (int status) {
            std::cout << "Failure reading input file, would you like to try again? (y/n)" << "\n";
            std::cin >> inputFileRetryChoice;
            if (inputFileRetryChoice != "y") {exit(-1);}
            std::cout << "\n";

            std::cout << "Please enter input file name: " << "\n";
            std::cin >> inputFilename;
            std::cout << "\n";
            std::cin.ignore();
        }
    }
    std::cout << "Press enter to begin.";
    std::cin.ignore();
    std::cout << std::endl;
    auto inputData = inputData_optional.value();



    /*-------------------------------------------------------------------------------------*\
                                           Meshing
    \*-------------------------------------------------------------------------------------*/

    MeshStructure meshStructure(inputData);


    /*-------------------------------------------------------------------------------------*\
                                           Output
    \*-------------------------------------------------------------------------------------*/

    // Data to pass to writer
    VTK::dataType VTKDataType = VTK::DOUBLE;
    if (std::is_same<SIM::floatType, float>::value) {
        VTKDataType = VTK::FLOAT;
    } 
    VTK::VTKWriterConfig config( meshStructure.cellCenters_x.size(), meshStructure.cellCenters_y.size(), meshStructure.cellCenters_z.size(), VTKDataType);
        config.SetWriteMode("ascii");
        config.SetASCIIPrecision(8);
    VTK::gridVectorType gridVector = {meshStructure.cellCenters_x.data(), meshStructure.cellCenters_y.data(), meshStructure.cellCenters_z.data()};
    VTK::scalarMapType scalarMap = { };
    VTK::vectorMapType vectorMap = { };

    // Write output
    VTK::VTKWriter writer(gridVector, scalarMap, vectorMap, config);
    writer.WriteData("mesh.vtk", "3D Rectilinear Grid");


    /*-------------------------------------------------------------------------------------*\
                                           Testing
    \*-------------------------------------------------------------------------------------*/



    return 0;
}