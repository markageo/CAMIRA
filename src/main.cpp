/*---------------------------------------------------------------------------*\
    3D coupled Navier-Stokes solver on rectilinear grids

    Mark George
\*---------------------------------------------------------------------------*/

#include "SimulationParameters.h"
#include "InputProcessing.h"
#include "MeshStructure.h"


#include <iostream>
#include <optional>

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
                                           Testing
    \*-------------------------------------------------------------------------------------*/

    std::cout << inputData.domainSize_x << std::endl;
    std::cout << inputData.domainSize_y << std::endl;
    std::cout << inputData.domainSize_z << std::endl;
    std::cout << "\n";

    for (int i = 0; i != 3; i++) {
        std::cout << inputData.mesh.biasFactors_x[i]  << std::endl;
        std::cout << inputData.mesh.nCells_x[i]  << std::endl;
        std::cout << inputData.mesh.segmentBounds_x[i].first << " ";
        std::cout << inputData.mesh.segmentBounds_x[i].second  << std::endl;
    }
    std::cout << "\n";

    for (int i = 0; i != 1; i++) {
        std::cout << inputData.mesh.biasFactors_y[i]  << std::endl;
        std::cout << inputData.mesh.nCells_y[i]  << std::endl;
        std::cout << inputData.mesh.segmentBounds_y[i].first << " ";
        std::cout << inputData.mesh.segmentBounds_y[i].second  << std::endl;
    }
    std::cout << "\n";

    for (int i = 0; i != 1; i++) {
        std::cout << inputData.mesh.biasFactors_z[i]  << std::endl;
        std::cout << inputData.mesh.nCells_z[i]  << std::endl;
        std::cout << inputData.mesh.segmentBounds_z[i].first << " ";
        std::cout << inputData.mesh.segmentBounds_z[i].second  << std::endl;
    }
    std::cout << "\n";

    std::cout << "Mesh dimensions" << std::endl;
    std::cout << meshStructure.cellCenters_x.dimension(0) << std::endl;
    std::cout << meshStructure.cellCenters_y.dimension(0) << std::endl;
    std::cout << meshStructure.cellCenters_z.dimension(0) << std::endl;

    return 0;
}