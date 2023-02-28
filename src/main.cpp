/*---------------------------------------------------------------------------*\
    3D coupled Navier-Stokes solver on rectilinear grids

    Mark George
\*---------------------------------------------------------------------------*/

#include "SimulationParameters.h"
#include "InputProcessing.h"

#include <iostream>
#include <optional>

int main(int argc, char const *argv[]) 
{

    std::string inputFilename = "input.inp";
    auto inputData_optional = ReadInputData(inputFilename);
    if (!inputData_optional) 
        return -1;
    auto inputData = inputData_optional.value();

    std::cout << inputData.domainSizeX << std::endl;
    std::cout << inputData.domainSizeY << std::endl;
    std::cout << inputData.domainSizeZ << std::endl;
    std::cout << "\n";

    for (int i = 0; i != 3; i++) {
        std::cout << inputData.mesh.xBiasFactors[i]  << std::endl;
        std::cout << inputData.mesh.xCells[i]  << std::endl;
        std::cout << inputData.mesh.xSegmentBounds[i].first << " ";
        std::cout << inputData.mesh.xSegmentBounds[i].second  << std::endl;
    }
    std::cout << "\n";

    for (int i = 0; i != 1; i++) {
        std::cout << inputData.mesh.yBiasFactors[i]  << std::endl;
        std::cout << inputData.mesh.yCells[i]  << std::endl;
        std::cout << inputData.mesh.ySegmentBounds[i].first << " ";
        std::cout << inputData.mesh.ySegmentBounds[i].second  << std::endl;
    }
    std::cout << "\n";

    for (int i = 0; i != 1; i++) {
        std::cout << inputData.mesh.zBiasFactors[i]  << std::endl;
        std::cout << inputData.mesh.zCells[i]  << std::endl;
        std::cout << inputData.mesh.zSegmentBounds[i].first << " ";
        std::cout << inputData.mesh.zSegmentBounds[i].second  << std::endl;
    }
    std::cout << "\n";

    return 0;
}