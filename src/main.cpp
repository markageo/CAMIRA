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

    return 0;
}