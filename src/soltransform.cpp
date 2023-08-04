/*---------------------------------------------------------------------------*\
   Read and transform a legacy format vtk field from given file

   Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"

#include "IO/InputProcessing.h"
#include "Tools/SweepTransformations.h"

#include <iostream>
#include <tuple>


std::tuple<std::string, std::string, std::string> ReadCommandLineInputs(int argc, char const *argv[])
{
    std::string inputFilename, originalFieldFilename, transformedFieldFilename;
    if ( argc != 3 ) {
        throw std::invalid_argument("Invalid command line options.");
    } 
    return { argv[1], argv[2], argv[3] };
}



int main(int argc, char const *argv[])
{

    using namespace CFD;

    // Read command line data
    auto [ inputFilename, originalFieldFilename, transformedFieldFilename ] = ReadCommandLineInputs( argc, argv );

    // Read just the axis transformation from the input file
    auto [planeSweepDirection, lineSweepDirection] = ReadSweepDirections( inputFilename );
    AxisTransformationMap axisTransformation = CreateAxisTransformation( planeSweepDirection, lineSweepDirection );

    // Read and store the field


    // Transform the field

    // Write the transformed field


    return 0;
}