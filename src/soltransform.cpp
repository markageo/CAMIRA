/*---------------------------------------------------------------------------*\
   Read and transform a legacy format vtk field from given file

   Mark George
\*---------------------------------------------------------------------------*/

#include "Types.h"

#include "IO/InputProcessing.h"
#include "Tools/SweepTransformations.h"

#include <iostream>
#include <tuple>


#include <vtkRectilinearGridReader.h>
#include <vtkAOSDataArrayTemplate.h>


std::tuple<std::string, std::string, std::string> ReadCommandLineInputs(int argc, char const *argv[])
{
    std::string inputFilename, originalFieldFilename, transformedFieldFilename;
    if ( argc != 4 ) {
        throw std::invalid_argument("Invalid command line arguments.");
    } 
    return { argv[1], argv[2], argv[3] };
}



int main(int argc, char const *argv[])
{

    // Read command line data
    auto [ inputFilename, originalFieldFilename, transformedFieldFilename ] = ReadCommandLineInputs( argc, argv );

    // Read just the axis transformation from the input file
    auto [planeSweepDirection, lineSweepDirection] = CFD::ReadSweepDirections( inputFilename );
    CFD::AxisTransformationMap axisTransformation = CFD::CreateAxisTransformation( planeSweepDirection, lineSweepDirection );

    // Read and store the field
    vtkNew< vtkRectilinearGridReader > gridReader;
    gridReader->SetFileName( originalFieldFilename.c_str() );
    gridReader->Update();

    // Transform the field

    // Write the transformed field


    // ---------------------------------- TESTING

    std::string testFilename = "test_field.vtk";
    CFD::intType ni{10}, nj{11}, nk{12};
    // CFD::array3D arr(ni, nj, nk);
    CFD::array1D arr(ni);
    arr.setRandom();

    vtkNew< vtkAOSDataArrayTemplate< double > > vtkArr;
    vtkArr->SetArray( arr.data(), arr.size(), 1);

    for ( CFD::intType i = 0; i != arr.size(); i++ ) {
        std::cout << "Eigen: " << arr(i) << "   "
                  << "VTK  : " << vtkArr->GetValue(i)
                  << "\n";

    }


    // ------------------------------------------


    return 0;
}