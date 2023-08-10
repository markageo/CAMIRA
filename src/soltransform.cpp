/*---------------------------------------------------------------------------*\
   Read and transform a legacy format vtk field from given file

   Mark George
\*---------------------------------------------------------------------------*/

#include "soltransform.h"

#include "Types.h"
#include "IO/VTKWriter.h"

#include "IO/InputProcessing.h"
#include "Tools/SweepTransformations.h"
#include "Tools/FVTools.h"
#include <iostream>
#include <tuple>

#include <vtkType.h>
#include <vtkCellData.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridReader.h>
#include <vtkAOSDataArrayTemplate.h>


int main(int argc, char const *argv[])
{

    // Read command line data
    auto [ inputFilename, originalFieldFilename, transformedFieldFilename ] = ReadCommandLineInputs( argc, argv );

    // Read just the axis transformation from the input file
    auto [planeSweepDirection, lineSweepDirection] = CFD::ReadSweepDirections( inputFilename );
    CFD::AxisTransformationMap axisTransformation = CFD::CreateAxisTransformation( planeSweepDirection, lineSweepDirection );

    // Read and store the field
    vtkNew< vtkRectilinearGridReader > vtkReader;
    vtkReader->SetFileName( originalFieldFilename.c_str() );
    vtkReader->Update();
    vtkRectilinearGrid* vtkGrid = vtkReader->GetOutput();
    
    using namespace CFD;

    // Make sure the data type in the file is the same type as the code
    if ( vtkGrid->GetScalarType() == VTK_DOUBLE ) {

        if ( !std::is_same<CFD::floatType, double>::value )
            throw std::runtime_error( "Type mismatch. VTK file to transform must be in double precision" );

    } else if ( vtkGrid->GetScalarType() == VTK_FLOAT ) {

        if ( !std::is_same<CFD::floatType, float>::value )
            throw std::runtime_error( "Type mismatch. VTK file to transform must be in single precision" );

    }

    // Copy fields into Eigen Tensors
    EnumVector<Axis, array1D > cellFaces = GetCellFaces( vtkGrid );
    FieldData<array3D> cellFields        = GetCellFields( vtkGrid );
    FieldData<array3D> vertexFields      = GetVertexFields( vtkGrid );

    // Transform the field
    // TODO

    // Put the transformed field into the Rectilinear grid object
    

    // Write the transformed field
    vtkNew< vtkRectilinearGridWriter > vtkWriter;
    vtkWriter->SetFileName( transformedFieldFilename.c_str() );
    vtkWriter->SetInputData( vtkGrid );
    vtkWriter->Write();


    // ---------------------------------- TESTING

    // std::string testFilename = "test_field.vtk";
    // CFD::intType ni{10}, nj{11}, nk{12};
    // // CFD::array3D arr(ni, nj, nk);
    // CFD::array1D arr(ni);
    // arr.setRandom();

    // vtkNew< vtkAOSDataArrayTemplate< CFD::floatType > > vtkArr;
    // vtkArr->SetArray( arr.data(), arr.size(), 1);

    // for ( CFD::intType i = 0; i != arr.size(); i++ ) {
    //     std::cout << "Eigen: " << arr(i) << "   "
    //               << "VTK  : " << vtkArr->GetValue(i)
    //               << "\n";

    // }


    VTK::VTKWriterConfig config(cellFaces[Axis::X].size(), cellFaces[Axis::Y].size(), cellFaces[Axis::Z].size());
    config.SetWriteMode(VTK::WriteModes::ASCII);
    VTK::gridVectorType<CFD::floatType> gridVector = {cellFaces[Axis::X].data(), cellFaces[Axis::Y].data(), cellFaces[Axis::Z].data()};

    VTK::scalarCollectionType<CFD::floatType> scalarMap = {{"Pressure", VTK::GridTypes::CELL_DATA, cellFields.P.data()}};

    VTK::vectorCollectionType<CFD::floatType> vectorMap = {{"Velocity", VTK::GridTypes::CELL_DATA, {cellFields.U[Axis::X].data(), cellFields.U[Axis::Y].data(), cellFields.U[Axis::Z].data()}}};
    
    VTK::VTKWriter writer(gridVector, scalarMap, vectorMap, config);


    writer.WriteData("output/test_rewrite.vtk", "CFD simulation output");


    // ------------------------------------------


    return 0;
}