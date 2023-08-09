#ifndef SOLTRANSFORM
#define SOLTRANSFORM

#include "Types.h"
#include "IO/InputProcessing.h"
#include "Tools/SweepTransformations.h"
#include "Tools/FVTools.h"

#include <iostream>
#include <tuple>

#include <vtkType.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridReader.h>
#include <vtkAOSDataArrayTemplate.h>


// Reads expected command line arguments
inline std::tuple<std::string, std::string, std::string> ReadCommandLineInputs(int argc, char const *argv[])
{
    if ( argc != 4 ) {
        throw std::invalid_argument("Invalid command line arguments.");
    } 
    std::string inputFilename            = argv[1], 
                originalFieldFilename    = argv[2], 
                transformedFieldFilename = argv[3];

    return { inputFilename, originalFieldFilename, transformedFieldFilename };
}



// Copy cell face data from vtkRectilinearGrid into Eigen Tensors
inline CFD::EnumVector<CFD::Axis, CFD::array1D> GetCellFaces( vtkRectilinearGrid *vtkGrid )
{
    using namespace CFD;

    int nFaces[3];
    vtkGrid->GetDimensions( nFaces );
    floatType *gridPointers[3] = { static_cast< floatType* >( vtkGrid->GetXCoordinates()->GetVoidPointer(0) ),
                                   static_cast< floatType* >( vtkGrid->GetYCoordinates()->GetVoidPointer(0) ),
                                   static_cast< floatType* >( vtkGrid->GetZCoordinates()->GetVoidPointer(0) ) };
    EnumVector<Axis, array1D > cellFaces;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        cellFaces[axis] = array1D( nFaces[axis] );
        std::memcpy( cellFaces[axis].data(), gridPointers[axis], nFaces[axis] * sizeof( floatType ) );
    } );

    return cellFaces;
}




CFD::EnumVector<CFD::Axis, CFD::array3D> GetVectorFieldFromVTKArray( vtkDataArray *dataArray,
                                                                     const int nCells[3] )
{
    using namespace CFD;

    EnumVector<Axis, array3D> vectorField( array3D( nCells[0], nCells[1], nCells[2] ) );
    intType nCellsTotal = static_cast<intType>( nCells[0] ) 
                        * static_cast<intType>( nCells[1] ) 
                        * static_cast<intType>( nCells[2] );

    for ( intType i = 0; i != nCellsTotal; i++ ) {

        // X component
        vectorField[Axis::X].data()[i] = dataArray->GetComponent( i, 0 );

        // Y component
        vectorField[Axis::Y].data()[i] = dataArray->GetComponent( i, 1 );

        // Z component
        vectorField[Axis::Z].data()[i] = dataArray->GetComponent( i, 2 );

    }

    return vectorField;
}



CFD::array3D GetScalarFieldFromVTKArray( vtkDataArray *dataArray,
                                         const int nCells[3] )
{
    using namespace CFD;

    array3D scalarField( nCells[0], nCells[1], nCells[2] );
    intType nCellsTotal = static_cast<intType>( nCells[0] ) 
                        * static_cast<intType>( nCells[1] ) 
                        * static_cast<intType>( nCells[2] );
    floatType *vtkArrayPointer = static_cast< floatType* >( dataArray->GetVoidPointer(0) );

    std::memcpy( scalarField.data(), vtkArrayPointer, nCellsTotal * sizeof( floatType ) );

    return scalarField;
}



// Copy vertex fields from vtkRectilinearGrid into Eigen Tensors
inline CFD::FieldData<CFD::array3D> GetFieldData( vtkDataSetAttributes *vtkDataSet,
                                                  const int nPoints[3] )
{
    using namespace CFD;
    
    FieldData<array3D> fieldData( array3D( nPoints[0], nPoints[1], nPoints[2] ) );

    for ( int i = 0; i < vtkDataSet->GetNumberOfArrays(); i++ ) {

        if ( strcmp(vtkDataSet->GetArrayName(i), "Velocity") ){
            std::cout << "Velocity found" << "\n";
            fieldData.U = GetVectorFieldFromVTKArray( vtkDataSet->GetArray("Velocity"), nPoints );
            continue;
        }
            

        if ( strcmp(vtkDataSet->GetArrayName(i), "Pressure")  ) {
            std::cout << "Pressure found" << "\n";
            fieldData.P = GetScalarFieldFromVTKArray( vtkDataSet->GetArray("Pressure"), nPoints );
            continue;
        }
            

    }

    return fieldData;
}



// Copy cell fields from vtkRectilinearGrid into Eigen Tensors
inline CFD::FieldData<CFD::array3D> GetCellFields( vtkRectilinearGrid *vtkGrid )
{
    using namespace CFD;

    int nCells[3];
    vtkGrid->GetCellDims( nCells );
    vtkCellData *cellData  = vtkGrid->GetCellData();
    FieldData<array3D> cellFields = GetFieldData( cellData, nCells );

    return cellFields;
}



// Copy vertex fields from vtkRectilinearGrid into Eigen Tensors
inline CFD::FieldData<CFD::array3D> GetVertexFields( vtkRectilinearGrid *vtkGrid )
{
    using namespace CFD;

    int nPoints[3];
    vtkGrid->GetDimensions( nPoints );
    vtkPointData *pointData = vtkGrid->GetPointData();
    FieldData<array3D> vertexFields = GetFieldData( pointData, nPoints );
    return vertexFields;
}


#endif  // SOLTRANSFORM