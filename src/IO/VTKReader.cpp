#ifdef CFD_HAS_VTK_LIB
#include "VTKReader.h"

#include <vtkType.h>
#include <vtkDataSetAttributes.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridReader.h>
#include <vtkAOSDataArrayTemplate.h>

namespace VTK
{

namespace 
{

// Copy cell face data from vtkRectilinearGrid into Eigen Tensors
inline CFD::EnumVector<CFD::Axis, CFD::Tensor1D> GetCellFaces( vtkRectilinearGrid *vtkGrid )
{
    using namespace CFD;

    int nFaces[3];
    vtkGrid->GetDimensions( nFaces );
    floatType *gridPointers[3] = { static_cast< floatType* >( vtkGrid->GetXCoordinates()->GetVoidPointer(0) ),
                                   static_cast< floatType* >( vtkGrid->GetYCoordinates()->GetVoidPointer(0) ),
                                   static_cast< floatType* >( vtkGrid->GetZCoordinates()->GetVoidPointer(0) ) };
    EnumVector<Axis, Tensor1D > cellFaces;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        cellFaces[axis] = Tensor1D( nFaces[axis] );
        std::memcpy( cellFaces[axis].data(), gridPointers[axis], nFaces[axis] * sizeof( floatType ) );
    } );

    return cellFaces;
}




CFD::EnumVector<CFD::Axis, CFD::Tensor3D> GetVectorFieldFromVTKArray( vtkDataArray *dataArray,
                                                                     const int nCells[3] )
{
    using namespace CFD;

    EnumVector<Axis, Tensor3D> vectorField( Tensor3D( nCells[0], nCells[1], nCells[2] ) );
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



CFD::Tensor3D GetScalarFieldFromVTKArray( vtkDataArray *dataArray,
                                         const int nCells[3] )
{
    using namespace CFD;

    Tensor3D scalarField( nCells[0], nCells[1], nCells[2] );
    intType nCellsTotal = static_cast<intType>( nCells[0] ) 
                        * static_cast<intType>( nCells[1] ) 
                        * static_cast<intType>( nCells[2] );
    floatType *vtkArrayPointer = static_cast< floatType* >( dataArray->GetVoidPointer(0) );

    std::memcpy( scalarField.data(), vtkArrayPointer, nCellsTotal * sizeof( floatType ) );

    return scalarField;
}



// Copy vertex fields from vtkRectilinearGrid into Eigen Tensors
inline CFD::FieldData<CFD::Tensor3D> GetFieldData( vtkDataSetAttributes *vtkDataSet,
                                                   const int nPoints[3] )
{
    using namespace CFD;
    
    FieldData<Tensor3D> fieldData( Tensor3D( nPoints[0], nPoints[1], nPoints[2] ) );

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
inline CFD::FieldData<CFD::Tensor3D> GetCellFields( vtkRectilinearGrid *vtkGrid )
{
    using namespace CFD;

    int nCells[3];
    vtkGrid->GetCellDims( nCells );
    vtkCellData *cellData  = vtkGrid->GetCellData();
    FieldData<Tensor3D> cellFields = GetFieldData( cellData, nCells );

    return cellFields;
}



// Copy vertex fields from vtkRectilinearGrid into Eigen Tensors
inline CFD::FieldData<CFD::Tensor3D> GetVertexFields( vtkRectilinearGrid *vtkGrid )
{
    using namespace CFD;

    int nPoints[3];
    vtkGrid->GetDimensions( nPoints );
    vtkPointData *pointData = vtkGrid->GetPointData();
    FieldData<Tensor3D> vertexFields = GetFieldData( pointData, nPoints );
    return vertexFields;
}


}   // end anonymous namespace


FieldFileData ReadVTKFields( const std::string &filename )
{
    using namespace CFD;
    
    FieldFileData fieldFileData;

    vtkNew< vtkRectilinearGridReader > vtkGridReader;
    vtkGridReader->SetFileName( filename.c_str() );
    vtkGridReader->Update();
    vtkRectilinearGrid* vtkGrid = vtkGridReader->GetOutput();
    
    // Make sure the data type in the file is the same type as the code
    if ( vtkGrid->GetScalarType() == VTK_DOUBLE ) {

        if ( !std::is_same<CFD::floatType, double>::value )
            throw std::runtime_error( "Type mismatch. VTK files to be read must be in double precision" );

    } else if ( vtkGrid->GetScalarType() == VTK_FLOAT ) {

        if ( !std::is_same<CFD::floatType, float>::value )
            throw std::runtime_error( "Type mismatch. VTK files to be read must be in single precision" );

    }

    fieldFileData.cellFaces    = GetCellFaces( vtkGrid );
    fieldFileData.cellFields   = GetCellFields( vtkGrid );
    fieldFileData.vertexFields = GetVertexFields( vtkGrid );

    return fieldFileData;
}


}   // end namespace VTK;

#endif // CFD_HAS_VTK_LIB