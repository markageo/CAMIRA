#include "FiniteVolumeFunctions.h"
#include "../IO/InputProcessing.h"
#include "../Core/FVTools.h"
#include "../Core/FVLookups.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"
#include "../IO/VTKReader.h"

#include <cmath>
#include <stdexcept>

namespace CAMIRA
{

#ifdef CAMIRA_HAS_VTK_LIB
FieldData<Tensor3D> SetInitialConditionFromVTKFile( const std::string &filename,
                                                    const Mesh &mesh,
                                                    const AxisTransformationMap &axisTransformation,
                                                    const InputData &inputData )
{
    VTK::FieldFileData fieldFileData = VTK::ReadVTKFields( filename );
    
    Mesh inputMesh( fieldFileData.cellFaces, inputData );
    TransformMeshToCodeCoordinates( inputMesh, axisTransformation);
    TransformFieldToCodeCoordinates( fieldFileData.cellFields, axisTransformation);

    TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                  extents = {mesh.nCells(0), mesh.nCells(1), mesh.nCells(2)};

    FieldData<Tensor3D> fields( Tensor3D( mesh.nCells(0) + 2*CAMIRA::nGhost, 
                                          mesh.nCells(1) + 2*CAMIRA::nGhost, 
                                          mesh.nCells(2) + 2*CAMIRA::nGhost).setZero() );

    // Careful! Just checking the mesh is the same size, however it is possible that cell centers are at different locations.
    // Ideally should add the ability to do an interpolation.
    if ( ( inputMesh.nCells(0) != mesh.nCells(0) ) ||
         ( inputMesh.nCells(1) != mesh.nCells(1) ) ||
         ( inputMesh.nCells(2) != mesh.nCells(2) )  ) {
            throw std::runtime_error( "Mesh dimensions for initial condition do not match!" );
    }
    
    ForAllFieldData( [&] (intType f) {
        fields[f].slice( offsets, extents ) = fieldFileData.cellFields[f];
    } );
    
    return fields;
}
#endif



FieldData<Tensor3D> SetInitialConditionUniform( const FieldData<floatType> &constantInitialConditions,
                                                const Mesh &mesh )
{
    TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                  extents = {mesh.nCells(0), mesh.nCells(1), mesh.nCells(2)};

    FieldData<Tensor3D> fields( Tensor3D( mesh.nCells(0) + 2*CAMIRA::nGhost, 
                                          mesh.nCells(1) + 2*CAMIRA::nGhost, 
                                          mesh.nCells(2) + 2*CAMIRA::nGhost).setZero() );

    ForAllFieldData( [&] (intType i) { 
        fields[i].slice( offsets, extents ).setConstant( constantInitialConditions[i] );  
    } );
    
    return fields;
}



FieldData<Tensor3D> InitialiseFields( const Mesh &mesh, 
                                      const InputData &inputData,
                     [[maybe_unused]] const AxisTransformationMap &axisTransformation )
{
    FieldData<Tensor3D> fields;

    #if defined( CAMIRA_HAS_VTK_LIB )
        switch ( inputData.initialConditionType ) {
            case InputData::InitialConditionTypes::uniform:
                fields = SetInitialConditionUniform( inputData.constantInitialConditions, mesh );
                break;

            case InputData::InitialConditionTypes::vtkFile:
                fields = SetInitialConditionFromVTKFile( inputData.initialConditionsFieldFilename, mesh, axisTransformation, inputData );
                break;
        }
    #else
        fields = SetInitialConditionUniform( inputData.constantInitialConditions, mesh );
    #endif

    return fields;
}


}   // end namespace CAMIRA