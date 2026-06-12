#include "FiniteVolumeFunctions.h"

#include "Core/FVTools.h"
#include "Core/FVLookups.h"
#include "Core/IO/VTKReader.h"
#include "Flow/InputProcessing/InputProcessing.h"
#include "Flow/CoordinateTransformations/AxisTransformationFunctions.h"


#include <cmath>
#include <stdexcept>

namespace CAMIRA
{

using namespace CORE;

namespace FLOW
{

#ifdef CAMIRA_HAS_VTK_LIB
void SetInitialConditionFromVTKFile( FieldData<Tensor3D> &fields,
                                     const std::string &filename,
                                     const Mesh &mesh,
                                     const AxisTransformationMap &axisTransformation )
{
    VTK::FieldFileData fieldFileData = VTK::ReadVTKFields( filename );
    
    Mesh inputMesh( fieldFileData.cellFaces );
    TransformMeshToCodeCoordinates( inputMesh, axisTransformation);
    TransformFieldToCodeCoordinates( fieldFileData.cellFields, axisTransformation);

    TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                  extents = {mesh.nCells(0), mesh.nCells(1), mesh.nCells(2)};

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
}
#endif



void SetInitialConditionUniform( FieldData<Tensor3D> &fields,
                                 const FieldData<floatType> &constantInitialConditions,
                                 const Mesh &mesh )
{
    TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                  extents = {mesh.nCells(0), mesh.nCells(1), mesh.nCells(2)};

    ForAllFieldData( [&] (intType i) { 
        fields[i].slice( offsets, extents ).setConstant( constantInitialConditions[i] );  
    } );
}



void InitialiseFields( FieldData<Tensor3D> &fields,
                       const Mesh &mesh, 
                       const InputData &inputData,
      [[maybe_unused]] const AxisTransformationMap &axisTransformation )
{
    #if defined( CAMIRA_HAS_VTK_LIB )
        switch ( inputData.initialConditionType ) {
            case InputData::InitialConditionTypes::uniform:
                SetInitialConditionUniform( fields, inputData.constantInitialConditions, mesh );
                break;

            case InputData::InitialConditionTypes::vtkFile:
                SetInitialConditionFromVTKFile( fields, inputData.initialConditionsFieldFilename, mesh, axisTransformation );
                break;
        }
    #else
        SetInitialConditionUniform( fields, inputData.constantInitialConditions, mesh );
    #endif
}


}   // end namespace FLOW

}   // end namespace CAMIRA