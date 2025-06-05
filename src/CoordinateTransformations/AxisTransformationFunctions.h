#ifndef AXIS_TRANSFORMATION_FUNCTIONS
#define AXIS_TRANSFORMATION_FUNCTIONS

#include "AxisTransformationMap.h"
#include "../Core/Types.h"
#include "../IO/InputProcessing.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolumeStructures.h"

// The problem is remapped so that the plane sweeping direction is always in the z direction and
// the line sweeping direction is always in the y direction (in the code). This is more memory
// efficient and is simpler to implement.

namespace CAMIRA
{

// Return a transformation map from a plane and line sweeping direction
AxisTransformationMap CreateAxisTransformation( BoundaryPatches::ENUMDATA, BoundaryPatches::ENUMDATA ); 

// Transform user input data so that sweeping is consistent with solver
void TransformUserInputData( InputData &, const AxisTransformationMap & );

// Transform fields and mesh to code coordinates for when solution data is read into the solver
void TransformMeshToCodeCoordinates( Mesh &, const AxisTransformationMap &);
void TransformFieldToCodeCoordinates( FieldData<Tensor3D> &, const AxisTransformationMap &);

// Transform back to the coordinates consistentent with the user input
void TransformMeshToUserCoordinates( Mesh &, const AxisTransformationMap &);
void TransformScalarFieldToUserCoordinates( Tensor3D &, const AxisTransformationMap &);
void TransformVectorFieldToUserCoordinates( EnumVector<Axis, Tensor3D> &, const AxisTransformationMap &);
void TransformBCDataToUserCoordinates( BoundaryConditionData &, const AxisTransformationMap & );


} // end namespace CAMIRA
 
#endif  // AXIS_TRANSFORMATION_FUNCTIONS