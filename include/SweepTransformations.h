#ifndef SWEEP_TRANSFORMATIONS
#define SWEEP_TRANSFORMATIONS

#include "Types.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"

// The problem is remapped so that the plane sweeping direction is always in the z direction and
// the line sweeping direction is always in the y direction (in the code). This is more memory
// efficient and is simpler to implement.

namespace CFD
{

// Transform user input data so that sweeping is consistent with solver
InputData TransformUserInputData( const InputData& );

// Transform back to the coordinates consistentent with the input file 
void TransformToUserCoordinates(Mesh &, ArrayAllocator<Fields, array3D> &, const InputData::AxisTransformationMap &);


} // end namespace CFD

#endif  // SWEEP_TRANSFORMATIONS