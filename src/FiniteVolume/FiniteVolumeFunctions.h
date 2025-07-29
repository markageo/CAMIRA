#ifndef CAMIRA_FINITE_VOLUME_FUNCTIONS   
#define CAMIRA_FINITE_VOLUME_FUNCTIONS

#include "Mesh.h"
#include "FVCoefficients.h"
#include "../Core/Types.h"
#include "../IO/InputProcessing.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../TurbulenceModels/TurbulenceModelData.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"

namespace CAMIRA
{

// -------------------------------------- Definition in InitialConditions.cpp -------------------------------------- //

// Allocate and initialise the fields
FieldData<Tensor3D> InitialiseFields(const Mesh &, const InputData &, const AxisTransformationMap & );



// -------------------------------------- Definition in BoundaryConditions.cpp -------------------------------------- //

// Calculate and set boundary condition data for all fields
BoundaryConditionData SetBoundaryConditionData( const InputData &, const Mesh & );



// -------------------------------------- Definition in FaceVelocities.cpp -------------------------------------- //

// Allocate and initialise face velocities
EnumVector<Axis, Tensor3D> InitialiseFaceFluxes( const Mesh &, 
                                                 const EnumVector<Axis, Tensor3D> &, 
                                                 const BoundaryConditionData &);


// Update face velocities
void UpdateFaceFluxes( EnumVector<Axis, Tensor3D> &, 
                       const Mesh &, 
                       const EnumVector<Axis, Tensor3D> &, 
                       const BoundaryConditionData &);

void UpdateFaceFluxesWithMWI( EnumVector<Axis, Tensor3D> &, 
                              const Mesh &, 
                              const FieldData<Tensor3D> &,
                              const FVCoefficients &, 
                              const BoundaryConditionData &);


void SetIBFaceFluxes( EnumVector<Axis, Tensor3D> &, 
                      const IBData & );



// ---------------------------------- Definition in FiniteVolumeCoefficients.cpp --------------------------------- //

// Allocate and initialise finite volume coefficients
FVCoefficients InitialiseFVCoefficients( const Mesh &,
                                         const InputData &);

// Update finite volume coefficients 
void UpdateFVCoefficients( FVCoefficients &, 
                           const Mesh &, 
                           const FieldData< Tensor3D > &, 
                           const FieldData< Tensor3D > &, 
                           const FieldData< Tensor3D > &, 
                           const EnumVector< Axis, Tensor3D > &,
                           const IBData &,
                           const TurbulenceModelData & );


// ---------------------------------------- Definition in VertexValues.cpp -------------------------------------- //


FieldData<Tensor3D> GetVertexFields( const FieldData<Tensor3D> &, const Mesh &, const BoundaryConditionData & );



} // end namespace CAMIRA

#endif // CAMIRA_FINITE_VOLUME_FUNCTIONS