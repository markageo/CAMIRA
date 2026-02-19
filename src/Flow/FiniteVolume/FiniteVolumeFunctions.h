#ifndef CAMIRA_FINITE_VOLUME_FUNCTIONS   
#define CAMIRA_FINITE_VOLUME_FUNCTIONS

#include "FVCoefficients.h"
#include "Core/Mesh/Mesh.h"
#include "Core/Types.h"
#include "Flow/InputProcessing/InputProcessing.h"
#include "Flow/ImmersedBoundary/ImmersedBoundary.h"
#include "Flow/TurbulenceModels/TurbulenceModelData.h"
#include "Flow/CoordinateTransformations/AxisTransformationFunctions.h"

namespace CAMIRA
{

using namespace CORE;

namespace FLOW
{

// -------------------------------------- Definition in InitialConditions.cpp -------------------------------------- //

// Allocate and initialise the fields
void InitialiseFields( FieldData<Tensor3D> &, const Mesh &, const InputData &, const AxisTransformationMap & );



// -------------------------------------- Definition in BoundaryConditions.cpp -------------------------------------- //

// Calculate and set boundary condition data for all fields
void SetBoundaryConditionData( BoundaryConditionData &, const InputData &, const Mesh & );



// -------------------------------------- Definition in FaceVelocities.cpp -------------------------------------- //

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
void InitialiseFVCoefficients( FVCoefficients &,
                               const Mesh &,
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


FieldData<Tensor3D> GetVertexFields( const FieldData<Tensor3D> &, 
                                     const Mesh &, 
                                     const BoundaryConditionData &, 
                                     const Tensor3D & );


}   // end namespace FLOW

}   // end namespace CAMIRA

#endif // CAMIRA_FINITE_VOLUME_FUNCTIONS