#ifndef CFD_FINITE_VOLUME_FUNCTIONS   
#define CFD_FINITE_VOLUME_FUNCTIONS

#include "Mesh.h"
#include "FiniteVolumeStructures.h"
#include "../Core/Types.h"
#include "../IO/InputProcessing.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"

namespace CFD
{

// -------------------------------------- Definition in FiniteVolumeStructures.cpp -------------------------------------- //

// Allocate and initialise the fields
FieldData<Tensor3D> InitialiseFields(const Mesh &, const InputData &, const AxisTransformationMap & );

// Calculate and set boundary condition data for all fields
BoundaryConditionData SetBoundaryConditionData( const InputData &, const Mesh & );



// -------------------------------------- Definition in FaceVelocities.cpp -------------------------------------- //

// Allocate and initialise face velocities
EnumVector<Axis, Tensor3D> InitialiseFaceFluxes( const Mesh &, 
                                                 const EnumVector<Axis, Tensor3D> &, 
                                                 const BoundaryConditionData &);

EnumVector< Axis, EnumVector<Axis, Tensor3D> > InitialiseFaceAdvectedVelocities( const Mesh &, 
                                                                                 const FVCoefficients &,
                                                                                 const EnumVector<Axis, Tensor3D> &, 
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

void UpdateFaceAdvectedVelocities( EnumVector< Axis, EnumVector<Axis, Tensor3D> > &, 
                                   const Mesh &, 
                                   const FVCoefficients &,
                                   const EnumVector<Axis, Tensor3D> &, 
                                   const EnumVector<Axis, Tensor3D> &, 
                                   const BoundaryConditionData &);

void SetIBFaceFluxes( EnumVector<Axis, Tensor3D> &, 
                      const IBData & );

template<AdvectionSchemes>
void SetIBFaceAdvectedVelocities( EnumVector< Axis, EnumVector<Axis, Tensor3D> > &,
                                  const EnumVector< Axis, Tensor3D > &,
                                  const FieldData<Tensor3D> &,
                                  const FVCoefficients &,
                                  const Mesh &,
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
                           const IBData & );



void ApplyImplicitRelaxation( FVCoefficients &, 
                              const FieldData<Tensor3D> &,
                              const Mesh & );


// ---------------------------------------- Definition in VertexValues.cpp -------------------------------------- //


FieldData<Tensor3D> GetVertexFields( const FieldData<Tensor3D> &, const Mesh &, const BoundaryConditionData & );



} // end namespace CFD

#endif // CFD_FINITE_VOLUME_FUNCTIONS