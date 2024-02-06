#ifndef CFD_FINITE_VOLUME   
#define CFD_FINITE_VOLUME

#include "Mesh.h"
#include "../Types.h"
#include "../IO/InputProcessing.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include <vector>

namespace CFD
{

// -------------------------------------- Definition in FiniteVolumeStructures.cpp -------------------------------------- //


struct MomentumEquation {
    MomentumEquation(const Axis::ENUMDATA, const iArray3 &, Linearisation);
    EnumVector< Axis, EnumVector< TransportCoefficients, Tensor3D > > AU;     // Velocity coefficients (LHS)
    EnumVector<TransportCoefficients, Tensor1D> AP;                           // Pressure coefficients (LHS)
    Tensor3D B;                                                               // Source Term (RHS)
    Tensor3D diagCoeffInv;                                                    // Inverse of diagonal coefficient
    EnumVector< Axis, EnumVector<TransportCoefficients, Tensor1D> > diff;     // Diffusion coefficients (LHS)
    EnumVector< BoundaryPatches, floatType > diffBoundary;                    // Diffusion coefficients for constant boundary conditions (LHS)
    EnumVector< BoundaryPatches, Tensor2D   > BUBoundary, BPBoundary;         // Constant terms that come from fixed BC (LHS)
    floatType relaxation;
    Axis::ENUMDATA component;                                                 // The momentum component
    Linearisation linearisation;
};


struct ContinuityEquation {
    ContinuityEquation(const iArray3 &, MomentumInterpolation);
    EnumVector< Axis, EnumVector< TransportCoefficients, Tensor1D > > AU;    // Velocity coefficients (LHS)
    EnumVector<TransportCoefficients, Tensor3D> AP;                          // Pressure coefficients (LHS)
    Tensor3D B;                                                              // Source term (RHS)
    EnumVector< Axis, std::array< Tensor1D, 4 > > mwiSparseCoeffs;           // Unweighted MWI coefficients from the sparse pressure gradient (LHS)
    EnumVector< Axis, std::array< Tensor1D, 2 > > mwiCompactCoeffs;          // Unweighted MWI coefficients from the compact pressure gradient (LHS)
    EnumVector< BoundaryPatches, Tensor2D   > BUBoundary, BPBoundary;        // Constant terms that come from fixed BC (LHS)
    floatType relaxation; 
    MomentumInterpolation momentumInterpolation;
};


// Structure to store all domain boundary condition information
struct BoundaryConditionData {
    struct Patch {
        BoundaryConditions::ENUMDATA type;
        Tensor2D value;
    };
    using Patches = EnumVector< BoundaryPatches, Patch >;
    FieldData< Patches > fields; 
    bool pressureFieldIsFloating;
};


// Structure to store finite volume discrete equation coefficients
struct FVCoefficients
{
    FVCoefficients(const iArray3 &, Linearisation, MomentumInterpolation);
    
    EnumVector<Axis, MomentumEquation > Mom;
    ContinuityEquation Cont;
    iArray3 nCells;
};


// Allocate and initialise the fields
FieldData<Tensor3D> InitialiseFields(const Mesh &, const InputData &);

// Calculate and set boundary condition data for all fields
BoundaryConditionData SetBoundaryConditionData( const InputData &, const Mesh & );



// -------------------------------------- Definition in FaceVelocities.cpp -------------------------------------- //

// Allocate and initialise face velocities
EnumVector<Axis, Tensor3D> InitialiseFaceFluxes( const Mesh &, 
                                                 const EnumVector<Axis, Tensor3D> &, 
                                                 const BoundaryConditionData &);

EnumVector< Axis, EnumVector<Axis, Tensor3D> > InitialiseAdvectedFaceVelocities( const Mesh &, 
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
                                   const EnumVector<Axis, Tensor3D> &, 
                                   const EnumVector<Axis, Tensor3D> &, 
                                   const BoundaryConditionData &);

void SetIBFaceFluxes( EnumVector<Axis, Tensor3D> &, 
                      const IBData & );




// ---------------------------------- Definition in FiniteVolumeCoefficients.cpp --------------------------------- //

// Allocate and initialise finite volume coefficients
FVCoefficients InitialiseFVCoefficients( const Mesh &, 
                                         const FieldData< Tensor3D > &, 
                                         const EnumVector< Axis, EnumVector< Axis, Tensor3D> > &,
                                         const EnumVector< Axis, Tensor3D > &,
                                         const IBData &, 
                                         const BoundaryConditionData &, 
                                         const InputData &);

// Update finite volume coefficients 
void UpdateFVCoefficients( FVCoefficients &, 
                           const Mesh &, 
                           const FieldData< Tensor3D > &, 
                           const EnumVector< Axis, EnumVector< Axis, Tensor3D> > &,
                           const EnumVector< Axis, Tensor3D > &,
                           const IBData &,
                           const BoundaryConditionData &);


// ---------------------------------------- Definition in VertexValues.cpp -------------------------------------- //


FieldData<Tensor3D> GetVertexFields( const FieldData<Tensor3D> &, const Mesh &, const BoundaryConditionData & );



} // end namespace CFD

#endif // CFD_FINITE_VOLUME