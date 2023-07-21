#ifndef FINITE_VOLUME   
#define FINITE_VOLUME

#include "Types.h"
#include "InputProcessing.h"
#include <vector>

namespace CFD
{

// -------------------------------------- Definition in FiniteVolumeStructures.cpp -------------------------------------- //

// Recitlinear mesh structure and mesher (on construction)
struct Mesh
{
    Mesh(const InputData &);
    iVector3 nCells;
    EnumVector<Axis, array1D> cellCenters, 
                              cellFaces,            // cellFaces[axis](i) -> cellFaces[axis](i - 1/2)
                              cellLengths, 
                              cellLengthsInv,       // inverse of cell lengths
                              cellCenterDiffInv,    // inverse of distance between cell centers, same convention as cell faces
                              interpFactors;        // faceValue(i) = (1 - interpFactor(i))*cellValue(i-1) + interpFactor(i)*cellValue(i)
    EnumVector<Axis, array2D> cellFaceAreas;        // Index by X, Y, Z order, not right hand rule.

    struct ExtrapFactorsStruct {
        floatType p,    // Boundary cell 
                  a;    // One from boundary cell
    };
    EnumVector< BoundaryPatches, ExtrapFactorsStruct > extrapFactors;
};



struct MomentumEquation {
    MomentumEquation(const Axis::ENUMDATA, const iVector3 &);
    EnumVector< Axis, EnumVector< TransportCoefficients, array3D > > AU;     // Velocity coefficients (LHS)
    EnumVector<TransportCoefficients, array1D> AP;                           // Pressure coefficients (LHS)
    array3D B;                                                               // Source Term (RHS)
    array3D diagCoeffInv;                                                    // Inverse of diagonal coefficient
    EnumVector< Axis, EnumVector<TransportCoefficients, array1D> > diff;     // Diffusion coefficients (LHS)
    EnumVector< BoundaryPatches, floatType > diffBoundary;                   // Diffusion coefficients for constant boundary conditions (LHS)
    EnumVector< BoundaryPatches, array2D   > BUBoundary, BPBoundary;         // Constant terms that come from fixed BC (LHS)
    floatType relaxation;
    Axis::ENUMDATA component;                                                // The momentum component
};



template< MomentumInterpolation MI >
struct ContinuityEquation {
    ContinuityEquation(const iVector3 &);
    EnumVector< Axis, EnumVector< TransportCoefficients, array1D > > AU;    // Velocity coefficients (LHS)
    EnumVector<TransportCoefficients, array3D> AP;                          // Pressure coefficients (LHS)
    array3D B;                                                              // Source term (RHS)
    EnumVector< Axis, std::array< array1D, 4 > > mwiSparseCoeffs;           // Unweighted MWI coefficients from the sparse pressure gradient (LHS)
    EnumVector< Axis, std::array< array1D, 2 > > mwiCompactCoeffs;          // Unweighted MWI coefficients from the compact pressure gradient (LHS)
    EnumVector< BoundaryPatches, array2D   > BUBoundary, BPBoundary;        // Constant terms that come from fixed BC (LHS)
    floatType relaxation; 
};


struct BoundaryConditionConfig {
    BoundaryConditions::ENUMDATA type;
    array2D value;
};
using BoundaryConditionData = EnumVector< BoundaryPatches, BoundaryConditionConfig >;



// Structure to store finite volume discrete equation coefficients (Picard linearisation)
template< MomentumInterpolation MI >
struct FVCoefficients
{
    // In the finite volume formulation, all equations are divided by the cell volume. This 
    // means that the pressure coefficients in the momentum equations and the velocity 
    // coefficients in the continuity equations can be stored in 1D arrays when using a 
    // rectilinear grid.

    FVCoefficients(const iVector3 &);
    
    EnumVector<Axis, MomentumEquation> Mom;
    ContinuityEquation< MI > Cont;
    iVector3 nCells;
    
};


// Allocate and initialise the fields
FieldData<array3D> InitialiseFields(const Mesh &, const InputData &);

// Remove ghost cells from a 3D array
void RemoveGhostCells( array3D &, const intType);

// Calculate and set boundary condition data for all fields
FieldData< BoundaryConditionData > SetBoundaryConditionData( const InputData &, const Mesh & );



// -------------------------------------- Definition in FaceVelocities.cpp -------------------------------------- //

// Allocate and initialise face velocities
EnumVector<Axis, array3D> InitialiseFaceFluxes(const Mesh &, const EnumVector<Axis, array3D> &, const FieldData< BoundaryConditionData > &);

// Update face velocities
void UpdateFaceFluxes( EnumVector<Axis, array3D> &, const Mesh &, const EnumVector<Axis, array3D> &, const FieldData< BoundaryConditionData > &);





// ---------------------------------- Definition in FiniteVolumeCoefficients.cpp --------------------------------- //

// Allocate and initialise finite volume coefficients
template< MomentumInterpolation MI >
FVCoefficients<MI> InitialiseFVCoefficients( const Mesh &, 
                                             const FieldData< array3D > &, 
                                             const EnumVector< Axis, array3D > &, 
                                             const FieldData< BoundaryConditionData > &, 
                                             const InputData &);

// Update finite volume coefficients 
template< MomentumInterpolation MI >
void UpdateFVCoefficients( FVCoefficients<MI> &, 
                           const Mesh &, 
                           const FieldData< array3D > &, 
                           const EnumVector< Axis, array3D > &,
                           const FieldData< BoundaryConditionData > &);


// ---------------------------------------- Definition in VertexValues.cpp -------------------------------------- //


FieldData<array3D> GetVertexFields( const FieldData<array3D> &, const Mesh &, const FieldData< BoundaryConditionData > & );



} // end namespace CFD

#endif // FINITE_VOLUME