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
    Mesh(const CFD::InputData &);
    indexVector3 nCells;
    ArrayAllocator<Axis, array1D> cellCenters, 
                                  cellFaces,            // cellFaces[axis](i) -> cellFaces[axis](i - 1/2)
                                  cellLengths, 
                                  cellLengthsInv,       // inverse of cell lengths
                                  cellCenterDiffInv,    // inverse of distance between cell centers, same convention as cell faces
                                  interpFactors;        // faceValue(i) = (1 - interpFactor(i))*cellValue(i-1) + interpFactor(i)*cellValue(i)

    struct ExtrapFactorsStruct {
        floatType p,    // Boundary cell 
                  a;    // One from boundary cell
    };
    EnumVector< BoundaryPatches, ExtrapFactorsStruct > extrapFactors;
};


// Structure to store finite volume discrete equation coefficients (Picard linearisation)
struct FVCoefficients
{
    // In the finite volume formulation, all equations are divided by the cell volume. This 
    // means that the pressure coefficients in the momentum equations and the velocity 
    // coefficients in the continuity equations can be stored in 1D arrays when using a 
    // rectilinear grid.

    FVCoefficients(const indexVector3 &);

    struct MomentumEquation {
        MomentumEquation(const Fields::ENUMDATA, const indexVector3 &);
        ArrayAllocator<TransportCoefficients, array3D> AU, AV, AW;          // Velocity coefficients (LHS)
        ArrayAllocator<TransportCoefficients, array1D> AP;                  // Pressure coefficients (LHS)
        array3D B;                                                          // Source Term (RHS)
        EnumVector< Axis, ArrayAllocator<TransportCoefficients, array1D> > diff;    // Diffusion coefficients (LHS)
        EnumVector< BoundaryPatches, floatType > boundaryDiff, boundaryP;           // Constant terms that come from uniform BC (LHS)
        EnumVector< BoundaryPatches, array2D > boundaryVel;       
    };

    struct ContinuityEquation {
        ContinuityEquation(const indexVector3 &);
        ArrayAllocator<TransportCoefficients, array1D> AU, AV, AW;          // Velocity coefficients (LHS)
        ArrayAllocator<TransportCoefficients, array3D> AP;                  // Pressure coefficients (LHS)
        array3D B;                                                          // Source term (RHS)
        EnumVector< BoundaryPatches, array2D > boundaryP;                   // Constant terms that come from uniform BC (LHS)
        EnumVector< BoundaryPatches, floatType > boundaryVel;
    };

    MomentumEquation Umom, Vmom, Wmom;
    ContinuityEquation Cont;
};


// Allocate and initialise the fields
ArrayAllocator<Fields, array3D> InitialiseFields(const Mesh &, const InputData &);


// Transform back to the coordinates consistentent with the input file 
void TransformToUserCoordinates(Mesh &, ArrayAllocator<Fields, array3D> &, const InputData::AxisTransformationMap &);


// Remove ghost cells from a 3D array
void RemoveGhostCells( array3D &, const intType);





// -------------------------------------- Definition in FaceVelocities.cpp -------------------------------------- //

// Allocate and initialise face velocities
ArrayAllocator<Fields, array3D> InitialiseFaceVelocities(const Mesh &, const ArrayAllocator<Fields, array3D> &, const InputData &);

// Update face velocities
void UpdateFaceVelocities( ArrayAllocator<Fields, array3D> &, const Mesh &, const ArrayAllocator<Fields, array3D> &, const InputData &);





// ---------------------------------- Definition in FiniteVolumeCoefficients.cpp --------------------------------- //

// Allocate and initialise finite volume coefficients
FVCoefficients InitialiseFVCoefficients(const Mesh &, const ArrayAllocator<Fields, array3D> &, const ArrayAllocator<Fields, array3D> &, const InputData &);

// Update finite volume coefficients (Picard linearisation)
void UpdateFVCoefficients(FVCoefficients &, const Mesh &, const ArrayAllocator<Fields, array3D> &, const ArrayAllocator<Fields, array3D> &, const InputData &);





} // end namespace CFD

#endif // FINITE_VOLUME