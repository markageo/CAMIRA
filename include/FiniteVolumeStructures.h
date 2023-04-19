#ifndef FV_STRUCTURES
#define FV_STRUCTURES

#include "Types.h"
#include "InputProcessing.h"
#include <vector>

namespace CFD
{

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
    ArrayAllocator<Axis, array2D> cellFaceAreas;        // Index by right hand rule

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
        EnumVector< BoundaryPatches, floatType > boundaryDiff, boundaryP;             // Constant terms that come from uniform BC (LHS)
        EnumVector< BoundaryPatches, array2D > boundaryVel;       

    };

    struct ContinuityEquation {
        ContinuityEquation(const indexVector3 &);
        ArrayAllocator<TransportCoefficients, array1D> AU, AV, AW;          // Velocity coefficients (LHS)
        ArrayAllocator<TransportCoefficients, array3D> AP;                  // Pressure coefficients (LHS)
        array3D B;                                                          // Source term (RHS)
        EnumVector< BoundaryPatches, floatType > boundaryVel, boundaryP;    // Constant terms that come from uniform BC (LHS)
    };

    MomentumEquation Umom, Vmom, Wmom;
    ContinuityEquation Cont;
};


} // end namespace CFD

#endif // FV_STRUCTURES