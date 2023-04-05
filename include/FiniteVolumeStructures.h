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
    std::vector< ExtrapFactorsStruct > extrapFactors;    // extrapFactors[boundaryPatch]
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
        ArrayAllocator<TransportCoefficients, array3D> AU, AV, AW;          // Velocity coefficients
        ArrayAllocator<TransportCoefficients, array1D> AP;                  // Pressure coefficients
        array3D B;                                                          // Source Term
        std::vector< ArrayAllocator<TransportCoefficients, array1D> > diff; // Diffusion coefficients diff[Axis][TransportCoefficient]
        std::vector< floatType > boundaryDiff, boundaryP;                   // Constant terms that come from uniform BC boundaryDiff[BoundaryPatch]
        std::vector< array2D > boundaryVel;
    };

    struct ContinuityEquation {
        ContinuityEquation(const indexVector3 &);
        ArrayAllocator<TransportCoefficients, array1D> AU, AV, AW;          // Velocity coefficients
        ArrayAllocator<TransportCoefficients, array3D> AP;                  // Pressure coefficients
        array3D B;                                                          // Source term
        std::vector< floatType > boundaryVel;                               // Constant terms that come from uniform BC boundaryVel[BoundaryPatch]
        std::vector< floatType > boundaryP;                                 // These may need to be 2D arrays   
    };

    MomentumEquation Umom, Vmom, Wmom;
    ContinuityEquation Cont;
};


} // end namespace CFD

#endif // FV_STRUCTURES