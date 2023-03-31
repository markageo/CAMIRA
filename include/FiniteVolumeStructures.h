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
                                //   cellCenterDiffInv,    // inverse of distance between cell centers, same convention as cell faces
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
    // Naming convention:
    //  aev - coefficient for equation 'e' multiplying with variable 'v'
    //  be  - source term for equation 'e'
    // 
    // 'e' can take the values
    //  u: U momentum, v: V momentum, w: W momentum, c: Conitnuity
    //
    // 'v' can take the values
    //  u: x velocity, v: v velocity, w: w velocity, p: pressure

    // In the finite volume formulation, all equations are divided by the cell volume. This 
    // means that the pressure coefficients in the momentum equations and the velocity 
    // coefficients in the continuity equations can be stored in 1D arrays when using a 
    // rectilinear grid.

    FVCoefficients(const CFD::indexVector3 &);
    ArrayAllocator<TransportCoefficients, array3D> auu, avv, aww;          // Momentum velocity coefficients
    ArrayAllocator<TransportCoefficients, array1D> aup, avp, awp;          // Momentum pressure coefficients
    ArrayAllocator<TransportCoefficients, array1D> acu, acv, acw;          // Continuity velocity coefficients
    ArrayAllocator<TransportCoefficients, array3D> acp;                    // Continuity pressure coefficients
    array3D                                        bu, bv, bw, bc;         // Momentum and continuity source terms (appear on the right hand size)

    std::vector< ArrayAllocator<TransportCoefficients, array1D> > diffu, diffv, diffw;    // Diffusion coefficients, diffu[Axis][TransportCoefficient]
};


} // end namespace CFD

#endif // FV_STRUCTURES