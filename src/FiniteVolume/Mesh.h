#ifndef CFD_MESH  
#define CFD_MESH

#include "../Types.h"
#include "../Tools/AxisTransformationMap.h"
#include "../IO/InputProcessing.h"

namespace CFD
{

// Recitlinear mesh structure and mesher (on construction)
struct Mesh
{
    Mesh() {};
    Mesh(const iArray3 &);
    Mesh(const InputData &);
    Mesh(const EnumVector<Axis, Tensor1D> &, const InputData &);
    iArray3 nCells;
    EnumVector<Axis, iArray3> nFacesNormal;
    EnumVector<Axis, Tensor1D> cellCenters, 
                               cellFaces,            // cellFaces[axis](i) -> cellFaces[axis](i - 1/2)
                               cellLengths, 
                               cellLengthsInv,       // inverse of cell lengths
                               cellCenterDiffInv,    // inverse of distance between cell centers, same convention as cell faces
                               interpFactors;        // faceValue(i) = (1 - interpFactor(i))*cellValue(i-1) + interpFactor(i)*cellValue(i)
    EnumVector<Axis, Tensor2D> cellFaceAreas;        // Index by X, Y, Z order, not right hand rule.

    struct ExtrapFactorsStruct {
        floatType p,    // Boundary cell 
                  a;    // One from boundary cell
    };
    EnumVector< BoundaryPatches, ExtrapFactorsStruct > extrapFactors;
};


// Create mesh and output mesh information to console
Mesh CreateMesh( const InputData &,
                 const AxisTransformationMap & );

// Determines if mesh can be coarsened
bool MeshCanBeCoarsened( const Mesh& );

// Coarsens mesh by agglomorating cells
Mesh CoarsenMesh( const Mesh&, 
                  const FaceInterpolationSchemes ); 


} // end namespace CFD

#endif // CFD_MESH