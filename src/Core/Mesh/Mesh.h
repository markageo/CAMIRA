#ifndef CAMIRA_MESH  
#define CAMIRA_MESH

#include "Core/Types.h"
#include "Core/AxisTransformationMap/AxisTransformationMap.h"

namespace CAMIRA
{

// Stores data for a segment of the mesh that would come from use input
struct MeshSegment {
    floatType startCoordinate;
    floatType endCoordinate;
    intType nCells;
    floatType biasFactor;
};


// Recitlinear mesh structure and mesher (on construction)
struct Mesh
{
    Mesh() {};
    Mesh(const iArray3 &);
    Mesh(const EnumVector< Axis, std::vector< MeshSegment > > & );
    Mesh(const EnumVector<Axis, Tensor1D> & );
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
Mesh CreateMesh( const EnumVector< Axis, std::vector< MeshSegment > > &,
                 const AxisTransformationMap & );

// Determines if mesh can be coarsened
bool MeshCanBeCoarsened( const Mesh & );

// Coarsens mesh by agglomorating cells
Mesh CoarsenMesh( const Mesh & ); 


} // end namespace CAMIRA

#endif // CAMIRA_MESH