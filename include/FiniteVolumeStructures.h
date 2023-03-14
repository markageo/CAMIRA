#ifndef FV_STRUCTURES
#define FV_STRUCTURES

#include "SimulationParameters.h"
#include "InputProcessing.h"
#include "Tensor"

namespace CFD 
{

namespace TC
{

enum TransportCoefficients {
    p,    // (i  , j  , k  ) 
    n,    // (i  , j+1, k  )
    e,    // (i+1, j  , k  ) 
    s,    // (i  , j-1, k  ) 
    w,    // (i-1, j  , k  ) 
    t,    // (i  , j  , k+1)
    b,    // (i  , j  , k-1)
    nn,   // (i  , j+2, k  ) 
    ee,   // (i+2, j  , k  )
    ss,   // (i  , j-2, k  )
    ww,   // (i-2, j  , k  )
    tt,   // (i  , j  , k+2)
    bb,   // (i  , j  , k-2)
    count
}; 

}   // end namespace TC


// Allocates and points to general coefficients in transport equation
class TransportEquation
{   
    public:

        TransportEquation(const intType n_x, const intType n_y, const intType n_z, const std::vector<TC::TransportCoefficients> coeffs) : 
            coeffsUsed(coeffs), coeffPointers(TC::count) 
        {
            for (const auto& index : coeffsUsed){
                coeffPointers[index] = new CFD::array3D(n_x, n_y, n_z);
            }
        }

        ~TransportEquation() {
            for (const auto& index : coeffsUsed) {
                delete coeffPointers[index];
            }
        }

        CFD::array3D& operator[](TC::TransportCoefficients idx){
            return *coeffPointers[idx];
        }

    private:

        std::vector<TC::TransportCoefficients> coeffsUsed;
        std::vector<CFD::array3D*> coeffPointers;
};


// Recitlinear mesh structure and mesher
struct Mesh
{
    Mesh(const CFD::InputData &);
    intType nCells_x, nCells_y, nCells_z;
    array1D cellCenters_x, cellCenters_y, cellCenters_z;
    array2D cellFaceAreas_x, cellFaceAreas_y, cellFaceAreas_z; // Index by right hand rule
};


// Storage of cell face velocities
struct FaceVelocities
{
    // Faces are staggered in the negative direction:
    //   cellFaceVelocity_x(i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceVelocity_y(i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceVelocity_z(i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    FaceVelocities(const intType n_x, const intType n_y, const intType n_z);
    array3D cellFaceVelocities_x, cellFaceVelocities_y, cellFaceVelocities_z;
};


// Velocity and pressure solution arrays
struct SolutionFields
{
    SolutionFields(const intType n_x, const intType n_y, const intType n_z);
    array3D U, V, W, P;
};


}   // end namespace CFD


#endif // FV_STRUCTURES