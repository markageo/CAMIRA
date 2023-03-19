#ifndef FV_STRUCTURES
#define FV_STRUCTURES

#include "SimulationParameters.h"
#include "InputProcessing.h"
#include "Tensor"
#include <memory>

namespace CFD
{

// Allocate 3D arrays using enums
template <typename arrayEnum>
class ArrayAllocator
{
    static_assert(std::is_enum<arrayEnum>::value, "Template parameter must be enum type.");

    public:
        ArrayAllocator(const intType n_x, const intType n_y, const intType n_z, const std::vector<arrayEnum> coeffs) : 
            coeffsUsed(coeffs), coeffPointers(arrayEnum::count)
        {
            for (const auto &index : coeffsUsed)
            {
                coeffPointers[index] = std::make_unique<CFD::array3D>( CFD::array3D(n_x, n_y, n_z) );
            }
        }

        CFD::array3D &operator[](arrayEnum idx)
        {
            return *coeffPointers[idx];
        }

    private:
        std::vector<arrayEnum> coeffsUsed;
        std::vector< std::unique_ptr<CFD::array3D> > coeffPointers;
};


// Recitlinear mesh structure and mesher
struct Mesh
{
    Mesh(const CFD::InputData &);
    intType nCells_x, nCells_y, nCells_z;
    array1D cellCenters_x, cellCenters_y, cellCenters_z;
    array2D cellFaceAreas_x, cellFaceAreas_y, cellFaceAreas_z; // Index by right hand rule
};


// Storage of cell face normal velocities
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

} // end namespace CFD

#endif // FV_STRUCTURES