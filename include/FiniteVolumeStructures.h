#ifndef FV_STRUCTURES
#define FV_STRUCTURES

#include "SimulationParameters.h"
#include "InputProcessing.h"
#include "Tensor"

namespace CFD
{

    namespace TransportCoefficients
    {
        enum ENUMDATA
        {
            p,  // (i  , j  , k  )
            n,  // (i  , j+1, k  )
            e,  // (i+1, j  , k  )
            s,  // (i  , j-1, k  )
            w,  // (i-1, j  , k  )
            t,  // (i  , j  , k+1)
            b,  // (i  , j  , k-1)
            nn, // (i  , j+2, k  )
            ee, // (i+2, j  , k  )
            ss, // (i  , j-2, k  )
            ww, // (i-2, j  , k  )
            tt, // (i  , j  , k+2)
            bb, // (i  , j  , k-2)
            count
        };
    }

    namespace Fields
    {
        enum ENUMDATA
        {
            U, // x velocity
            V, // y velocity
            W, // z velocity
            P, // Pressure
            count
        };
    }

    // Allocate 3D arrays using enums
    template <typename arrayEnum = TransportCoefficients::ENUMDATA>
    class ArrayAllocator
    {
        static_assert(std::is_enum<arrayEnum>::value, "Template parameter must be enum type.");

    public:
        ArrayAllocator(const intType n_x, const intType n_y, const intType n_z, const std::vector<arrayEnum> coeffs) : 
            coeffsUsed(coeffs), coeffPointers(arrayEnum::count)
        {
            for (const auto &index : coeffsUsed)
            {
                coeffPointers[index] = new CFD::array3D(n_x, n_y, n_z);
            }
        }

        ~ArrayAllocator()
        {
            for (const auto &index : coeffsUsed)
            {
                delete coeffPointers[index];
            }
        }

        CFD::array3D &operator[](arrayEnum idx)
        {
            return *coeffPointers[idx];
        }

    private:
        std::vector<arrayEnum> coeffsUsed;
        std::vector<CFD::array3D *> coeffPointers;
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