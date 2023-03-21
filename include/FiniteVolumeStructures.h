#ifndef FV_STRUCTURES
#define FV_STRUCTURES

#include "Types.h"
#include "InputProcessing.h"
#include "Tensor"
#include <memory>
#include <utility>

namespace CFD
{

// Allocate 3D arrays using enums
template <typename arrayEnum>
class ArrayAllocator
{
    static_assert(std::is_enum<arrayEnum>::value, "Template parameter must be enum type.");

    public:

        // All arrays have same dimensions
        ArrayAllocator(const std::vector<arrayEnum> coeffs, const indexVector &dims) : 
            coeffPointers(arrayEnum::count)
        {
            for (const auto &index : coeffs) {
                coeffPointers[index] = std::make_unique<CFD::array3D>( CFD::array3D( dims(0), dims(1), dims(2) ) );
            }
        }

        // Arrays can have different dimensions
        ArrayAllocator( const std::vector< std::pair< arrayEnum, CFD::indexVector > > &arraySpec) : 
            coeffPointers(arrayEnum::count)
        {
            CFD::indexVector dims;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dims = arraySpec[i].second;
                coeffPointers[ arraySpec[i].first ] = std::make_unique<CFD::array3D>( CFD::array3D( dims(0), dims(1), dims(2) ) );
            }
        }

        CFD::array3D &operator[](arrayEnum idx)
        {
            return *coeffPointers[idx];
        }

        CFD::array3D &operator[](arrayEnum idx) const
        {
            return *coeffPointers[idx];
        }

    private:
        std::vector< std::unique_ptr<CFD::array3D> > coeffPointers;
};


// Recitlinear mesh structure and mesher (on construction)
struct Mesh
{
    Mesh(const CFD::InputData &);
    indexVector nCells;
    array1D cellCenters_x, cellCenters_y, cellCenters_z;
    array1D cellFaces_x, cellFaces_y, cellFaces_z;
    array1D interpFactors_x, interpFactors_y, interpFactors_z; // faceValue(i) = (1 - interpFactor(i))*cellValue(i-1) + interpFactor(i)*cellValue(i)
    array2D cellFaceAreas_x, cellFaceAreas_y, cellFaceAreas_z; // Index by right hand rule
};


} // end namespace CFD

#endif // FV_STRUCTURES