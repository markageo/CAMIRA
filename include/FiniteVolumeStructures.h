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
template <typename arrayEnum, typename arrayType = CFD::array3D>
class ArrayAllocator
{
    static_assert(std::is_enum<arrayEnum>::value, "Template parameter must be enum type.");

    public:

        // ----------------------------- Constructor: All arrays have same dimensions ----------------------------- //

        // 3D
        ArrayAllocator(const std::vector<arrayEnum> &coeffs, const indexVector3 &dims) requires ( std::is_same< arrayType, CFD::array3D >::value ) : 
            coeffPointers(arrayEnum::count)
        {
            for (const auto &index : coeffs) {
                coeffPointers[index] = std::make_unique<CFD::array3D>( CFD::array3D( dims(0), dims(1), dims(2) ) );
            }
        }

        // 2D
        ArrayAllocator(const std::vector<arrayEnum> &coeffs, const indexVector2 &dims) requires ( std::is_same< arrayType, CFD::array2D >::value ) : 
            coeffPointers(arrayEnum::count)
        {
            for (const auto &index : coeffs) {
                coeffPointers[index] = std::make_unique<CFD::array2D>( CFD::array2D( dims(0), dims(1) ) );
            }
        }

        // 1D
        ArrayAllocator(const std::vector<arrayEnum> &coeffs, const intType &dim) requires ( std::is_same< arrayType, CFD::array1D >::value ) : 
            coeffPointers(arrayEnum::count)
        {
            for (const auto &index : coeffs) {
                coeffPointers[index] = std::make_unique<CFD::array1D>( CFD::array1D( dim ) );
            }
        }


        // --------------------------- Constructor: Arrays can have different dimensions --------------------------- //

        // 3D
        ArrayAllocator( const std::vector< std::pair< arrayEnum, CFD::indexVector3 > > &arraySpec) requires ( std::is_same< arrayType, CFD::array3D >::value ) : 
            coeffPointers(arrayEnum::count)
        {
            CFD::indexVector3 dims;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dims = arraySpec[i].second;
                coeffPointers[ arraySpec[i].first ] = std::make_unique<CFD::array3D>( CFD::array3D( dims(0), dims(1), dims(2) ) );
            }
        }


        // 2D
        ArrayAllocator( const std::vector< std::pair< arrayEnum, CFD::indexVector2 > > &arraySpec) requires ( std::is_same< arrayType, CFD::array2D >::value ) : 
            coeffPointers(arrayEnum::count)
        {
            CFD::indexVector2 dims;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dims = arraySpec[i].second;
                coeffPointers[ arraySpec[i].first ] = std::make_unique<CFD::array2D>( CFD::array2D( dims(0), dims(1) ) );
            }
        }


        // 1D
        ArrayAllocator( const std::vector< std::pair< arrayEnum, CFD::intType > > &arraySpec) requires ( std::is_same< arrayType, CFD::array1D >::value ) : 
            coeffPointers(arrayEnum::count)
        {
            CFD::intType dim;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dim = arraySpec[i].second;
                coeffPointers[ arraySpec[i].first ] = std::make_unique<CFD::array1D>( CFD::array1D( dim ) );
            }
        }



        // ----------------------------------- Array reference return operators ----------------------------------- //

        arrayType &operator[](arrayEnum idx)
        {
            return *coeffPointers[idx];
        }

        arrayType &operator[](arrayEnum idx) const
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
    indexVector3 nCells;
    array1D cellCenters_x, cellCenters_y, cellCenters_z;
    array1D cellFaces_x, cellFaces_y, cellFaces_z;
    array1D cellLengths_x, cellLengths_y, cellLengths_z;
    array1D interpFactors_x, interpFactors_y, interpFactors_z; // faceValue(i) = (1 - interpFactor(i))*cellValue(i-1) + interpFactor(i)*cellValue(i)
    array2D cellFaceAreas_x, cellFaceAreas_y, cellFaceAreas_z; // Index by right hand rule
};


} // end namespace CFD

#endif // FV_STRUCTURES