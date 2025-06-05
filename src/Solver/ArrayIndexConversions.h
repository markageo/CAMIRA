#ifndef INDEX_CONVERSIONS
#define INDEX_CONVERSIONS

#include "../Core/Types.h"
#include <utility>

namespace CAMIRA
{


// Convert list of subscript array indices into linear index
template< int ndims,
          Eigen::StorageOptions layout = Eigen::ColMajor, 
          class ...Args > 
__attribute__((always_inline)) 
inline constexpr intType Sub2Ind( const Eigen::Array<intType, ndims, 1> &dims, Args... args )
{
    static_assert( ndims == sizeof...(Args) );
    static_assert( layout == Eigen::ColMajor || 
                   layout == Eigen::RowMajor  );

    std::array<std::common_type_t<Args...>, ndims> indices = {args...};

    intType idx = 0;
    if constexpr ( layout == Eigen::ColMajor ) {

        for ( intType i = 0; i != ndims; i++ ) {

            intType stride = 1;
            for ( intType j = 0; j != i; j++ ) {
                stride *= dims(j);
            }

            idx += stride * indices[i];
        }

    } else if constexpr ( layout == Eigen::RowMajor ) {

        for ( intType i = 0; i != ndims; i++ ) {

            intType stride = 1;
            for ( intType j = ndims; j > i; j-- ) {
                stride *= dims(j);
            }

            idx += stride * indices[i];
        }

    }

    return idx;

}



// Convert linear index into list of subscripts
template< int ndims,
          Eigen::StorageOptions layout = Eigen::ColMajor > 
__attribute__((always_inline)) 
inline constexpr Eigen::array<Eigen::Index, ndims> Ind2Sub( const Eigen::Array<intType, ndims, 1> &dims, intType idx )
{
    static_assert( layout == Eigen::ColMajor || 
                   layout == Eigen::RowMajor  );
    static_assert( ndims > 1 );

    std::array<intType, ndims> indices;
    if constexpr ( layout == Eigen::ColMajor ) {

        for ( intType i = ndims-1; i != -1; i-- ) {

            intType stride = 1;
            for ( intType j = 0; j < i; j++ ) {
                stride *= dims(j);
            }

            indices[i] = idx / stride;
            idx %= stride;
        }

    } else if constexpr ( layout == Eigen::RowMajor ) {

        for ( intType i = 0; i != ndims; i++ ) {

            intType stride = 1;
            for ( intType j = ndims-1; j != i; j-- ) {
                stride *= dims(j);
            }

            indices[i] = idx / stride;
            idx %= stride;
        }

    }

    return indices;

}


}   // end namespace CAMIRA

#endif // INDEX_CONVERSIONS