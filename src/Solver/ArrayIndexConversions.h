#ifndef INDEX_CONVERSIONS
#define INDEX_CONVERSIONS

#include "../Core/Types.h"

namespace CFD
{

// Convert list of subscript array indices into linear 
template< int ndims,
          Eigen::StorageOptions layout = Eigen::ColMajor, 
          class ...Args > 
__attribute__((always_inline)) 
inline intType Sub2Ind( const Eigen::Array<intType, ndims, 1> &dims, Args... args )
{
    static_assert( ndims == sizeof...(Args) );
    static_assert( layout == Eigen::ColMajor || 
                   layout == Eigen::RowMajor  );

    intType idx = 0;
    if constexpr ( layout == Eigen::ColMajor ) {

        // idx = 0;
        // idx = idx * m_dims[ndims - 1] + args[ndims - 1]
        // idx = idx * m_dims[ndims - 2] + args[ndims - 2]
        // idx = idx * m_dims[ndims - 3] + args[ndims - 3]
        //              :
        //              :
        // idx = idx * m_dims[0] + args[0]
        intType i = ndims-1;
        intType dummy;
        ( dummy = ... = (idx = idx*dims[i--] + args ) );  // Trick to evaluate fold expression in reverse

    } else if constexpr ( layout == Eigen::RowMajor ) {

        // idx = 0;
        // idx = idx * m_dims[0] + args[0]
        // idx = idx * m_dims[1] + args[1]
        // idx = idx * m_dims[2] + args[2]
        //              :
        //              :
        // idx = idx * m_dims[ndims - 1] + args[ndims - 1]
        intType i = 0;
        ( ... , (idx = idx*dims[i++] + args ));
    }

    return idx;

}



// Convert linear index into list of subscripts
// TODO: Optimise this so that the loops are unrolled at compile time e.g. use fold expressions like above
template< int ndims,
          Eigen::StorageOptions layout = Eigen::ColMajor > 
__attribute__((always_inline)) 
Eigen::array<Eigen::Index, ndims> Ind2Sub( const Eigen::Array<intType, ndims, 1> &dims, intType idx )
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


}   // end namespace CFD

#endif // INDEX_CONVERSIONS