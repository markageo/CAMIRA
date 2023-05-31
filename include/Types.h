#ifndef TYPES
#define TYPES

#include "Eigen/unsupported/CXX11/Tensor"
#include <type_traits>
#include <vector>
#include <memory>
#include <utility>
#include <algorithm>

namespace CFD 
{

// ----------------------------------------------------- Type Aliases ----------------------------------------------------- //

#ifdef DOUBLE_PRECISION
    using floatType = double;
#else
    using floatType = float;
#endif

using intType = Eigen::Index;
using array0D = Eigen::Tensor<floatType, 0>;    // Colum major
using array1D = Eigen::Tensor<floatType, 1>;    // Column major
using array2D = Eigen::Tensor<floatType, 2>;    // Column major
using array3D = Eigen::Tensor<floatType, 3>;    // Column major
using indexVector3 = Eigen::Array<intType, 3, 1>;
using indexVector2 = Eigen::Array<intType, 2, 1>;
using floatVector3 = Eigen::Array<floatType, 3, 1>;
using floatVector2 = Eigen::Array<floatType, 2, 1>;





// -------------------------------------------------- Ghost Cell Indexing ------------------------------------------------- //

// Number of ghost cells in solution field
constexpr intType nGhost = 2;


// Eigen::Tensor ghost cell indexing - array indexing
template<std::size_t dim>
inline Eigen::array<Eigen::Index, dim> G( Eigen::array<Eigen::Index, dim> idx )
{
    for (auto it = idx.begin(); it != idx.end(); it++) {
        *it += nGhost;                               
    }
    return idx;
}


// Eigen::Tensor ghost cell indexing - list indexing
template<class ...Args> 
requires( sizeof...(Args) > 1 )
inline Eigen::array< Eigen::Index, sizeof...(Args) > G( Args... args )
{  
    Eigen::array< Eigen::Index, sizeof...(Args) > idx( {args...} );
    return G( idx ); 
}

// Special case for 1D, should just return a number
inline Eigen::Index G( Eigen::Index idx )
{
    return idx + nGhost;
}



// --------------------------------------------------------- Enums -------------------------------------------------------- //
// Place inside structs to avoid name conflicts with "count"

struct Axis
{
    enum ENUMDATA
    {
        X = 0,
        Y = 1,
        Z = 2,
    };
    const static int count =  3;
};


struct Fields
{
    enum ENUMDATA
    {
        U, // x velocity
        V, // y velocity
        W, // z velocity
        P, // Pressure
    };
    const static int count = 4;
};


struct BoundaryConditions 
{
    enum ENUMDATA 
    {
        zeroGradient,
        uniform,
        extrapolated,
    };
    const static int count = 3;
};


struct BoundaryPatches 
{   
    enum ENUMDATA
    {
        xPositive,
        xNegative,
        yPositive,
        yNegative,
        zPositive,
        zNegative,
    };
    const static int count = 6; 
};


struct TransportCoefficients
{
    enum ENUMDATA
    {
        tt, // (i  , j  , k+2)
        nn, // (i  , j+2, k  )
        ee, // (i+2, j  , k  )
        t,  // (i  , j  , k+1)
        n,  // (i  , j+1, k  )
        e,  // (i+1, j  , k  )
        p,  // (i  , j  , k  )
        w,  // (i-1, j  , k  )
        s,  // (i  , j-1, k  )
        b,  // (i  , j  , k-1)
        ww, // (i-2, j  , k  )
        ss, // (i  , j-2, k  )
        bb, // (i  , j  , k-2)
    };
    const static int count = 13;
};


// For looping through enum and applying lambda to each element
template<typename enumStruct, typename L>
void EnumFor(L&& f)
{
    static_assert(std::is_same<enumStruct, Axis                 >::value ||
                  std::is_same<enumStruct, Fields               >::value ||
                  std::is_same<enumStruct, BoundaryConditions   >::value ||
                  std::is_same<enumStruct, BoundaryPatches      >::value ||
                  std::is_same<enumStruct, TransportCoefficients>::value);

    typename enumStruct::ENUMDATA enumName;
    for (int i = 0; i != enumStruct::count; i++) {
        enumName = static_cast<typename enumStruct::ENUMDATA>(i);
        f(enumName);
    }
}





// ----------------------------------------------------- Enum Lookups ----------------------------------------------------- //


// Get index offset from TransportCoeffienct
constexpr std::array<intType, TransportCoefficients::count> CoeffIndex = {2,     // tt
                                                                          2,     // nn    
                                                                          2,     // ee
                                                                          1,     // t
                                                                          1,     // n
                                                                          1,     // e
                                                                          0,     // p
                                                                          -1,    // w
                                                                          -1,    // s
                                                                          -1,    // b
                                                                          -2,    // ww
                                                                          -2,    // ss
                                                                          -2};   // bb


// Get BoundaryPatches from Axis
constexpr std::array<BoundaryPatches::ENUMDATA, 3> PositivePatch = {BoundaryPatches::ENUMDATA::xPositive,
                                                                    BoundaryPatches::ENUMDATA::yPositive,
                                                                    BoundaryPatches::ENUMDATA::zPositive};

constexpr std::array<BoundaryPatches::ENUMDATA, 3> NegativePatch{BoundaryPatches::ENUMDATA::xNegative,
                                                                 BoundaryPatches::ENUMDATA::yNegative,
                                                                 BoundaryPatches::ENUMDATA::zNegative};

// Get TransportCoefficient from Axis
constexpr std::array<TransportCoefficients::ENUMDATA, 3> PositiveCoeff{TransportCoefficients::ENUMDATA::e,
                                                                       TransportCoefficients::ENUMDATA::n,
                                                                       TransportCoefficients::ENUMDATA::t};


constexpr std::array<TransportCoefficients::ENUMDATA, 3> NegativeCoeff{TransportCoefficients::ENUMDATA::w,
                                                                       TransportCoefficients::ENUMDATA::s,
                                                                       TransportCoefficients::ENUMDATA::b};


// Get Axis from BoundaryPatches
constexpr std::array<Axis::ENUMDATA, 6> BoundaryPatchAxis{Axis::ENUMDATA::X,    // xPositive
                                                          Axis::ENUMDATA::X,    // xNegative
                                                          Axis::ENUMDATA::Y,    // yPositive
                                                          Axis::ENUMDATA::Y,    // yNegative
                                                          Axis::ENUMDATA::Z,    // zPositive
                                                          Axis::ENUMDATA::Z};   // zNegative





// ---------------------------------------------------- Solver Parameters -------------------------------------------------- //

// Solver settings
enum class LinearSolvers {
    SUGS
};

enum class PlaneSolvers {
    SUGS
};

enum class LineSolvers {
    SUGS
};

enum class Linearisation {
    Picard, Newton
};





// ------------------------------------------------------- Containers ----------------------------------------------------- //

// Namespace for internal implementation.
namespace Internal
{

// Contain the types needed when constructing multidimensional arrays. 1D arrays are not constructed using arrays for
// their indices.
template<class B>
struct dimTypes
{ 
    using dimsArray = intType;
    using dimsArrayInternal = intType;

    static dimsArrayInternal ConvertDimsArrayInternal( dimsArray dims ) 
    { return dims; }
};

// Specialisation for multidimensional arrays, which are constructed using arrays for thier dimensions
template<class B>
requires( std::is_same< B, array2D >::value || std::is_same< B, array3D >::value )
struct dimTypes<B>
{ 
    using dimsArray = Eigen::Array<intType, B::NumDimensions, 1>;
    using dimsArrayInternal = Eigen::array<intType, B::NumDimensions>;

    // For converting Eigen::Array to Eigen::array. 
    // Only Eigen::array<Eigen::Index, ...> can be used to construct tensors.
    static dimsArrayInternal ConvertDimsArrayInternal( const dimsArray &dims ) 
    {
        dimsArrayInternal dimsInternal;
        for ( Eigen::Index i = 0; i != B::NumDimensions; i++ ) {
            dimsInternal[ static_cast<size_t>( i ) ] = dims(i);
        }
        return dimsInternal;
    }
    
};

}   //  end namespace Internal


// Wrapper for std::vector that can only be indexed using enums
template <typename enumStruct, typename T>
class EnumVector
{
    static_assert(std::is_same<enumStruct, CFD::Axis                 >::value ||
                  std::is_same<enumStruct, CFD::Fields               >::value ||
                  std::is_same<enumStruct, CFD::BoundaryConditions   >::value ||
                  std::is_same<enumStruct, CFD::BoundaryPatches      >::value ||
                  std::is_same<enumStruct, CFD::TransportCoefficients>::value,
                  "Template parameter must be struct containing ENUMDATA type.");

    typedef typename enumStruct::ENUMDATA ENUMDATA;

    static constexpr bool isArray = std::is_same< T, CFD::array3D >::value || 
                                    std::is_same< T, CFD::array2D >::value || 
                                    std::is_same< T, CFD::array1D >::value;


    // Construction for 1D arrays is handeled a bit differently
    static constexpr bool isArrayND = std::is_same< T, CFD::array2D >::value ||
                                      std::is_same< T, CFD::array3D >::value;

    using dimsArray =  typename Internal::dimTypes<T>::dimsArray;
    using dimsArrayInternal =  typename Internal::dimTypes<T>::dimsArrayInternal;

    public:

        // Constructors for general types
        EnumVector() : m_dataVector( enumStruct::count ) {};
        EnumVector( const T &data ) : m_dataVector( enumStruct::count, data )  { };
        EnumVector( const std::array<T, enumStruct::count> &arr ) : m_dataVector( arr.begin(), arr.end() ) { };

        // Special constructors for array objects, all having same dimenions
        EnumVector(const std::vector< ENUMDATA > &coeffs, const dimsArray &dims) 
        requires( isArray ) :
            m_dataVector(enumStruct::count)
        {
            dimsArrayInternal dimsInternal = Internal::dimTypes<T>::ConvertDimsArrayInternal( dims );
            for (const auto &index : coeffs) {
                m_dataVector[index] = T( dimsInternal ).setZero();
            }
        }

        // Special constructors for array objects, can have different dimensions
        EnumVector( const std::vector< std::pair< ENUMDATA, dimsArray > > &arraySpec)  
        requires ( isArray ) :
            m_dataVector(enumStruct::count)
        {
            dimsArrayInternal dimsInternal;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dimsInternal = Internal::dimTypes<T>::ConvertDimsArrayInternal( arraySpec[i].second );
                m_dataVector[ arraySpec[i].first ] = T( dimsInternal ).setZero();
            }
        }

        // Strong type indexing
        T &operator[](const typename enumStruct::ENUMDATA idx)
        { return m_dataVector[idx]; }

        const T &operator[](const typename enumStruct::ENUMDATA idx) const 
        { return m_dataVector[idx]; }

        // Get underlying data vector
        std::vector<T> &get()
        { return m_dataVector; }

    private:
        std::vector<T> m_dataVector;

};





// Allocate arrays using enums. Arrays are initialised to zero.
template <typename enumStruct, typename arrayType>
class ArrayAllocator
{
    static_assert(std::is_same<enumStruct, CFD::Axis                 >::value ||
                  std::is_same<enumStruct, CFD::Fields               >::value ||
                  std::is_same<enumStruct, CFD::BoundaryConditions   >::value ||
                  std::is_same<enumStruct, CFD::BoundaryPatches      >::value ||
                  std::is_same<enumStruct, CFD::TransportCoefficients>::value,
                  "Template parameter must be struct containing ENUMDATA type.");

    static_assert(std::is_same< arrayType, CFD::array1D >::value ||
                  std::is_same< arrayType, CFD::array2D >::value ||
                  std::is_same< arrayType, CFD::array3D >::value,
                  "Array type invalid.");

    typedef typename enumStruct::ENUMDATA ENUMDATA;

    // Construction for 1D arrays is handeled a bit differently
    static constexpr bool isArrayND = std::is_same< arrayType, CFD::array2D >::value ||
                                      std::is_same< arrayType, CFD::array3D >::value;

    using dimsArray =  typename Internal::dimTypes<arrayType>::dimsArray;
    using dimsArrayInternal = typename Internal::dimTypes<arrayType>::dimsArrayInternal;

    public:

        // Constructor, All arrays have same dimensions 
        ArrayAllocator(const std::vector< ENUMDATA > &coeffs, const dimsArray &dims ) :
            m_coeffPointers(enumStruct::count)
        {
            dimsArrayInternal dimsInternal = Internal::dimTypes<arrayType>::ConvertDimsArrayInternal( dims );

            for (const auto &index : coeffs) {
                m_coeffPointers[index] = std::make_unique<arrayType>( arrayType( dimsInternal ).setZero() );
            }
        }


        // Constructor, arrays can have different dimensions
        ArrayAllocator( const std::vector< std::pair< ENUMDATA, dimsArray > > &arraySpec) : 
            m_coeffPointers(enumStruct::count)
        {
            dimsArrayInternal dimsInternal;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dimsInternal =  Internal::dimTypes<arrayType>::ConvertDimsArrayInternal( arraySpec[i].second );
                m_coeffPointers[ arraySpec[i].first ] = std::make_unique<arrayType>( arrayType( dimsInternal ).setZero() );
            }
        }


        // Copy Constructor
        ArrayAllocator(const ArrayAllocator &that) :
            m_coeffPointers(enumStruct::count)
        {
            // Allocate a new array object only if it was allocated in the original
            for (size_t i = 0; i != m_coeffPointers.size(); i++) {
                if (that.m_coeffPointers[i]) {
                    m_coeffPointers[i] = std::make_unique<arrayType>( *that.m_coeffPointers[i] );
                }
            }
        }


        // Copy assignment
        ArrayAllocator &operator=(ArrayAllocator that)
        {   
            std::swap( this->m_coeffPointers, that.m_coeffPointers );
            return *this;
        }


        // Move constructor
        ArrayAllocator(ArrayAllocator&& that) noexcept :
            m_coeffPointers( std::move( that.m_coeffPointers ) )
        {}


        // Move assignment
        // No need due to copy-swap 


        // Indexing operators
        arrayType &operator[](const typename enumStruct::ENUMDATA idx)
        { return *m_coeffPointers[idx]; }

        const arrayType &operator[](const typename enumStruct::ENUMDATA idx) const 
        { return *m_coeffPointers[idx]; }


        // Return the pointer to the enumed object
        std::unique_ptr<arrayType> &get(const typename enumStruct::ENUMDATA idx)
        { return m_coeffPointers[idx]; }

        const std::unique_ptr<arrayType> &get(const typename enumStruct::ENUMDATA idx) const
        { return m_coeffPointers[idx]; }



    private:
        std::vector< std::unique_ptr<arrayType> > m_coeffPointers;

};


}   // end namespace CFD

#endif // TYPES