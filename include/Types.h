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

using floatType = double;
using intType = int;
using iterType = int;
using array1D = Eigen::Tensor<floatType, 1>;    // Column major
using array2D = Eigen::Tensor<floatType, 2>;    // Column major
using array3D = Eigen::Tensor<floatType, 3>;    // Column major
using indexVector3 = Eigen::Array<intType, 3, 1>;
using indexVector2 = Eigen::Array<intType, 2, 1>;
using floatVector3 = Eigen::Array<floatType, 3, 1>;
using floatVector2 = Eigen::Array<floatType, 2, 1>;


// Number of ghost cells in solution field
constexpr intType nGhost = 2;

// Convert regular indexing to ghost cell indexing, to use with Eigen::Tensor
inline Eigen::array<Eigen::Index, 1> G(const Eigen::Index i) 
    { return {i + nGhost}; };

inline Eigen::array<Eigen::Index, 2> G(const Eigen::Index i, const Eigen::Index j) 
    { return {i + nGhost, j + nGhost}; };

inline Eigen::array<Eigen::Index, 3> G(const Eigen::Index i, const Eigen::Index j, const Eigen::Index k) 
    { return {i + nGhost, j + nGhost, k + nGhost}; };



// Enums to be used as indices for containers
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


// Lookup arrays for coefficients and patches coerresponding to each axis
// *** Should maybe move these to the FVCoefficients.cpp implementation file if they are not used anywhere else ***
constexpr std::array<BoundaryPatches::ENUMDATA, 3> positivePatches = {BoundaryPatches::ENUMDATA::xPositive,
                                                                      BoundaryPatches::ENUMDATA::yPositive,
                                                                      BoundaryPatches::ENUMDATA::zPositive};

constexpr std::array<BoundaryPatches::ENUMDATA, 3> negativePatches{BoundaryPatches::ENUMDATA::xNegative,
                                                                   BoundaryPatches::ENUMDATA::yNegative,
                                                                   BoundaryPatches::ENUMDATA::zNegative};

constexpr std::array<TransportCoefficients::ENUMDATA, 3> eastCoefficients{TransportCoefficients::ENUMDATA::e,
                                                                          TransportCoefficients::ENUMDATA::n,
                                                                          TransportCoefficients::ENUMDATA::t};


constexpr std::array<TransportCoefficients::ENUMDATA, 3> westCoefficients{TransportCoefficients::ENUMDATA::w,
                                                                          TransportCoefficients::ENUMDATA::s,
                                                                          TransportCoefficients::ENUMDATA::b};


// Lookup array for determining Axis based on BoundaryPatch
constexpr std::array<Axis::ENUMDATA, 6> boundaryPatchAxis{Axis::ENUMDATA::X,    // xPositive
                                                          Axis::ENUMDATA::X,    // xNegative
                                                          Axis::ENUMDATA::Y,    // yPositive
                                                          Axis::ENUMDATA::Y,    // yNegative
                                                          Axis::ENUMDATA::Z,    // zPositive
                                                          Axis::ENUMDATA::Z};   // zNegative



// Wrapper for std::array that can only be indexed using enums
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

    public:

        // Constructors for general types
        EnumVector() : m_dataVector( enumStruct::count ) {};
        EnumVector( const T &data ) : m_dataVector( enumStruct::count, data )  { };
        EnumVector( const std::array<T, enumStruct::count> &arr ) : m_dataVector( arr.begin(), arr.end() ) { };

        // Special constructors for array3D, array2D, and array1D objects, all having same dimenions
        EnumVector(const std::vector< ENUMDATA > &coeffs, const indexVector3 &dims) // 3D
        requires ( std::is_same< T, CFD::array3D >::value ) :
            m_dataVector(enumStruct::count)
        {
            for (const auto &index : coeffs) {
                m_dataVector[index] = CFD::array3D( dims(0), dims(1), dims(2) ).setZero();
            }
        }

        EnumVector(const std::vector< ENUMDATA > &coeffs, const indexVector2 &dims) // 2D
        requires ( std::is_same< T, CFD::array2D >::value ) :
            m_dataVector(enumStruct::count)
        {
            for (const auto &index : coeffs) {
                m_dataVector[index] = CFD::array2D( dims(0), dims(1) ).setZero();
            }
        }

        EnumVector(const std::vector< ENUMDATA > &coeffs, const intType &dim) // 1D
        requires ( std::is_same< T, CFD::array1D >::value ) :
            m_dataVector(enumStruct::count)
        {
            for (const auto &index : coeffs) {
                m_dataVector[index] = CFD::array1D( dim ).setZero();
            }
        }


        // Special constructors for array3D, array2D, and array1D objects, can have different dimensions
        EnumVector( const std::vector< std::pair< ENUMDATA, CFD::indexVector3 > > &arraySpec)   // 3D
        requires ( std::is_same< T, CFD::array3D >::value ) : 
            m_dataVector(enumStruct::count)
        {
            CFD::indexVector3 dims;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dims = arraySpec[i].second;
                m_dataVector[ arraySpec[i].first ] = CFD::array3D( dims(0), dims(1), dims(2) ).setZero();
            }
        }

        EnumVector( const std::vector< std::pair< ENUMDATA, CFD::indexVector2 > > &arraySpec)   // 2D
        requires ( std::is_same< T, CFD::array2D >::value ) : 
            m_dataVector(enumStruct::count)
        {
            CFD::indexVector2 dims;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dims = arraySpec[i].second;
                m_dataVector[ arraySpec[i].first ] = CFD::array2D( dims(0), dims(1) ).setZero();
            }
        }

        EnumVector( const std::vector< std::pair< ENUMDATA, CFD::intType > > &arraySpec)   // 1D
        requires ( std::is_same< T, CFD::array1D >::value ) : 
            m_dataVector(enumStruct::count)
        {
            CFD::intType dim;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dim = arraySpec[i].second;
                m_dataVector[ arraySpec[i].first ] = CFD::array1D( dim ).setZero();
            }
        }


        // Strong type indexing
        T &operator[](const enumStruct::ENUMDATA idx)
        {
            return m_dataVector[idx];
        }

        const T &operator[](const enumStruct::ENUMDATA idx) const 
        {
            return m_dataVector[idx];
        }

        // Get underlying data vector
        std::vector<T> &get()
        {
            return m_dataVector;
        }

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

    public:

        // Constructor, all arrays have the same dimensions
        // 3D
        ArrayAllocator(const std::vector< ENUMDATA > &coeffs, const indexVector3 &dims) 
        requires ( std::is_same< arrayType, CFD::array3D >::value ) : 
            m_coeffPointers(enumStruct::count)
        {
            for (const auto &index : coeffs) {
                m_coeffPointers[index] = std::make_unique<CFD::array3D>( CFD::array3D( dims(0), dims(1), dims(2) ).setZero() );
            }
        }

        // 2D
        ArrayAllocator(const std::vector< ENUMDATA > &coeffs, const indexVector2 &dims) 
        requires ( std::is_same< arrayType, CFD::array2D >::value ) : 
            m_coeffPointers(enumStruct::count)
        {
            for (const auto &index : coeffs) {
                m_coeffPointers[index] = std::make_unique<CFD::array2D>( CFD::array2D( dims(0), dims(1) ).setZero() );
            }
        }

        // 1D
        ArrayAllocator(const std::vector< ENUMDATA > &coeffs, const intType &dim) 
        requires ( std::is_same< arrayType, CFD::array1D >::value ) : 
            m_coeffPointers(enumStruct::count)
        {
            for (const auto &index : coeffs) {
                m_coeffPointers[index] = std::make_unique<CFD::array1D>( CFD::array1D( dim ).setZero() );
            }
        }


        // Constructor, array can have difference dimensions
        // 3D
        ArrayAllocator( const std::vector< std::pair< ENUMDATA, CFD::indexVector3 > > &arraySpec) 
        requires ( std::is_same< arrayType, CFD::array3D >::value ) : 
            m_coeffPointers(enumStruct::count)
        {
            CFD::indexVector3 dims;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dims = arraySpec[i].second;
                m_coeffPointers[ arraySpec[i].first ] = std::make_unique<CFD::array3D>( CFD::array3D( dims(0), dims(1), dims(2) ).setZero() );
            }
        }

        // 2D
        ArrayAllocator( const std::vector< std::pair< ENUMDATA, CFD::indexVector2 > > &arraySpec) 
        requires ( std::is_same< arrayType, CFD::array2D >::value ) : 
            m_coeffPointers(enumStruct::count)
        {
            CFD::indexVector2 dims;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dims = arraySpec[i].second;
                m_coeffPointers[ arraySpec[i].first ] = std::make_unique<CFD::array2D>( CFD::array2D( dims(0), dims(1) ).setZero() );
            }
        }

        // 1D
        ArrayAllocator( const std::vector< std::pair< ENUMDATA, CFD::intType > > &arraySpec) 
        requires ( std::is_same< arrayType, CFD::array1D >::value ) : 
            m_coeffPointers(enumStruct::count)
        {
            CFD::intType dim;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dim = arraySpec[i].second;
                m_coeffPointers[ arraySpec[i].first ] = std::make_unique<CFD::array1D>( CFD::array1D( dim ).setZero() );
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
        ArrayAllocator &operator=(ArrayAllocator &&that) noexcept
        {
            return ArrayAllocator( std::move( that ) );
        }    

        // Destructor
        ~ArrayAllocator() = default;


        // Indexing operators
        arrayType &operator[](const enumStruct::ENUMDATA idx)
        {
            return *m_coeffPointers[idx];
        }

        const arrayType &operator[](const enumStruct::ENUMDATA idx) const 
        {
            return *m_coeffPointers[idx];
        }


        // Return the pointer to the enumed object
        std::unique_ptr<arrayType> &get(const enumStruct::ENUMDATA idx)
        {
            return m_coeffPointers[idx];
        }

        const std::unique_ptr<arrayType> &get(const enumStruct::ENUMDATA idx) const
        {
            return m_coeffPointers[idx];
        }



    private:
        std::vector< std::unique_ptr<arrayType> > m_coeffPointers;

};


}   // end namespace CFD

#endif // TYPES