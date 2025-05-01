#ifndef TYPES
#define TYPES

#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <type_traits>
#include <vector>
#include <utility>
#include <algorithm>


namespace CFD 
{


// ----------------------------------------------------- Type Aliases ----------------------------------------------------- //

#ifdef CFD_DOUBLE_PRECISION
    using floatType = double;
#else
    using floatType = float;
#endif

using intType = Eigen::Index;

using Tensor0D = Eigen::Tensor<floatType, 0>;    // Column major
using Tensor1D = Eigen::Tensor<floatType, 1>;    // Column major
using Tensor2D = Eigen::Tensor<floatType, 2>;    // Column major
using Tensor3D = Eigen::Tensor<floatType, 3>;    // Column major
using TensorIndex3D = Eigen::array<Eigen::Index, 3>;
using TensorIndex2D = Eigen::array<Eigen::Index, 2>;

using iArray3 = Eigen::Array<intType, 3, 1>;
using iArray2 = Eigen::Array<intType, 2, 1>;
using fArray3 = Eigen::Array<floatType, 3, 1>;
using fArray2 = Eigen::Array<floatType, 2, 1>;

using fVector3 = Eigen::Matrix<floatType, 3, 1>;
using fVector2 = Eigen::Matrix<floatType, 2, 1>;



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


struct BoundaryConditions 
{
    enum ENUMDATA 
    {
        zeroGradient,
        fixed,
        extrapolated,
        periodic
    };
    const static int count = 4;
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
inline void EnumFor( L&& f )
{
    static_assert(std::is_same<enumStruct, Axis                 >::value ||
                  std::is_same<enumStruct, BoundaryConditions   >::value ||
                  std::is_same<enumStruct, BoundaryPatches      >::value ||
                  std::is_same<enumStruct, TransportCoefficients>::value);

    typename enumStruct::ENUMDATA enumName;
    for (int i = 0; i != enumStruct::count; i++) {
        enumName = static_cast<typename enumStruct::ENUMDATA>(i);
        f(enumName);
    }
}





// ---------------------------------------------------- Solver Parameters -------------------------------------------------- //

// Solver settings
enum class Smoothers {
    nestedLineSymmetricSerial, domainSymmetricSerial, domainSymmetricParallel
};

enum class GeometryBoundaryTreatement {
    Staircase, DirectionalImmersedBoundary
};

enum class MomentumInterpolation {
    Implicit, SemiExplicit
};

enum class AdvectionSchemes {
    Upwind, Central, SOU, QUICK
};

enum class FaceInterpolationSchemes {
    WeightedLinear, Average
};

enum class TimeSchemes {
    Steady, BackwardsEuler, BackwardsThreeLevel
};

enum class MultigridCycleType {
    V, F, W
};



// ------------------------------------------------------- Containers ----------------------------------------------------- //

// Namespace for internal implementation.
namespace CFD_INTERNAL
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
    requires( std::is_same< B, Tensor2D >::value || std::is_same< B, Tensor3D >::value )
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

}   //  end namespace CFD_INTERNAL




// Wrapper for std::array that can only be indexed using enums
template <typename enumStruct, typename T>
class EnumVector
{
    static_assert(std::is_same<enumStruct, CFD::Axis                 >::value ||
                  std::is_same<enumStruct, CFD::BoundaryConditions   >::value ||
                  std::is_same<enumStruct, CFD::BoundaryPatches      >::value ||
                  std::is_same<enumStruct, CFD::TransportCoefficients>::value,
                  "Template parameter must be struct containing ENUMDATA type.");

    typedef typename enumStruct::ENUMDATA ENUMDATA;

    static constexpr bool isArray = std::is_same< T, CFD::Tensor3D >::value || 
                                    std::is_same< T, CFD::Tensor2D >::value || 
                                    std::is_same< T, CFD::Tensor1D >::value;


    // Construction for 1D arrays is handeled a bit differently
    static constexpr bool isArrayND = std::is_same< T, CFD::Tensor2D >::value ||
                                      std::is_same< T, CFD::Tensor3D >::value;

    using dimsArray =  typename CFD_INTERNAL::dimTypes<T>::dimsArray;
    using dimsArrayInternal =  typename CFD_INTERNAL::dimTypes<T>::dimsArrayInternal;

    public:

        // Constructors for general types
        EnumVector() {};
        EnumVector( const T &data ) { std::fill( m_dataVector.begin(), m_dataVector.end(), data ); };
        constexpr EnumVector( const std::array<T, enumStruct::count> &arr ) : m_dataVector( arr ) {};

        // Special constructors for array objects, all having same dimenions
        EnumVector(const std::vector< ENUMDATA > &coeffs, const dimsArray &dims) 
        requires( isArray ) 
        {
            dimsArrayInternal dimsInternal = CFD_INTERNAL::dimTypes<T>::ConvertDimsArrayInternal( dims );
            for (const auto &index : coeffs) {
                m_dataVector[index] = T( dimsInternal ).setZero();
            }
        }

        // Special constructors for array objects, can have different dimensions
        EnumVector( const std::vector< std::pair< ENUMDATA, dimsArray > > &arraySpec)  
        requires ( isArray ) 
        {
            dimsArrayInternal dimsInternal;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dimsInternal = CFD_INTERNAL::dimTypes<T>::ConvertDimsArrayInternal( arraySpec[i].second );
                m_dataVector[ arraySpec[i].first ] = T( dimsInternal ).setZero();
            }
        }

        // Strong type indexing
        constexpr T &operator[](const typename enumStruct::ENUMDATA idx)
        { return m_dataVector[idx]; }

        constexpr const T &operator[](const typename enumStruct::ENUMDATA idx) const 
        { return m_dataVector[idx]; }


        // Only allowed when the object is holding axis data since the enum values should not change
        constexpr T &operator[](const size_t idx) requires( std::is_same<enumStruct, Axis>::value )
        { return m_dataVector[idx]; }

        constexpr const T &operator[](const size_t idx) const requires( std::is_same<enumStruct, Axis>::value )
        { return m_dataVector[idx]; }


        // Get underlying data
        constexpr std::array<T, enumStruct::count> &get()
        { return m_dataVector; }

    private:
        std::array<T, enumStruct::count> m_dataVector;

};



// A general struct for holding values corresponding to different fields
template < typename dataType >
struct FieldData {
    EnumVector<Axis, dataType> U;
    dataType P;
    
    FieldData() { SetPointers(); };
    FieldData( const dataType &data ) : U( data ), P( data ) 
    { SetPointers(); };

    static constexpr intType nData = Axis::count + 1;   // This doesn't depend on datatype 

    // For assigning all data to the same value
    FieldData &operator=( dataType &data )
    {
        this->U = data;
        this->P = data;
        return *this;
    }

    dataType &operator[]( const intType idx )
    { return *m_dataPointers[ static_cast<size_t>( idx ) ]; }

    const dataType &operator[]( const intType idx ) const
    { return *m_dataPointers[ static_cast<size_t>( idx ) ]; }

    // Copy constructor
    FieldData( const FieldData &that ) : U( that.U ), P( that.P )
    { SetPointers(); }

    // Copy assignment
    FieldData &operator=( FieldData that )
    { 
        std::swap( this->U, that.U );
        std::swap( this->P, that.P );
        SetPointers();
        return *this;
    }

    // Move constructor, don't need move assignment due to copy and swap
    FieldData( FieldData&& that ) noexcept : U( std::move( that.U ) ), P( std::move( that.P ) ) 
    { SetPointers(); };


    private:
        std::array< dataType*, nData > m_dataPointers;
        
        void SetPointers()
        { m_dataPointers = { &U[Axis::X], &U[Axis::Y], &U[Axis::Z], &P }; }
};


// To allow iterating through values of a FielData object
template< typename L >
void ForAllFieldData( L&& f )
{
    for ( intType i = 0; i != FieldData<int>::nData; i++ ) {
        f( i );
    }
}


}   // end namespace CFD

#endif // TYPES