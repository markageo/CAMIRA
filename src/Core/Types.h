#ifndef TYPES
#define TYPES

#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include "RAJA/RAJA.hpp"
#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"
#include "umpire/Umpire.hpp"

#include <type_traits>
#include <vector>
#include <utility>
#include <algorithm>


namespace CAMIRA 
{


// ----------------------------------------------------- Type Aliases ----------------------------------------------------- //

#ifdef CAMIRA_DOUBLE_PRECISION
    using floatType = double;
#else
    using floatType = float;
#endif

using intType = Eigen::Index;

template<int DIM>
using Tensor = Eigen::Tensor<floatType, DIM>;

using Tensor0D = Tensor<0>;    // Column major
using Tensor1D = Tensor<1>;    // Column major
using Tensor2D = Tensor<2>;    // Column major
using Tensor3D = Tensor<3>;    // Column major

using TensorIndex3D = Eigen::array<Eigen::Index, 3>;
using TensorIndex2D = Eigen::array<Eigen::Index, 2>;

using iArray3 = Eigen::Array<intType, 3, 1>;
using iArray2 = Eigen::Array<intType, 2, 1>;
using fArray3 = Eigen::Array<floatType, 3, 1>;
using fArray2 = Eigen::Array<floatType, 2, 1>;

using fVector3 = Eigen::Matrix<floatType, 3, 1>;
using fVector2 = Eigen::Matrix<floatType, 2, 1>;


template<int DIM, typename T>
class ArrayND;

using Array1D = ArrayND<1, floatType>;
using Array2D = ArrayND<2, floatType>;
using Array3D = ArrayND<3, floatType>;



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


struct Fields
{
    enum ENUMDATA
    {
        U, 
        V, 
        W, 
        P
    };
    const static int count = 3;
};


// For looping through enum and applying lambda to each element
template<typename enumStruct, typename L>
inline void EnumFor( L&& f )
{
    static_assert(std::is_same<enumStruct, Axis                 >::value ||
                  std::is_same<enumStruct, BoundaryConditions   >::value ||
                  std::is_same<enumStruct, BoundaryPatches      >::value ||
                  std::is_same<enumStruct, TransportCoefficients>::value ||
                  std::is_same<enumStruct, Fields               >::value   );

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


// ------------------------------------------------------ ArrayND Class --------------------------------------------------- //


namespace CAMIRA_INTERNAL 
{

    // Structs that are specialised to hold permutation array for RAJA::Views. These need to be std::arrays, 
    // but do not want to make them members of the ArrayND class.
    template<int>
    struct ViewPermutation;

    template<>
    struct ViewPermutation<1>
    { static constexpr std::array< intType, 1 > perm{{0}}; };

    template<>
    struct ViewPermutation<2>
    { static constexpr std::array< intType, 2 > perm{{1, 0}}; };

    template<>
    struct ViewPermutation<3>
    { static constexpr std::array< intType, 3 > perm{{2, 1, 0}}; };

}   // end namespace CAMIRA_INTERNAL


// Class that wraps a RAJA::View and allocates its own memory to either DEVICE or HOST using umpire. 
template<int DIM, typename T>
class ArrayND {

static_assert( DIM <= 3, "ArrayND dimensions must be 3 or less." );

using ViewType = RAJA::View<T, RAJA::Layout<2, intType, 0>>;

private:

    floatType* m_data;
    int m_allocatorId;  // Only hold the allocator id, not the allocator itself since it is not safe to use in device code
    ViewType m_view;

public:

    template<typename... Args>
    ArrayND( std::string resource,
             Args... dim_sizes) :
        m_data( nullptr ),
        m_allocatorId( -1 ),
        m_view( m_data, RAJA::make_permuted_layout( {{dim_sizes...}}, 
                                                    CAMIRA_INTERNAL::ViewPermutation<DIM>::perm ) )    // Need to set dimensions of view at construction
    {
        static_assert( sizeof...(Args) == DIM, "Incorrect number of indices in ArrayND." );

        if ( resource != "HOST" && resource != "DEVICE" ) {
            throw std::runtime_error( "Resource for ArrayND must either be \"HOST\" or \"DEVICE\"." );
        }

        // Umpire allocator
        auto &rm = umpire::ResourceManager::getInstance();
        auto allocator = rm.getAllocator( resource ); 
        m_allocatorId = allocator.getId();

        // Allocate memory
        intType size = (dim_sizes * ...);
        m_data = static_cast<T*>( allocator.allocate( size * sizeof(T) ) );

        // Reset view pointer
        m_view.set_data( m_data );
    }


    ArrayND()                             = default;
    ArrayND( ArrayND const & )            = default;
    ArrayND( ArrayND && )                 = default;
    ArrayND &operator=( ArrayND const & ) = default;
    ArrayND &operator=( ArrayND && )      = default;
    ~ArrayND()                            = default;


    void CleanUp() {
        if ( m_data ) {
            auto &rm = umpire::ResourceManager::getInstance();
            auto allocator = rm.getAllocator( m_allocatorId );
            allocator.deallocate(m_data);
            m_data = nullptr;
        }     
        m_view.set_data(m_data);
    }


    // Access the underlying view
    RAJA_HOST_DEVICE
    ViewType& view()
    { return m_view; }

    RAJA_HOST_DEVICE
    ViewType& view() const
    { return m_view; }


    // Access to raw pointer
    RAJA_HOST_DEVICE
    T* data() const
    { return static_cast<T*>( m_data ); }


    // Total number of elements
    RAJA_HOST_DEVICE
    intType size() const 
    { return m_view.size(); }


    // Extent of a particular dimension
    RAJA_HOST_DEVICE
    intType dimension(intType dim) const { 
        assert( (void("Multiarray dimension outside extents of array."), dim >= 0 && dim < DIM) );
        return m_view.layout.extent(dim); 
    }


    // Multidimensional access
    template<typename... Args>
    RAJA_HOST_DEVICE
    T& operator()(Args... indices) {
        static_assert(sizeof...(Args) == DIM, "Incorrect number of indices in ArrayND.");
        return m_view(indices...);
    }

    template<typename... Args>
    RAJA_HOST_DEVICE
    T& operator()(Args... indices) const {
        static_assert(sizeof...(Args) == DIM, "Incorrect number of indices in ArrayND.");
        return m_view(indices...);
    }

    // Copy underlying data from given array to this one
    void CopyDataFrom( const ArrayND<DIM, T> &sourceArray ) const 
    { 
        auto &rm = umpire::ResourceManager::getInstance();
        rm.copy( m_data, sourceArray.data() );  
    }


    // Set all elements to constant value - host only
    void SetConstant( T value ) const {
        auto &rm = umpire::ResourceManager::getInstance();
        rm.memset( m_data, value );
    }


    // Set all elements to zero - host only
    void SetZero() const {
        static_assert( std::is_arithmetic<T>::value, "SetZero member function only supported for arithmetic types." );
        auto &rm = umpire::ResourceManager::getInstance();
        rm.memset( m_data, static_cast<T>(0) );
    }

};


// ---------------------------------------------------- EnumVector Class -------------------------------------------------- //

// Simple wrapper for c-style array that is safe to use in device code
template <typename enumStruct, typename T>
struct EnumVector
{
    static_assert(std::is_same<enumStruct, CAMIRA::Axis                 >::value ||
                  std::is_same<enumStruct, CAMIRA::BoundaryConditions   >::value ||
                  std::is_same<enumStruct, CAMIRA::BoundaryPatches      >::value ||
                  std::is_same<enumStruct, CAMIRA::TransportCoefficients>::value || 
                  std::is_same<enumStruct, CAMIRA::Fields               >::value,
                  "Template parameter must be struct containing ENUMDATA type.");

    T m_data[enumStruct::count];


    // No explicit constructor/copy/destructor/move for aggregate type


    // Strong type indexing
    RAJA_HOST_DEVICE
    constexpr T &operator[](const typename enumStruct::ENUMDATA idx)
    { return m_data[idx]; }

    RAJA_HOST_DEVICE
    constexpr const T &operator[](const typename enumStruct::ENUMDATA idx) const 
    { return m_data[idx]; }


    // Only allowed when the object is holding axis data since the enum values should not change and are fairly standard
    RAJA_HOST_DEVICE
    constexpr T &operator[](const size_t idx)
    { 
        static_assert( std::is_same<enumStruct, CAMIRA::Axis>::value, "EnumVector indexed with integral value only allowed for Axis enums" );
        return m_data[idx]; 
    }

    RAJA_HOST_DEVICE
    constexpr const T &operator[](const size_t idx) const 
    { 
        static_assert( std::is_same<enumStruct, CAMIRA::Axis>::value, "EnumVector indexed with integral value only allowed for Axis enums" );
        return m_data[idx]; 
    }


    // Get pointer to data
    RAJA_HOST_DEVICE
    constexpr T* data()
    { return m_data; }

    RAJA_HOST_DEVICE
    constexpr const T* data() const
    { return m_data; }

};


// ---------------------------------------------------- FieldData Class -------------------------------------------------- //


// A general struct for holding values corresponding to different fields
template < typename dataType >
struct FieldData {
    EnumVector<Axis, dataType> U;
    dataType P;
    
    FieldData() { SetPointers(); };
    FieldData( const dataType &data ) : U{ data, data, data }, P( data ) 
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
        { m_dataPointers = { &U[Axis::X], &U[Axis::Y], &U[Axis::Z], &P }; } // Must match the enums in the Fields struct
};


// To allow iterating through values of a FielData object
template< typename L >
void ForAllFieldData( L&& f )
{
    for ( intType i = 0; i != FieldData<int>::nData; i++ ) {
        f( i );
    }
}



}   // end namespace CAMIRA

#endif // TYPES