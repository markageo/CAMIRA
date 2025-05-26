#ifndef TYPES
#define TYPES

#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Dense>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include "umpire/Allocator.hpp"
#include "umpire/ResourceManager.hpp"
#include "umpire/Umpire.hpp"
#include "RAJA/RAJA.hpp"


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


using RAJALayout1D = RAJA::Layout<1>;
using RAJALayout2D = RAJA::Layout<2, RAJA::PERM_JI >;   // Column major
using RAJALayout3D = RAJA::Layout<3, RAJA::PERM_KJI>;   // Column major

template <int, typename T>
struct View; 

// Specialization for 1D
template <typename T>
struct View<1, T> {
    using type = RAJA::View<T, RAJALayout1D>;
};

// Specialization for 2D
template <typename T>
struct View<2, T> {
    using type = RAJA::View<T, RAJALayout2D>;
};

// Specialization for 3D
template <typename T>
struct View<3, T> {
    using type = RAJA::View<T, RAJALayout3D>;
};

using View1D = View<1, floatType>;
using View2D = View<2, floatType>;
using View3D = View<3, floatType>;


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


// ---------------------------------------------------- EnumVector Class -------------------------------------------------- //

// Simple wrapper for c-style array that is safe to use in device code
template <typename enumStruct, typename T>
struct EnumVector
{
    static_assert(std::is_same<enumStruct, CFD::Axis                 >::value ||
                  std::is_same<enumStruct, CFD::BoundaryConditions   >::value ||
                  std::is_same<enumStruct, CFD::BoundaryPatches      >::value ||
                  std::is_same<enumStruct, CFD::TransportCoefficients>::value || 
                  std::is_same<enumStruct, CFD::Fields               >::value,
                  "Template parameter must be struct containing ENUMDATA type.");

    T m_data[enumStruct::count];


    // No explicit constructor/copy/destructor/move for aggregate type


    // Strong type indexing
    __host__ __device__
    constexpr T &operator[](const typename enumStruct::ENUMDATA idx)
    { return m_data[idx]; }

    __host__ __device__
    constexpr const T &operator[](const typename enumStruct::ENUMDATA idx) const 
    { return m_data[idx]; }


    // Only allowed when the object is holding axis data since the enum values should not change and are fairly standard
    __host__ __device__
    constexpr T &operator[](const size_t idx)
    { 
        static_assert( std::is_same<enumStruct, CFD::Axis>::value, "EnumVector indexed with integral value only allowed for Axis enums" );
        return m_data[idx]; 
    }

    __host__ __device__
    constexpr const T &operator[](const size_t idx) const 
    { 
        static_assert( std::is_same<enumStruct, CFD::Axis>::value, "EnumVector indexed with integral value only allowed for Axis enums" );
        return m_data[idx]; 
    }


    // Get pointer to data
    __host__ __device__
    constexpr T* data()
    { return m_data; }

    __host__ __device__
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


// ---------------------------------------------------- ArrayND Class -------------------------------------------------- //

// Class that wraps a RAJA::View and allocates its own memory to either DEVICE or HOST using umpire. 
template<int DIM, typename T = floatType>
class ArrayND {

static_assert( DIM <= 3, "ArrayND dimensions must be 3 or less." );

using ViewType = View<DIM, T>::type;

private:

    void* m_data;
    int m_allocatorId;
    ViewType m_view;

public:

    template<typename... Args>
    ArrayND( std::string resource,
                Args... dim_sizes)
    {
        static_assert( sizeof...(Args) == DIM, "Incorrect number of indices in ArrayND." );

        if ( resource != "HOST" || resource != "DEVICE" ) {
            throw std::runtime_error( "Resource for ArrayND must either be \"HOST\" or \"DEVICE\"." );
        }

        // Umpire allocator
        auto &rm = umpire::ResourceManager::getInstance();
        auto allocator = rm.getAllocator( resource ); 
        m_allocatorId = allocator.getId();

        // Allocate memory
        intType size = (dim_sizes * ...);
        m_data = static_cast<T*>( allocator.allocate( size * sizeof(T) ) );

        // Create View
        m_view = ViewType( m_data, dim_sizes... );
    }


    ArrayND()                                = default;
    ArrayND( ArrayND const & )            = default;
    ArrayND( ArrayND && )                 = default;
    ArrayND &operator=( ArrayND const & ) = default;
    ArrayND &operator=( ArrayND && )      = default;
   

    ~ArrayND() {
        if ( m_data ) {
            auto &rm = umpire::ResourceManager::getInstance();
            auto allocator = rm.getAllocator( m_allocatorId );
            allocator.deallocate(m_data);
        } 
    }


    // Access the underlying view
    __host__ __device__
    ViewType& view()
    { return m_view; }

    __host__ __device__
    ViewType& view() const
    { return m_view; }


    // Access to raw pointer
    __host__ __device__
    T* data()
    { return static_cast<T*>( m_data ); }


    // Total number of elements
    __host__ __device__
    intType size() const 
    { return m_view.size(); }


    // Extent of a particular dimension
    __host__ __device__
    intType dimension(intType dim) const { 
        assert( (void("Multiarray dimension outside extents of array."), dim >= 0 && dim < DIM) );
        return m_view.layout.extent(dim); 
    }


    // Multidimensional access
    template<typename... Args>
    __host__ __device__
    T& operator()(Args... indices) {
        static_assert(sizeof...(Args) == DIM, "Incorrect number of indices in ArrayND.");
        return m_view(indices...);
    }

    template<typename... Args>
    __host__ __device__
    const T& operator()(Args... indices) const {
        static_assert(sizeof...(Args) == DIM, "Incorrect number of indices in ArrayND.");
        return m_view(indices...);
    }

    // Set all elements to constant value - host only
    void setConstant( T value ) const {
        auto &rm = umpire::ResourceManager::getInstance();
        rm.memset( m_data, value );
    }


    // Set all elements to zero - host only
    void setZero() const {
        static_assert( std::is_arithmetic<T>::value, "SetZero member function only support for arithmetic types." );
        auto &rm = umpire::ResourceManager::getInstance();
        rm.memset( m_data, static_cast<T>(0) );
    }

};


}   // end namespace CFD

#endif // TYPES