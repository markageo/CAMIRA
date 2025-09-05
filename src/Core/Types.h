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
    Staircase, ImmersedBoundary
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

enum class TurbulenceModels {
    Null, Laminar, PrandtlZeroEquation, ZEQ0, ZEQ1, ZEQ2, ZEQ3, ZEQ4
};

enum class MultigridCycleType {
    V, F, W
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
    constexpr T &operator[](const typename enumStruct::ENUMDATA idx)
    { return m_data[idx]; }

    constexpr const T &operator[](const typename enumStruct::ENUMDATA idx) const 
    { return m_data[idx]; }


    // Only allowed when the object is holding axis data since the enum values should not change and are fairly standard
    constexpr T &operator[](const size_t idx)
    { 
        static_assert( std::is_same<enumStruct, CAMIRA::Axis>::value, "EnumVector indexed with integral value only allowed for Axis enums" );
        return m_data[idx]; 
    }

    constexpr const T &operator[](const size_t idx) const 
    { 
        static_assert( std::is_same<enumStruct, CAMIRA::Axis>::value, "EnumVector indexed with integral value only allowed for Axis enums" );
        return m_data[idx]; 
    }


    // Get pointer to data
    constexpr T* data()
    { return m_data; }

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



// ---------------------------------------------- Parallel Tensor SetConstant ------------------------------------------------- //

// These functions initialise the data in parallel, analogous to the setZero member function for Tensors.
// These are useful when running code on NUMA machines to make use of the first touch policy.
template< typename T >
[[maybe_unused]]
inline void SetTensorConstantParallel( T &tensor, 
                                 const floatType &value )
{
    static_assert(std::is_same<T, Tensor0D >::value ||
                  std::is_same<T, Tensor1D >::value ||
                  std::is_same<T, Tensor2D >::value ||
                  std::is_same<T, Tensor3D >::value );


    #pragma omp parallel for
    for ( intType i = 0; i != tensor.size(); i++ ) {
        tensor.data()[i] = value;
    }

}



template< typename T >
[[maybe_unused]]
inline void SetTensorZeroParallel( T &tensor )
{
    SetTensorConstantParallel( tensor, 0.0f );
}


// Overloads for FieldData objects of Tensors
template< typename T >
[[maybe_unused]]
inline void SetTensorConstantParallel( FieldData<T> &fieldData, 
                                       const floatType &value )
{
    ForAllFieldData( [&] (intType f) {
        SetTensorConstantParallel( fieldData[f], value );
    } );
}


template< typename T >
[[maybe_unused]]
inline void SetTensorZeroParallel( FieldData<T> &fieldData )
{
    SetTensorConstantParallel( fieldData, 0.0f );
}



// ------------------------------------------------ Parallel Tensor Copy ---------------------------------------------------- //

// Copies contents of one tensor to another in parallel.
// These are useful when running code on NUMA machines to make use of the first touch policy.
// To use this, need to allocate tensor being copied into, but not initialised to make use of first touch policy
template< typename T >
[[maybe_unused]]
inline void CopyTensorParallel( T &tensorCopy, 
                                const T &tensor )
{
    static_assert(std::is_same<T, Tensor0D >::value ||
                  std::is_same<T, Tensor1D >::value ||
                  std::is_same<T, Tensor2D >::value ||
                  std::is_same<T, Tensor3D >::value );


    #pragma omp parallel for
    for ( intType i = 0; i != tensor.size(); i++ ) {
        tensorCopy.data()[i] = tensor.data()[i];
    }
    
}



// This function will allocate the tensor being copied to such that it is the same size as the initial tensor
template< typename T >
[[maybe_unused]]
inline void SetAndCopyTensorParallel( T &tensorCopy, 
                                      const T &tensor )
{
    static_assert(std::is_same<T, Tensor0D >::value ||
                  std::is_same<T, Tensor1D >::value ||
                  std::is_same<T, Tensor2D >::value ||
                  std::is_same<T, Tensor3D >::value );

    tensorCopy = T( tensor.dimensions() );

    CopyTensorParallel( tensorCopy, tensor );
}



// Overloads for FieldData objects of Tensors
template< typename T >
[[maybe_unused]]
inline void CopyTensorParallel( FieldData<T> &fieldDataCopy, 
                                const FieldData<T> &fieldData )
{
    ForAllFieldData( [&] (intType f) {
        CopyTensorParallel( fieldDataCopy[f], fieldData[f] );
    } );
}

template< typename T >
[[maybe_unused]]
inline void SetAndCopyTensorParallel( FieldData<T> &fieldDataCopy, 
                                      const FieldData<T> &fieldData )
{
    ForAllFieldData( [&] (intType f) {
        fieldDataCopy[f] = T( fieldData[f].dimensions() );
    } );
    
    CopyTensorParallel( fieldDataCopy, fieldData );
}



}   // end namespace CAMIRA

#endif // TYPES