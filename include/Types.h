#ifndef TYPES
#define TYPES

#include "Eigen/unsupported/CXX11/Tensor"
#include <type_traits>
#include <vector>
#include <memory>
#include <utility>
#include <iostream>

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
        p,  // (i  , j  , k  )
        n,  // (i  , j+1, k  )
        e,  // (i+1, j  , k  )
        s,  // (i  , j-1, k  )
        w,  // (i-1, j  , k  )
        t,  // (i  , j  , k+1)
        b,  // (i  , j  , k-1)
        nn, // (i  , j+2, k  )
        ee, // (i+2, j  , k  )
        ss, // (i  , j-2, k  )
        ww, // (i-2, j  , k  )
        tt, // (i  , j  , k+2)
        bb, // (i  , j  , k-2)
    };
    const static int count = 13;
};


// Lookup arrays for coefficients and patches coerresponding to each axis
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
constexpr std::array<Axis::ENUMDATA, 6> BoundaryPatchAxis{Axis::ENUMDATA::X,    // xPositive
                                                          Axis::ENUMDATA::X,    // xNegative
                                                          Axis::ENUMDATA::Y,    // yPositive
                                                          Axis::ENUMDATA::Y,    // yNegative
                                                          Axis::ENUMDATA::Z,    // zPositive
                                                          Axis::ENUMDATA::Z};   // zNegative

                                                

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

        // ----------------------------- Constructor: All arrays have same dimensions ----------------------------- //

        // 3D
        ArrayAllocator(const std::vector< ENUMDATA > &coeffs, const indexVector3 &dims) 
        requires ( std::is_same< arrayType, CFD::array3D >::value ) : 
            coeffPointers(enumStruct::count)
        {
            for (const auto &index : coeffs) {
                coeffPointers[index] = std::make_unique<CFD::array3D>( CFD::array3D( dims(0), dims(1), dims(2) ).setZero() );
            }
        }

        // 2D
        ArrayAllocator(const std::vector< ENUMDATA > &coeffs, const indexVector2 &dims) 
        requires ( std::is_same< arrayType, CFD::array2D >::value ) : 
            coeffPointers(enumStruct::count)
        {
            for (const auto &index : coeffs) {
                coeffPointers[index] = std::make_unique<CFD::array2D>( CFD::array2D( dims(0), dims(1) ).setZero() );
            }
        }

        // 1D
        ArrayAllocator(const std::vector< ENUMDATA > &coeffs, const intType &dim) 
        requires ( std::is_same< arrayType, CFD::array1D >::value ) : 
            coeffPointers(enumStruct::count)
        {
            for (const auto &index : coeffs) {
                coeffPointers[index] = std::make_unique<CFD::array1D>( CFD::array1D( dim ).setZero() );
            }
        }


        // --------------------------- Constructor: Arrays can have different dimensions --------------------------- //

        // 3D
        ArrayAllocator( const std::vector< std::pair< ENUMDATA, CFD::indexVector3 > > &arraySpec) 
        requires ( std::is_same< arrayType, CFD::array3D >::value ) : 
            coeffPointers(enumStruct::count)
        {
            CFD::indexVector3 dims;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dims = arraySpec[i].second;
                coeffPointers[ arraySpec[i].first ] = std::make_unique<CFD::array3D>( CFD::array3D( dims(0), dims(1), dims(2) ).setZero() );
            }
        }


        // 2D
        ArrayAllocator( const std::vector< std::pair< ENUMDATA, CFD::indexVector2 > > &arraySpec) 
        requires ( std::is_same< arrayType, CFD::array2D >::value ) : 
            coeffPointers(enumStruct::count)
        {
            CFD::indexVector2 dims;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dims = arraySpec[i].second;
                coeffPointers[ arraySpec[i].first ] = std::make_unique<CFD::array2D>( CFD::array2D( dims(0), dims(1) ).setZero() );
            }
        }


        // 1D
        ArrayAllocator( const std::vector< std::pair< ENUMDATA, CFD::intType > > &arraySpec) 
        requires ( std::is_same< arrayType, CFD::array1D >::value ) : 
            coeffPointers(enumStruct::count)
        {
            CFD::intType dim;
            for (size_t i = 0; i != arraySpec.size(); i++) {
                dim = arraySpec[i].second;
                coeffPointers[ arraySpec[i].first ] = std::make_unique<CFD::array1D>( CFD::array1D( dim ).setZero() );
            }
        }

        // ------------------------------------------- Copy Constructor ------------------------------------------- //

        ArrayAllocator(const ArrayAllocator &that) :
            coeffPointers(enumStruct::count)
        {
            // Allocate a new array object only if it was allocated in the original
            for (size_t i = 0; i != coeffPointers.size(); i++) {
                if (that.coeffPointers[i]) {
                    coeffPointers[i] = std::make_unique<arrayType>( *that.coeffPointers[i] );
                }
            }
        }


        // --------------------------------------- Copy Assignment Operator --------------------------------------- //

        ArrayAllocator &operator=(ArrayAllocator that)
        {   
            std::swap( this->coeffPointers, that.coeffPointers );
            return *this;
        }


        // ------------------------------------------- Move Constructor ------------------------------------------- //

        ArrayAllocator(ArrayAllocator&& that) noexcept :
            coeffPointers( std::move( that.coeffPointers ) )
        {}


        // --------------------------------------- Move Assignment Operator --------------------------------------- //

        ArrayAllocator &operator=(ArrayAllocator &&that) noexcept
        {
            return ArrayAllocator( std::move( that ) );
        }    


        // ---------------------------------------------- Destructor ---------------------------------------------- //

        ~ArrayAllocator() = default;


        // ----------------------------------- Array reference return operators ----------------------------------- //

        // For int parameters
        arrayType &operator[](const intType idx)
        {
            return *coeffPointers[idx];
        }

        arrayType &operator[](const intType idx) const 
        {
            return *coeffPointers[idx];
        }


        // ------------------------------------------- Container access ------------------------------------------- //

        std::vector< std::unique_ptr<arrayType> > &container()
        {
            return coeffPointers;
        }


    private:
        std::vector< std::unique_ptr<arrayType> > coeffPointers;

};

}   // end namespace CFD

#endif // TYPES