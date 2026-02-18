#ifndef CAMIRA_FV_LOOKUPS
#define CAMIRA_FV_LOOKUPS

#include "Types.h"


namespace CAMIRA 
{

// Lookup Tables
namespace LUT
{

    // ----------------------------------------------------- Enum Lookups ----------------------------------------------------- //

    // Get index offset from TransportCoeffienct
    constexpr EnumVector<TransportCoefficients, intType> CoeffIndex( { 2,      // tt
                                                                       2,      // nn    
                                                                       2,      // ee
                                                                       1,      // t
                                                                       1,      // n
                                                                       1,      // e
                                                                       0,      // p
                                                                       -1,     // w
                                                                       -1,     // s
                                                                       -1,     // b
                                                                       -2,     // ww
                                                                       -2,     // ss
                                                                       -2} );  // bb


    // Get BoundaryPatches from Axis
    constexpr EnumVector<Axis, BoundaryPatches::ENUMDATA> PositivePatch( { BoundaryPatches::xPositive,
                                                                           BoundaryPatches::yPositive,
                                                                           BoundaryPatches::zPositive} );

    constexpr EnumVector<Axis, BoundaryPatches::ENUMDATA> NegativePatch( { BoundaryPatches::xNegative,
                                                                           BoundaryPatches::yNegative,
                                                                           BoundaryPatches::zNegative} );

    // Get TransportCoefficient from Axis
    constexpr EnumVector<Axis, TransportCoefficients::ENUMDATA> HiCoeff( { TransportCoefficients::e,
                                                                           TransportCoefficients::n,
                                                                           TransportCoefficients::t} );
                                                                    
    constexpr EnumVector<Axis, TransportCoefficients::ENUMDATA> HiHiCoeff( { TransportCoefficients::ee,
                                                                             TransportCoefficients::nn,
                                                                             TransportCoefficients::tt} );


    constexpr EnumVector<Axis, TransportCoefficients::ENUMDATA> LoCoeff( { TransportCoefficients::w,
                                                                           TransportCoefficients::s,
                                                                           TransportCoefficients::b} );

    constexpr EnumVector<Axis, TransportCoefficients::ENUMDATA> LoLoCoeff( { TransportCoefficients::ww,
                                                                             TransportCoefficients::ss,
                                                                             TransportCoefficients::bb} );                                                                 


    // Get Axis from BoundaryPatches
    constexpr EnumVector<BoundaryPatches, Axis::ENUMDATA> BoundaryPatchAxis( { Axis::X,      // xPositive
                                                                               Axis::X,      // xNegative
                                                                               Axis::Y,      // yPositive
                                                                               Axis::Y,      // yNegative
                                                                               Axis::Z,      // zPositive
                                                                               Axis::Z} );   // zNegative


    // Axis orthogonal to a given one
    constexpr EnumVector<Axis, Axis::ENUMDATA> LoOrthogonalAxis( { Axis::Y,
                                                                   Axis::X,
                                                                   Axis::X } );

    constexpr EnumVector<Axis, Axis::ENUMDATA> HiOrthogonalAxis( { Axis::Z,
                                                                   Axis::Z,
                                                                   Axis::Y } );
 
}   // end namespace LUT


}   // end namespace CAMIRA

#endif // CAMIRA_FV_LOOKUPS