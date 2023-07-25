#include "../Types.h"
#include <type_traits>

namespace CFD
{

template< Axis::ENUMDATA field, 
          TransportCoefficients::ENUMDATA staggeredCoeff >
class StaggerIndexing
{
    using TC = TransportCoefficients::ENUMDATA;

    static_assert( ( field == Axis::X && staggeredCoeff == TC::e ) ||
                   ( field == Axis::X && staggeredCoeff == TC::w ) ||
                   ( field == Axis::X && staggeredCoeff == TC::p ) ||
                   
                   ( field == Axis::Y && staggeredCoeff == TC::n ) ||
                   ( field == Axis::Y && staggeredCoeff == TC::s ) ||
                   ( field == Axis::Y && staggeredCoeff == TC::p ) ||
                   
                   ( field == Axis::Z && staggeredCoeff == TC::t ) ||
                   ( field == Axis::Z && staggeredCoeff == TC::b ) ||
                   ( field == Axis::Z && staggeredCoeff == TC::p )  );
    
    static constexpr TC OppositeCoeff(TC coeff)
    {
        using enum TransportCoefficients::ENUMDATA;

        if      ( coeff == p )
            return p;
        else if ( coeff == n )
            return s;
        else if ( coeff == s )
            return n;
        else if ( coeff == e)
            return w;
        else if ( coeff == w )
            return e;
        else if ( coeff == t )
            return b;
        else if ( coeff == b )
            return t;   
        else    
            return p;
    }


    static constexpr TC LeftCoeff(TC coeff)
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Axis::ENUMDATA;

        if     ( coeff == p ) { 
            if      ( field == X )
                return w;
            else if ( field == Y )
                return s;
            else if ( field == Z )
                return b;
        } 
        else if ( coeff == n )
            return s;
        else if ( coeff == s )
            return p;
        else if ( coeff == e )
            return w;
        else if ( coeff == w )
            return p;
        else if ( coeff == t )
            return b;
        else if ( coeff == b )
            return p;
        else 
            return p;
    }


    static constexpr TC RightCoeff(TC coeff)
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Axis::ENUMDATA;

        if     ( coeff == p ) { 
            if      ( field == X )
                return e;
            else if ( field == Y )
                return n;
            else if ( field == Z )
                return t;
        } 
        else if ( coeff == n )
            return p;
        else if ( coeff == s )
            return n;
        else if ( coeff == e )
            return p;
        else if ( coeff == w )
            return e;
        else if ( coeff == t )
            return p;
        else if ( coeff == b )
            return t;
        else 
            return p;
    }


    public:

        // Pressure coupled in the momentum equation
        struct MomentumPressure
        {
            static constexpr TC cCoupled = StaggerIndexing::OppositeCoeff( staggeredCoeff ),
                                cLeft    = StaggerIndexing::LeftCoeff( StaggerIndexing::OppositeCoeff( staggeredCoeff ) ), 
                                cRight   = StaggerIndexing::RightCoeff( StaggerIndexing::OppositeCoeff( staggeredCoeff ) );

            static constexpr intType iCoupled = LUT::CoeffIndex[ cCoupled ],
                                     iLeft    = LUT::CoeffIndex[ cLeft ], 
                                     iRight   = LUT::CoeffIndex[ cRight ];
        };


        // Velocity coupled in the continuity equation
        struct ContinuityVelocity
        {
            static constexpr TC cCoupled = staggeredCoeff,
                                cLeft    = StaggerIndexing::LeftCoeff( staggeredCoeff ),
                                cRight   = StaggerIndexing::RightCoeff( staggeredCoeff );

            static constexpr intType iCoupled = LUT::CoeffIndex[ cCoupled ],
                                     iLeft    = LUT::CoeffIndex[ cLeft ], 
                                     iRight   = LUT::CoeffIndex[ cRight ];
        };

};


}   // end namespace CFD