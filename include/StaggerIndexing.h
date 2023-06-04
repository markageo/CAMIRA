#include "Types.h"
#include <type_traits>

namespace CFD
{

template< Fields::ENUMDATA field, 
          TransportCoefficients::ENUMDATA staggeredCoeff >
class StaggerIndexing
{
    using TC = TransportCoefficients::ENUMDATA;
    using F = Fields::ENUMDATA;

    static_assert( ( field == F::U && staggeredCoeff == TC::e ) ||
                   ( field == F::U && staggeredCoeff == TC::w ) ||
                   ( field == F::U && staggeredCoeff == TC::p ) ||
                   
                   ( field == F::V && staggeredCoeff == TC::n ) ||
                   ( field == F::V && staggeredCoeff == TC::s ) ||
                   ( field == F::V && staggeredCoeff == TC::p ) ||
                   
                   ( field == F::W && staggeredCoeff == TC::t ) ||
                   ( field == F::W && staggeredCoeff == TC::b ) ||
                   ( field == F::W && staggeredCoeff == TC::p )  );
    
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
        using enum Fields::ENUMDATA;

        if     ( coeff == p ) { 
            if      ( field == U )
                return w;
            else if ( field == V)
                return s;
            else if ( field == W )
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
        using enum Fields::ENUMDATA;

        if     ( coeff == p ) { 
            if      ( field == U )
                return e;
            else if ( field == V)
                return n;
            else if ( field == W )
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

            static constexpr intType iCoupled = CoeffIndex[ cCoupled ],
                                     iLeft    = CoeffIndex[ cLeft ], 
                                     iRight   = CoeffIndex[ cRight ];
        };


        // Velocity coupled in the continuity equation
        struct ContinuityVelocity
        {
            static constexpr TC cCoupled = staggeredCoeff,
                                cLeft    = StaggerIndexing::LeftCoeff( staggeredCoeff ),
                                cRight   = StaggerIndexing::RightCoeff( staggeredCoeff );

            static constexpr intType iCoupled = CoeffIndex[ cCoupled ],
                                     iLeft    = CoeffIndex[ cLeft ], 
                                     iRight   = CoeffIndex[ cRight ];
        };

};


}   // end namespace CFD