#ifndef STAGGER_INDEXING2
#define STAGGER_INDEXING2

#include "../Core/Types.h"
#include "../Core/FVLookups.h"
#include <type_traits>

namespace CFD
{

class StaggerIndexing2
{
    using TC = TransportCoefficients::ENUMDATA;

public:

        // Pressure coupled in the momentum equation
        struct MomentumPressure
        {
            MomentumPressure( Axis::ENUMDATA axis,
                              TransportCoefficients::ENUMDATA staggeredCoeff ) :
                cCoupled( OppositeCoeff( staggeredCoeff ) ),
                cOpposite( staggeredCoeff ),
                cLeft( LeftCoeff( axis, OppositeCoeff( staggeredCoeff ) ) ),
                cRight( RightCoeff( axis, OppositeCoeff( staggeredCoeff ) ) ),
    
                iCoupled( LUT::CoeffIndex[ cCoupled ] ),
                iOpposite( - iCoupled ),
                iLeft( LUT::CoeffIndex[ cLeft ] ),
                iRight( LUT::CoeffIndex[ cRight ] )
            {};
    
            TC cCoupled,  
               cOpposite,
               cLeft, 
               cRight;
    
            intType iCoupled,
                    iOpposite,
                    iLeft, 
                    iRight;
        };
    
    
        // Velocity coupled in the continuity equation
        struct ContinuityVelocity
        {
            ContinuityVelocity( Axis::ENUMDATA axis,
                                TransportCoefficients::ENUMDATA staggeredCoeff ) :
                cCoupled( staggeredCoeff ),
                cOpposite( OppositeCoeff(staggeredCoeff) ),
                cLeft( LeftCoeff( axis, staggeredCoeff ) ),
                cRight( RightCoeff( axis, staggeredCoeff ) ),
    
                iCoupled( LUT::CoeffIndex[ cCoupled ] ),
                iOpposite( - iCoupled ),
                iLeft( LUT::CoeffIndex[ cLeft ] ),
                iRight( LUT::CoeffIndex[ cRight ] )
            {}
    
            TC cCoupled,
               cOpposite,
               cLeft,
               cRight;
    
            intType iCoupled,
                    iOpposite,
                    iLeft, 
                    iRight;
        };

    StaggerIndexing2( Axis::ENUMDATA axis,
                      TransportCoefficients::ENUMDATA staggeredCoeff ) :
                        momentumPressure( axis, staggeredCoeff ),
                        continuityVelocity( axis, staggeredCoeff )
                    {};

    MomentumPressure momentumPressure;
    ContinuityVelocity continuityVelocity;

    
private:
    
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


    static constexpr TC LeftCoeff( Axis::ENUMDATA axis,
                                   TC coeff)
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Axis::ENUMDATA;

        if     ( coeff == p ) { 
            if      ( axis == X )
                return w;
            else if ( axis == Y )
                return s;
            else if ( axis == Z )
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


    static constexpr TC RightCoeff( Axis::ENUMDATA axis,
                                    TC coeff )
    {
        using enum TransportCoefficients::ENUMDATA;
        using enum Axis::ENUMDATA;

        if     ( coeff == p ) { 
            if      ( axis == X )
                return e;
            else if ( axis == Y )
                return n;
            else if ( axis == Z )
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
};


}   // end namespace CFD

#endif // STAGGER_INDEXING2