#include "Parallel.h"

namespace CFD
{


namespace 
{

std::vector<intType> RangeStrideVector( intType start, 
                                        intType end,
                                        intType stride )
{
    std::vector<intType> rangeStrideVector;

    auto boundsCheck = [&] (intType i) -> bool {
        if ( start <= end ) {
            return i < end;
        } else {
            return i > end;
        }
    }; 

    for ( intType i = start; boundsCheck(i); i += stride ) {
        rangeStrideVector.push_back( i );
    }
    return rangeStrideVector;
}    

}   // end anonymous namespace




std::vector< std::vector<intType> > CreateForward1DColourSet( const intType n )
{
    std::vector< std::vector<intType> > colorSet;

    // Red nodes
    colorSet.push_back( RangeStrideVector( 0, n, 3 ) );

    // Green nodes
    colorSet.push_back( RangeStrideVector( 1, n, 3 ) );

    // Blue nodes
    colorSet.push_back( RangeStrideVector( 2, n, 3 ) );

    return colorSet;
}



std::vector< std::vector<intType> > CreateReverse1DColourSet( const intType n )
{
    std::vector< std::vector<intType> > colorSet;

    // Red nodes
    intType redStart = n - 3;   // If divisibly by 3
    if ( n % 3 == 2 ) {
        redStart = n - 2;
    } else if ( n % 3 == 1 ) {
        redStart = n - 1;
    }
    colorSet.push_back( RangeStrideVector( redStart, -1, -3 ) );


    // Green nodes
    intType greenStart = redStart + 1;   // If divisibly by 3
    if ( n % 3 == 2 ) {
        greenStart = redStart + 1;
    } else if ( n % 3 == 1 ) {
        greenStart = redStart - 2;
    }
    colorSet.push_back( RangeStrideVector( greenStart, -1, -3 ) );

    // Blue nodes
    intType blueStart = greenStart + 1;   // If divisibly by 3
    if ( n % 3 == 2 ) {
        blueStart = redStart - 1;
    } else if ( n % 3 == 1 ) {
        blueStart = redStart - 1;
    }
    colorSet.push_back( RangeStrideVector( blueStart, -1, -3 ) );

    return colorSet;
}



} // end namespace CFD