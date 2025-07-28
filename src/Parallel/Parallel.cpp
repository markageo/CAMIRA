#include "Parallel.h"
#include "../Core/ArrayIndexConversions.h"

namespace CAMIRA
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




std::vector< std::vector<intType> > CreateForward1DColorSet( const intType n )
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



std::vector< std::vector<intType> > CreateReverse1DColorSet( const intType n )
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



std::vector< std::vector<intType> > CreateForward2DColorSet( const iArray3 &nCells )
{
    std::vector< std::vector<intType> > colorSet;
    std::vector<intType> greenNodes, blueNodes, redNodes;

    const iArray2 nCellsPlane = { nCells(1), nCells(2) }; // j and k dimensions

    const intType nColors = 3;
    enum { GREEN = 0, 
           BLUE = 1, 
           RED = 2 };
    intType startingColor = GREEN,
            currentColor  = startingColor;

    for ( intType k = 0; k != nCells(2); k++ ) {

        currentColor = startingColor;

        for ( intType j = 0; j != nCells(1); j++ ) {

            intType idx = Sub2Ind( nCellsPlane, j, k );

            switch (currentColor)
            {
                case GREEN:
                    greenNodes.push_back( idx );
                    break;

                case BLUE:
                    blueNodes.push_back( idx );
                    break;

                case RED:
                    redNodes.push_back( idx );
                    break;
            }

            // Go to next colour
            currentColor = (currentColor + 1) % nColors;
        }

        startingColor = ( startingColor + 1 ) % nColors;
    }

    colorSet.push_back( redNodes );
    colorSet.push_back( greenNodes );
    colorSet.push_back( blueNodes );

    return colorSet;
}



std::vector< std::vector<intType> > CreateReverse2DColorSet( const iArray3 &nCells )
{
    std::vector< std::vector<intType> > colorSet = CreateForward2DColorSet( nCells );

    // Reverse the arrays
    for ( auto &color : colorSet ) {
        std::reverse( color.begin(), color.end() );
    }

    return colorSet;
}



std::vector< std::vector<intType> > CreateForward3DColorSet( const iArray3 &nCells )
{
    std::vector< std::vector<intType> > colorSet;
    std::vector<intType> greenNodes, blueNodes, redNodes;

    const intType nColors = 3;
    enum { GREEN = 0, 
           BLUE = 1, 
           RED = 2 };
    intType startingColor          = GREEN,
            prevPlaneStartingColor = startingColor,
            currentColor  = startingColor;

    for ( intType k = 0; k != nCells(2); k++ ) {

        for ( intType j = 0; j != nCells(1); j++ ) {

            currentColor = startingColor;

            for ( intType i = 0; i != nCells(0); i++ ) {

                intType idx = Sub2Ind( nCells, i, j, k );

                switch (currentColor)
                {
                    case GREEN:
                        greenNodes.push_back( idx );
                        break;

                    case BLUE:
                        blueNodes.push_back( idx );
                        break;

                    case RED:
                        redNodes.push_back( idx );
                        break;
                }

                // Go to next colour
                currentColor = (currentColor + 1) % nColors;
            }

            startingColor = ( startingColor + 1 ) % nColors;

        }

        startingColor = ( prevPlaneStartingColor + 1 ) % nColors;
        prevPlaneStartingColor = startingColor;
    }

    colorSet.push_back( redNodes );
    colorSet.push_back( greenNodes );
    colorSet.push_back( blueNodes );

    return colorSet;
}



std::vector< std::vector<intType> > CreateReverse3DColorSet( const iArray3 &nCells )
{
    std::vector< std::vector<intType> > colorSet = CreateForward3DColorSet( nCells );

    // Reverse the arrays
    for ( auto &color : colorSet ) {
        std::reverse( color.begin(), color.end() );
    }

    return colorSet;
}




} // end namespace CAMIRA