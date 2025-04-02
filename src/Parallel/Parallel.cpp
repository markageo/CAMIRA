#include "Parallel.h"
#include "../Solver/ArrayIndexConversions.h"

// DEBUGGING
#include "../IO/ArrayIO.h"

namespace CFD
{

RAJA::TypedIndexSet<RAJA::TypedRangeStrideSegment<intType>> CreateForward1DColourSet( const intType n )
{
    RAJA::TypedIndexSet<RAJA::TypedRangeStrideSegment<intType>> colorSet;

    // Red nodes
    colorSet.push_back( RAJA::TypedRangeStrideSegment<intType>( 0, n, 3 ) );

    // Green nodes
    colorSet.push_back( RAJA::TypedRangeStrideSegment<intType>( 1, n, 3 ) );

    // Blue nodes
    colorSet.push_back( RAJA::TypedRangeStrideSegment<intType>( 2, n, 3 ) );

    return colorSet;
}



RAJA::TypedIndexSet<RAJA::TypedRangeStrideSegment<intType>> CreateReverse1DColourSet( const intType n )
{
    RAJA::TypedIndexSet<RAJA::TypedRangeStrideSegment<intType>> colorSet;

    // Red nodes
    intType redStart = n - 3;   // If divisibly by 3
    if ( n % 3 == 2 ) {
        redStart = n - 2;
    } else if ( n % 3 == 1 ) {
        redStart = n - 1;
    }
    colorSet.push_back( RAJA::TypedRangeStrideSegment<intType>( redStart, -1, -3 ) );


    // Green nodes
    intType greenStart = redStart + 1;   // If divisibly by 3
    if ( n % 3 == 2 ) {
        greenStart = redStart + 1;
    } else if ( n % 3 == 1 ) {
        greenStart = redStart - 2;
    }
    colorSet.push_back( RAJA::TypedRangeStrideSegment<intType>( greenStart, -1, -3 ) );

    // Blue nodes
    intType blueStart = greenStart + 1;   // If divisibly by 3
    if ( n % 3 == 2 ) {
        blueStart = redStart - 1;
    } else if ( n % 3 == 1 ) {
        blueStart = redStart - 1;
    }
    colorSet.push_back( RAJA::TypedRangeStrideSegment<intType>( blueStart, -1, -3 ) );

    return colorSet;
}



RAJA::TypedIndexSet<RAJA::TypedListSegment<intType>> CreateForward3ColorSet( const iArray3 &nCells,
                                                                             camp::resources::Resource res )
{
    RAJA::TypedIndexSet<RAJA::TypedListSegment<intType>> colorSet;
    std::vector<intType> greenNodes, blueNodes, redNodes;

    const intType nColors = 3;
    enum { GREEN = 0, 
           BLUE = 1, 
           RED = 2 };
    intType startingColor          = GREEN,
            prevPlaneStartingColor = startingColor,
            currentColor  = startingColor;

    // DEBUGGING
    Tensor3D colorArray( nCells(0), nCells(1), nCells(2) );
    colorArray.setZero();

    for ( intType k = 0; k != nCells(2); k++ ) {

        for ( intType j = 0; j != nCells(1); j++ ) {

            currentColor = startingColor;

            for ( intType i = 0; i != nCells(0); i++ ) {

                intType idx = Sub2Ind( nCells, i, j, k );

                switch (currentColor)
                {
                    case GREEN:
                        greenNodes.push_back( idx );
                        colorArray(i, j, k) = GREEN;
                        break;

                    case BLUE:
                        blueNodes.push_back( idx );
                        colorArray(i, j, k) = BLUE;
                        break;

                    case RED:
                        redNodes.push_back( idx );
                        colorArray(i, j, k) = RED;
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

    // DEBUGGING
    WriteArray("color_array.dat", colorArray);


    colorSet.push_back( RAJA::TypedListSegment<intType>( redNodes.data()  , redNodes.size()  , res ) );
    colorSet.push_back( RAJA::TypedListSegment<intType>( greenNodes.data(), greenNodes.size(), res ) );
    colorSet.push_back( RAJA::TypedListSegment<intType>( blueNodes.data() , blueNodes.size() , res ) );

    return colorSet;
}



RAJA::TypedIndexSet<RAJA::TypedListSegment<intType>> CreateReverse3ColorSet( const iArray3 &nCells,
                                                                             camp::resources::Resource res )
{
    RAJA::TypedIndexSet<RAJA::TypedListSegment<intType>> colorSet;
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

    // Reverse the arrays
    std::reverse(redNodes.begin()  , redNodes.end());
    std::reverse(greenNodes.begin(), greenNodes.end());
    std::reverse(blueNodes.begin() , blueNodes.end());

    colorSet.push_back( RAJA::TypedListSegment<intType>( redNodes.data()  , redNodes.size()  , res ) );
    colorSet.push_back( RAJA::TypedListSegment<intType>( greenNodes.data(), greenNodes.size(), res ) );
    colorSet.push_back( RAJA::TypedListSegment<intType>( blueNodes.data() , blueNodes.size() , res ) );

    return colorSet;

}





RAJA::TypedIndexSet<RAJA::TypedListSegment<intType>> Create3ColorSetColumns( const iArray3 &nCells,
                                                                             camp::resources::Resource res )
{
    RAJA::TypedIndexSet<RAJA::TypedListSegment<intType>> colorSet;
    std::vector<intType> greenNodes, blueNodes, redNodes;

    const intType nColors = 3;
    enum { GREEN = 0, 
           BLUE = 1, 
           RED = 2 };
    intType startingColor          = GREEN,
            prevPlaneStartingColor = startingColor,
            currentColor  = startingColor;

    // DEBUGGING
    Tensor3D colorArray( nCells(0), nCells(1), nCells(2) );
    colorArray.setZero();

    for ( intType k = 0; k != nCells(2); k++ ) {

        for ( intType j = 0; j != nCells(1); j++ ) {

            currentColor = startingColor;

            for ( intType i = 0; i != nCells(0); i++ ) {

                intType idx = Sub2Ind( nCells, i, j, k );

                switch (currentColor)
                {
                    case GREEN:
                        greenNodes.push_back( idx );
                        colorArray(i, j, k) = GREEN;
                        break;

                    case BLUE:
                        blueNodes.push_back( idx );
                        colorArray(i, j, k) = BLUE;
                        break;

                    case RED:
                        redNodes.push_back( idx );
                        colorArray(i, j, k) = RED;
                        break;
                }

            }

            startingColor = ( startingColor + 1 ) % nColors;

        }

        startingColor = ( prevPlaneStartingColor + 1 ) % nColors;
        prevPlaneStartingColor = startingColor;
    }

    // DEBUGGING
    WriteArray("color_array.dat", colorArray);


    colorSet.push_back( RAJA::TypedListSegment<intType>( redNodes.data()  , redNodes.size()  , res ) );
    colorSet.push_back( RAJA::TypedListSegment<intType>( greenNodes.data(), greenNodes.size(), res ) );
    colorSet.push_back( RAJA::TypedListSegment<intType>( blueNodes.data() , blueNodes.size() , res ) );

    return colorSet;
}


RAJA::TypedIndexSet<RAJA::TypedListSegment<intType>> Create3ColorSetPlane( const iArray3 &nCells,
                                                                           camp::resources::Resource res )
{
    RAJA::TypedIndexSet<RAJA::TypedListSegment<intType>> colorSet;
    std::vector<intType> greenNodes, blueNodes, redNodes;

    const iArray2 nCellsPlane = { nCells(1), nCells(2) }; // j and k dimensions

    const intType nColors = 3;
    enum { GREEN = 0, 
           BLUE = 1, 
           RED = 2 };
    intType startingColor = GREEN,
            currentColor  = startingColor;

    // DEBUGGING
    Tensor2D colorArray( nCells(1), nCells(2) );
    colorArray.setZero();

    for ( intType k = 0; k != nCells(2); k++ ) {

        currentColor = startingColor;

        for ( intType j = 0; j != nCells(1); j++ ) {

            intType idx = Sub2Ind( nCellsPlane, j, k );

            switch (currentColor)
            {
                case GREEN:
                    greenNodes.push_back( idx );
                    colorArray(j, k) = GREEN;
                    break;

                case BLUE:
                    blueNodes.push_back( idx );
                    colorArray(j, k) = BLUE;
                    break;

                case RED:
                    redNodes.push_back( idx );
                    colorArray(j, k) = RED;
                    break;
            }

            // Go to next colour
            currentColor = (currentColor + 1) % nColors;

        }

        startingColor = ( startingColor + 1 ) % nColors;

    }

    // DEBUGGING
    WriteArray("color_array_plane.dat", colorArray);

    colorSet.push_back( RAJA::TypedListSegment<intType>( redNodes.data()  , redNodes.size()  , res ) );
    colorSet.push_back( RAJA::TypedListSegment<intType>( greenNodes.data(), greenNodes.size(), res ) );
    colorSet.push_back( RAJA::TypedListSegment<intType>( blueNodes.data() , blueNodes.size() , res ) );

    return colorSet;
}



} // end namespace CFD