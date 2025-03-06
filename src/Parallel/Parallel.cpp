#include "Parallel.h"
#include "../Solver/ArrayIndexConversions.h"



namespace CFD
{

RAJA::TypedIndexSet<RAJA::TypedListSegment<intType>> Create3ColorSet( const iArray3 &nCells,
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

    colorSet.push_back( RAJA::TypedListSegment<intType>( redNodes.data()  , redNodes.size()  , res ) );
    colorSet.push_back( RAJA::TypedListSegment<intType>( greenNodes.data(), greenNodes.size(), res ) );
    colorSet.push_back( RAJA::TypedListSegment<intType>( blueNodes.data() , blueNodes.size() , res ) );

    return colorSet;
}

} // end namespace CFD