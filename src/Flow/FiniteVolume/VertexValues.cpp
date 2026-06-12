#include "FiniteVolume.h"
#include "Core/FVLookups.h"
#include "Core/FVTools.h"

namespace CAMIRA
{

using namespace CORE;

namespace FLOW
{

namespace 
{
    // Weighted linear interpolation of the cell centers to the cell verticies of the interior points
    void InterpolateInteriorPoints( Tensor3D &vertexField, 
                                    const Tensor3D &field, 
                                    const Mesh &mesh )
    {
        using enum Axis::ENUMDATA;
        using FVT::G;

        for ( intType k = 0; k != vertexField.dimension(2); k++ ) {
            for ( intType j = 0; j != vertexField.dimension(1); j++ ) {
                for ( intType i = 0; i != vertexField.dimension(0); i++ ) {

                    // Points to interpolate from
                    floatType c000 = field( G(i-1, j-1, k-1) ),
                              c100 = field( G(i  , j-1, k-1) ),
                              c010 = field( G(i-1, j  , k-1) ),
                              c001 = field( G(i-1, j-1, k  ) ),
                              c101 = field( G(i  , j-1, k  ) ),
                              c011 = field( G(i-1, j  , k  ) ),
                              c110 = field( G(i  , j  , k-1) ),
                              c111 = field( G(i  , j  , k  ) );

                    // Linear interpolation in z direction
                    floatType lambdaZ = mesh.interpFactors[Z](k);
                    floatType c00 = ( 1-lambdaZ ) * c000  +  lambdaZ * c001,
                              c10 = ( 1-lambdaZ ) * c100  +  lambdaZ * c101,
                              c01 = ( 1-lambdaZ ) * c010  +  lambdaZ * c011,
                              c11 = ( 1-lambdaZ ) * c110  +  lambdaZ * c111;

                    
                    // Linear interpolation in y direction
                    floatType lambdaY = mesh.interpFactors[Y](j);
                    floatType c0 = ( 1-lambdaY ) * c00  +  lambdaY * c01,
                              c1 = ( 1-lambdaY ) * c10  +  lambdaY * c11;


                    // Linear interpolation in x direction
                    floatType lambdaX = mesh.interpFactors[X](i);
                    floatType c = ( 1-lambdaX ) * c0  +  lambdaX * c1;

                    vertexField(i, j, k) = c;

                }
            }
        }

    }



    floatType MaskInterpFactor( const floatType interpFactor,
                                const floatType loMask,
                                const floatType hiMask )
    {
        if ( loMask == 0.0f ) {
            return 1.0f;
        } else if ( hiMask == 0.0f ) {
            return 0.0f;
        } 

        return interpFactor;
    }



    // Weighted linear interpolation of the cell centers to the cell verticies of the interior points. Solid cells are not used
    // in the interpolation 
    void InterpolateInteriorPointsWithMasking( Tensor3D &vertexField, 
                                               const Tensor3D &field, 
                                               const Tensor3D &mask,
                                               const Mesh &mesh )
    {
        using enum Axis::ENUMDATA;
        using FVT::G;

        for ( intType k = 0; k != vertexField.dimension(2); k++ ) {
            for ( intType j = 0; j != vertexField.dimension(1); j++ ) {
                for ( intType i = 0; i != vertexField.dimension(0); i++ ) {

                    // Points to interpolate from
                    floatType c000 = field( G(i-1, j-1, k-1) ),
                              c100 = field( G(i  , j-1, k-1) ),
                              c010 = field( G(i-1, j  , k-1) ),
                              c001 = field( G(i-1, j-1, k  ) ),
                              c101 = field( G(i  , j-1, k  ) ),
                              c011 = field( G(i-1, j  , k  ) ),
                              c110 = field( G(i  , j  , k-1) ),
                              c111 = field( G(i  , j  , k  ) );

                    // Mask values for the corresponding cells
                    floatType m000 = mask( G(i-1, j-1, k-1) ),
                              m100 = mask( G(i  , j-1, k-1) ),
                              m010 = mask( G(i-1, j  , k-1) ),
                              m001 = mask( G(i-1, j-1, k  ) ),
                              m101 = mask( G(i  , j-1, k  ) ),
                              m011 = mask( G(i-1, j  , k  ) ),
                              m110 = mask( G(i  , j  , k-1) ),
                              m111 = mask( G(i  , j  , k  ) );

                    // Linear interpolation in z direction
                    floatType lambdaZ = mesh.interpFactors[Z](k);

                    floatType lambdaZ00 = MaskInterpFactor( lambdaZ, m000, m001 ),
                              lambdaZ10 = MaskInterpFactor( lambdaZ, m100, m101 ),
                              lambdaZ01 = MaskInterpFactor( lambdaZ, m010, m011 ),
                              lambdaZ11 = MaskInterpFactor( lambdaZ, m110, m111 );

                    floatType c00 = ( 1-lambdaZ00 ) * c000  +  lambdaZ00 * c001,
                              c10 = ( 1-lambdaZ10 ) * c100  +  lambdaZ10 * c101,
                              c01 = ( 1-lambdaZ01 ) * c010  +  lambdaZ01 * c011,
                              c11 = ( 1-lambdaZ11 ) * c110  +  lambdaZ11 * c111;

                    floatType m00 = m000 + m001 > 0.0f,
                              m10 = m100 + m101 > 0.0f,
                              m01 = m010 + m011 > 0.0f,
                              m11 = m110 + m111 > 0.0f;

                    // Linear interpolation in y direction
                    floatType lambdaY = mesh.interpFactors[Y](j);

                    floatType lambdaY0 = MaskInterpFactor( lambdaY, m00, m01 ),
                              lambdaY1 = MaskInterpFactor( lambdaY, m10, m11 );

                    floatType c0 = ( 1-lambdaY0 ) * c00  +  lambdaY0 * c01,
                              c1 = ( 1-lambdaY1 ) * c10  +  lambdaY1 * c11;

                    floatType m0 = m00 + m01 > 0.0f,
                              m1 = m10 + m11 > 0.0f;

                    // Linear interpolation in x direction
                    floatType lambdaX = mesh.interpFactors[X](i);

                    lambdaX = MaskInterpFactor( lambdaX, m0, m1 );

                    floatType c = ( 1-lambdaX ) * c0  +  lambdaX * c1;

                    vertexField(i, j, k) = c;

                }
            }
        }

    }



    // Does a weighted linear average to set the edge values
    void WeightedAverageEdges( Tensor3D &vertexField,
                               const Mesh &mesh )
    {

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
            Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[ axis ];
            Axis::ENUMDATA axis2 = LUT::HiOrthogonalAxis[ axis ];

            // Indexing for each edge
            std::array<intType, 2> iVals{ 0, vertexField.dimension(axis1)-1 },
                                   jVals{ 0, vertexField.dimension(axis2)-1 },
                                   iCellIndex{ 0, vertexField.dimension(axis1)-2 },
                                   jCellIndex{ 0, vertexField.dimension(axis2)-2 },
                                   nghbr{ +1, -1 };


            // Iterate each edge
            for ( size_t iIndex = 0; iIndex != 2; iIndex++ ) {
                for ( size_t jIndex = 0; jIndex != 2; jIndex++ ) {
                    
                    intType i = iVals[ iIndex ];
                    intType j = jVals[ jIndex ];

                    // Weighting factors
                    floatType length = ( mesh.cellLengths[axis1]( iCellIndex[ iIndex ] ) + mesh.cellLengths[axis2]( jCellIndex[ jIndex ] ) );
                    floatType wfN1 = ( length - mesh.cellLengths[axis1]( iCellIndex[ iIndex ] ) ) / length;
                    floatType wfN2 = ( length - mesh.cellLengths[axis2]( jCellIndex[ jIndex ] ) ) / length;

                    // Index of node at vertex
                    TensorIndex3D idx;
                    idx[axis]  = 0;
                    idx[axis1] = i;
                    idx[axis2] = j;

                    // Index neighbour in the i direction
                    TensorIndex3D idxN1;
                    idxN1[axis]  = 0;
                    idxN1[axis1] = i + nghbr[ iIndex ];
                    idxN1[axis2] = j;

                    // Index neighbour in the j direction
                    TensorIndex3D idxN2;
                    idxN2[axis]  = 0;
                    idxN2[axis1] = i;
                    idxN2[axis2] = j + nghbr[ jIndex ];


                    for ( intType k = 0; k != vertexField.dimension(axis); k++ ) {

                        idx[axis]   = k;
                        idxN1[axis] = k;
                        idxN2[axis] = k;

                        vertexField( idx ) =  wfN1 * vertexField( idxN1 )  +  wfN2 * vertexField( idxN2 );
                    }

                }
            }

        } );

    }



    // Does a weighted linear average to set the corner values
    void WeightedAverageCorners( Tensor3D &vertexField, 
                                 const Mesh &mesh )
    {
        using enum Axis::ENUMDATA;

        std::array<intType, 2> iVals{ 0, vertexField.dimension(X)-1 },
                               jVals{ 0, vertexField.dimension(Y)-1 },
                               kVals{ 0, vertexField.dimension(Z)-1 },
                               iCellIndex{ 0, vertexField.dimension(X)-2 },
                               jCellIndex{ 0, vertexField.dimension(Y)-2 },
                               kCellIndex{ 0, vertexField.dimension(Z)-2 },
                               nghbr{ +1, -1 };

        // Iterate each corner
        for ( size_t kIndex = 0; kIndex != 2; kIndex++ ) {
            for ( size_t jIndex = 0; jIndex != 2; jIndex++ ) {
                for ( size_t iIndex = 0; iIndex != 2; iIndex++ ) {

                    // Weighting factors
                    floatType length = mesh.cellLengths[X]( iCellIndex[ iIndex ] ) 
                                     + mesh.cellLengths[Y]( jCellIndex[ jIndex ] ) 
                                     + mesh.cellLengths[Z]( kCellIndex[ kIndex ] );
                    floatType wfNi = ( length - mesh.cellLengths[X]( iCellIndex[ iIndex ] ) ) / length;
                    floatType wfNj = ( length - mesh.cellLengths[Y]( jCellIndex[ jIndex ] ) ) / length;
                    floatType wfNk = ( length - mesh.cellLengths[Z]( kCellIndex[ kIndex ] ) ) / length;

                    // Index of corner point
                    TensorIndex3D idx{ iVals[ iIndex ], jVals[ jIndex ], kVals[ kIndex ] };

                    // Index of neighbouring points
                    TensorIndex3D idxNi{ idx }, idxNj{ idx }, idxNk{ idx };
                    idxNi[X] += nghbr[ iIndex ];
                    idxNj[Y] += nghbr[ jIndex ];
                    idxNk[Z] += nghbr[ kIndex ];

                    vertexField( idx ) = wfNi * vertexField( idxNi )
                                       + wfNj * vertexField( idxNj )
                                       + wfNk * vertexField( idxNk );
                }
            }
        }
    }


}   // end anonymous namespace




Tensor3D InterpolateToVertex( const Tensor3D &field,
                              const Mesh &mesh )
{
    Tensor3D vertexField( mesh.nFacesNormal[0][0], mesh.nFacesNormal[1][1], mesh.nFacesNormal[2][2] );
    vertexField.setZero();

    InterpolateInteriorPoints( vertexField, field, mesh );

    WeightedAverageEdges( vertexField, mesh );

    WeightedAverageCorners( vertexField, mesh );

    return vertexField;
}



Tensor3D InterpolateToVertexWithMasking( const Tensor3D &field,
                                         const Mesh &mesh, 
                                         const Tensor3D &mask )
{
    Tensor3D vertexField( mesh.nFacesNormal[0][0], mesh.nFacesNormal[1][1], mesh.nFacesNormal[2][2] );
    vertexField.setZero();

    InterpolateInteriorPointsWithMasking( vertexField, field, mask, mesh );

    WeightedAverageEdges( vertexField, mesh );

    WeightedAverageCorners( vertexField, mesh );

    return vertexField;
}



}   // end namespace FLOW

}   // end namespace CAMIRA