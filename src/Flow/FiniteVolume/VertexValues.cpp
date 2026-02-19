#include "FiniteVolume.h"
#include "Core/FVLookups.h"

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

        for ( intType k = 1; k != vertexField.dimension(2)-1; k++ ) {
            for ( intType j = 1; j != vertexField.dimension(1)-1; j++ ) {
                for ( intType i = 1; i != vertexField.dimension(0)-1; i++ ) {

                    // Points to interpolation from
                    floatType c000 = field( i-1, j-1, k-1 ),
                              c100 = field( i  , j-1, k-1 ),
                              c010 = field( i-1, j  , k-1 ),
                              c001 = field( i-1, j-1, k   ),
                              c101 = field( i  , j-1, k   ),
                              c011 = field( i-1, j  , k   ),
                              c110 = field( i  , j  , k-1 ),
                              c111 = field( i  , j  , k   );

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

        for ( intType k = 1; k != vertexField.dimension(2)-1; k++ ) {
            for ( intType j = 1; j != vertexField.dimension(1)-1; j++ ) {
                for ( intType i = 1; i != vertexField.dimension(0)-1; i++ ) {

                    // Points to interpolate from
                    floatType c000 = field( i-1, j-1, k-1 ),
                              c100 = field( i  , j-1, k-1 ),
                              c010 = field( i-1, j  , k-1 ),
                              c001 = field( i-1, j-1, k   ),
                              c101 = field( i  , j-1, k   ),
                              c011 = field( i-1, j  , k   ),
                              c110 = field( i  , j  , k-1 ),
                              c111 = field( i  , j  , k   );

                    // Mask values for the corresponding cells
                    floatType m000 = mask( i-1, j-1, k-1 ),
                              m100 = mask( i  , j-1, k-1 ),
                              m010 = mask( i-1, j  , k-1 ),
                              m001 = mask( i-1, j-1, k   ),
                              m101 = mask( i  , j-1, k   ),
                              m011 = mask( i-1, j  , k   ),
                              m110 = mask( i  , j  , k-1 ),
                              m111 = mask( i  , j  , k   );

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



    Tensor2D CellPlaneToVertexPlane( const Tensor2D &cellPlane, 
                                     const Tensor1D &interpFactors0,
                                     const Tensor1D &interpFactors1,
                   [[maybe_unused]]  const Tensor2D &maskPlane )
    {
        Tensor2D vertexPlane( cellPlane.dimension(0)+1, cellPlane.dimension(1)+1 );

        for ( intType j = 1; j != vertexPlane.dimension(1)-1; j++ ) {
            for ( intType i = 1; i != vertexPlane.dimension(0)-1; i++ ) {

                // Points to interpolation from
                floatType c00 = cellPlane( i-1, j-1 ),
                          c10 = cellPlane( i  , j-1 ),
                          c01 = cellPlane( i-1, j   ),
                          c11 = cellPlane( i  , j   );

                // Linear interpolation in j direction
                floatType lambdaY = interpFactors1(j);
                floatType c0 = ( 1-lambdaY ) * c00  +  lambdaY * c01,
                          c1 = ( 1-lambdaY ) * c10  +  lambdaY * c11;

                // Linear interpolation in i direction
                floatType lambdaX = interpFactors0(i);
                floatType c = ( 1-lambdaX ) * c0  +  lambdaX * c1;

                vertexPlane(i, j) = c;

            }
        }

        return vertexPlane;
    }



    Tensor2D CellPlaneToVertexPlaneWithMasking( const Tensor2D &cellPlane, 
                                                const Tensor1D &interpFactors0,
                                                const Tensor1D &interpFactors1,
                                                const Tensor2D &maskPlane )
    {
        Tensor2D vertexPlane( cellPlane.dimension(0)+1, cellPlane.dimension(1)+1 );

        for ( intType j = 1; j != vertexPlane.dimension(1)-1; j++ ) {
            for ( intType i = 1; i != vertexPlane.dimension(0)-1; i++ ) {

                // Points to interpolate from
                floatType c00 = cellPlane( i-1, j-1 ),
                          c10 = cellPlane( i  , j-1 ),
                          c01 = cellPlane( i-1, j   ),
                          c11 = cellPlane( i  , j   );

                // Mask values for corresponding cells
                floatType m00 = maskPlane( i-1, j-1 ),
                          m10 = maskPlane( i  , j-1 ),
                          m01 = maskPlane( i-1, j   ),
                          m11 = maskPlane( i  , j   );

                // Linear interpolation in j direction
                floatType lambdaY = interpFactors1(j);

                floatType lambdaY0 = MaskInterpFactor( lambdaY, m00, m01 ),
                          lambdaY1 = MaskInterpFactor( lambdaY, m10, m11 );

                floatType c0 = ( 1-lambdaY0 ) * c00  +  lambdaY0 * c01,
                          c1 = ( 1-lambdaY1 ) * c10  +  lambdaY1 * c11;

                floatType m0 = m00 + m01 > 0.0f,
                          m1 = m10 + m11 > 0.0f;

                // Linear interpolation in i direction
                floatType lambdaX = interpFactors0(i);

                lambdaX = MaskInterpFactor( lambdaX, m0, m1 );

                floatType c = ( 1-lambdaX ) * c0  +  lambdaX * c1;

                vertexPlane(i, j) = c;

            }
        }

        return vertexPlane;
    }



    // Set the face values of the cell verticies using the boundary conditions. Corners/edges will be incorrect, these will be set individually later
    void SetBoundaryFaces( Tensor3D &vertexField,
                           const Tensor3D &field,
                           const Tensor3D &mask,
                           const Mesh &mesh,
                           const BoundaryConditionData::Patches &boundaryConditions,
                           Tensor2D (*CellPlaneToVertexPlaneFunc)(const Tensor2D &, const Tensor1D &, const Tensor1D &, const Tensor2D & ) )
    {
        using BC = BoundaryConditions::ENUMDATA;

        
        EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA boundaryPatch) {
            Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[ boundaryPatch ];

            intType faceEndIndex;
            if ( boundaryPatch == LUT::PositivePatch[ axis ] ) {
                faceEndIndex = mesh.nCells(axis);
            } else {
                faceEndIndex = 0;
            }


            Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[ axis ];
            Axis::ENUMDATA axis2 = LUT::HiOrthogonalAxis[ axis ];
            switch ( boundaryConditions[boundaryPatch].type ) {
                case BC::zeroGradient: 
                {
                    intType k = ( boundaryPatch == LUT::PositivePatch[ axis ] ) ? mesh.nCells(axis)-1 : 0;
                    vertexField.chip(faceEndIndex, axis) = CellPlaneToVertexPlaneFunc( field.chip(k, axis), mesh.interpFactors[axis1], mesh.interpFactors[axis2], mask.chip(k, axis) );
                    break;
                }

                case BC::fixed: 
                {
                    intType k = ( boundaryPatch == LUT::PositivePatch[ axis ] ) ? mesh.nCells(axis)-1 : 0;
                    vertexField.chip(faceEndIndex, axis) = CellPlaneToVertexPlaneFunc( boundaryConditions[boundaryPatch].value, mesh.interpFactors[axis1], mesh.interpFactors[axis2], mask.chip(k, axis)  );
                    break;
                }
                    
                case BC::extrapolated: 
                {
                    intType k_p = ( boundaryPatch == LUT::PositivePatch[ axis ] ) ? mesh.nCells(axis)-1 : 0;
                    Tensor2D offBoundary_p = CellPlaneToVertexPlaneFunc( field.chip(k_p, axis), mesh.interpFactors[axis1], mesh.interpFactors[axis2], mask.chip(k_p, axis)  );

                    intType k_a = ( boundaryPatch == LUT::PositivePatch[ axis ] ) ? k_p-1 : k_p+1;
                    Tensor2D offBoundary_a = CellPlaneToVertexPlaneFunc( field.chip(k_a, axis), mesh.interpFactors[axis1], mesh.interpFactors[axis2], mask.chip(k_p, axis)  );

                    floatType extrapFactor_p = mesh.extrapFactors[boundaryPatch].p;
                    floatType extrapFactor_a = mesh.extrapFactors[boundaryPatch].a;

                    vertexField.chip(faceEndIndex, axis) = offBoundary_p * offBoundary_p.constant( extrapFactor_p )
                                                         + offBoundary_a * offBoundary_a.constant( extrapFactor_a );
                    break;
                }

                case BC::periodic:
                {
                    intType loCellIndex = 0;
                    Tensor2D loCellVertexPlane = CellPlaneToVertexPlaneFunc( field.chip(loCellIndex, axis), mesh.interpFactors[axis1], mesh.interpFactors[axis2], mask.chip(loCellIndex, axis)  );

                    intType hiCellIndex = mesh.nCells(axis) - 1;
                    Tensor2D hiCellVertexPlane = CellPlaneToVertexPlaneFunc( field.chip(hiCellIndex, axis), mesh.interpFactors[axis1], mesh.interpFactors[axis2], mask.chip(hiCellIndex, axis)  );

                    floatType interpFactor = mesh.cellLengths[axis]( loCellIndex ) / ( mesh.cellLengths[axis]( loCellIndex ) + mesh.cellLengths[axis]( hiCellIndex ) );

                    vertexField.chip(faceEndIndex, axis) = loCellVertexPlane * loCellVertexPlane.constant( 1 - interpFactor )
                                                         + hiCellVertexPlane * hiCellVertexPlane.constant( interpFactor );
                    break;
                }
                    
            }


        } );
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




// Assumes ghost cells have been removed
Tensor3D InterpolateToVertex( const Tensor3D &field,
                              const Mesh &mesh, 
                              const Tensor3D &mask,
                              const BoundaryConditionData::Patches &boundaryConditions )
{
    Tensor3D vertexField( field.dimension(0)+1, field.dimension(1)+1, field.dimension(2)+1 );
    vertexField.setZero();

    InterpolateInteriorPoints( vertexField, field, mesh );

    SetBoundaryFaces( vertexField, field, mask, mesh, boundaryConditions, CellPlaneToVertexPlane );

    WeightedAverageEdges( vertexField, mesh );

    WeightedAverageCorners( vertexField, mesh );

    return vertexField;
}



Tensor3D InterpolateToVertexWithMasking( const Tensor3D &field,
                                         const Mesh &mesh, 
                                         const Tensor3D &mask,
                                         const BoundaryConditionData::Patches &boundaryConditions )
{
    Tensor3D vertexField( field.dimension(0)+1, field.dimension(1)+1, field.dimension(2)+1 );
    vertexField.setZero();

    InterpolateInteriorPointsWithMasking( vertexField, field, mask, mesh );

    SetBoundaryFaces( vertexField, field, mask, mesh, boundaryConditions, CellPlaneToVertexPlaneWithMasking );

    WeightedAverageEdges( vertexField, mesh );

    WeightedAverageCorners( vertexField, mesh );

    return vertexField;
}



FieldData<Tensor3D> GetVertexFields( const FieldData<Tensor3D> &fields, 
                                     const Mesh &mesh,
                                     const BoundaryConditionData &bcData,
                                     const Tensor3D &mask )
{
    FieldData<Tensor3D> vertexFields;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        vertexFields.U[axis] = InterpolateToVertex( fields.U[axis], mesh, mask, bcData.fields.U[axis] );
    } );
    vertexFields.P = InterpolateToVertexWithMasking( fields.P, mesh, mask, bcData.fields.P );

    return vertexFields;
}


}   // end namespace FLOW

}   // end namespace CAMIRA