#include "FiniteVolume.h"
#include "../Tools/FVLookups.h"

namespace CFD
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


    Tensor2D CellPlaneToVertexPlane( const Tensor2D &cellPlane, 
                                    const Tensor1D &interpFactors0,
                                    const Tensor1D &interpFactors1 )
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



    // Set the face values of the cell verticies using the boundary conditions. Corners/edges will be incorrect, these will be set individually later
    void SetBoundaryFaces( Tensor3D &vertexField,
                           const Tensor3D &field,
                           const Mesh &mesh,
                           const BoundaryConditionData::Patches &boundaryConditions )
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
                    vertexField.chip(faceEndIndex, axis) = CellPlaneToVertexPlane( field.chip(k, axis), mesh.interpFactors[axis1], mesh.interpFactors[axis2] );
                    break;
                }

                case BC::fixed: 
                {
                    vertexField.chip(faceEndIndex, axis) = CellPlaneToVertexPlane( boundaryConditions[boundaryPatch].value, mesh.interpFactors[axis1], mesh.interpFactors[axis2] );
                    break;
                }
                    
                case BC::extrapolated: 
                {
                    intType k_p = ( boundaryPatch == LUT::PositivePatch[ axis ] ) ? mesh.nCells(axis)-1 : 0;
                    Tensor2D offBoundary_p = CellPlaneToVertexPlane( field.chip(k_p, axis), mesh.interpFactors[axis1], mesh.interpFactors[axis2] );

                    intType k_a = ( boundaryPatch == LUT::PositivePatch[ axis ] ) ? k_p-1 : k_p+1;
                    Tensor2D offBoundary_a = CellPlaneToVertexPlane( field.chip(k_a, axis), mesh.interpFactors[axis1], mesh.interpFactors[axis2] );;

                    floatType extrapFactor_p = mesh.extrapFactors[boundaryPatch].p;
                    floatType extrapFactor_a = mesh.extrapFactors[boundaryPatch].a;

                    vertexField.chip(faceEndIndex, axis) = offBoundary_p * offBoundary_p.constant( extrapFactor_p )
                                                         + offBoundary_a * offBoundary_a.constant( extrapFactor_a );
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
                              const BoundaryConditionData::Patches &boundaryConditions )
{
    Tensor3D vertexField( field.dimension(0)+1, field.dimension(1)+1, field.dimension(2)+1 );
    vertexField.setZero();

    InterpolateInteriorPoints( vertexField, field, mesh );

    SetBoundaryFaces( vertexField, field, mesh, boundaryConditions );

    WeightedAverageEdges( vertexField, mesh );

    WeightedAverageCorners( vertexField, mesh );

    return vertexField;
}



FieldData<Tensor3D> GetVertexFields( const FieldData<Tensor3D> &fields, 
                                    const Mesh &mesh,
                                    const BoundaryConditionData &bcData )
{
    FieldData<Tensor3D> vertexFields;
    ForAllFieldData( [&] (intType f) {
        vertexFields[f] = InterpolateToVertex( fields[f], mesh, bcData.fields[f] );
    } );
    return vertexFields;
}




} // end namespace CFD