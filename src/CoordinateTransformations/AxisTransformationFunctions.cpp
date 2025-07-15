#include "AxisTransformationFunctions.h"
#include "../Core/FVLookups.h"
#include <Eigen/Geometry>

namespace CAMIRA
{


/*-------------------------------------------------------------------------------------*\
                                    Helper Functions
\*-------------------------------------------------------------------------------------*/

namespace  
{

    // Returns a copy of a 1D array that has been reversed 
    Tensor1D ReversedArray1D( Tensor1D &array )
    {
        Eigen::array<bool, 1> rev({true});
        return array.reverse( rev );
    };

}   // end anonymous namespace





/*-------------------------------------------------------------------------------------*\
                            Input Data Transformations
\*-------------------------------------------------------------------------------------*/


namespace 
{

    Eigen::Matrix<CAMIRA::intType, 3, 1> Axis2Vector(const CAMIRA::BoundaryPatches::ENUMDATA axis)
    {
        using BP = CAMIRA::BoundaryPatches::ENUMDATA;
        switch (axis)
        {
            case (BP::xPositive):
                return {1, 0, 0};
                
            case (BP::xNegative):
                return {-1, 0, 0};

            case (BP::yPositive):
                return {0, 1, 0};

            case (BP::yNegative):
                return {0, -1, 0};

            case (BP::zPositive):
                return {0, 0, 1};

            case (BP::zNegative):
                return {0, 0, -1};

            default:
                // throw ERROR - invalud axis enum
                return {0, 0, 0};
        }
    }



    CAMIRA::BoundaryPatches::ENUMDATA Vector2Axis(const Eigen::Matrix<CAMIRA::intType, 3, 1> &vector)
    {
        using BP = CAMIRA::BoundaryPatches::ENUMDATA;
        if        ( ( vector.array() ==  Axis2Vector(BP::xPositive).array() ).all() ) {
            return BP::xPositive;

        } else if ( ( vector.array() ==  Axis2Vector(BP::xNegative).array() ).all() ) {
            return BP::xNegative;

        } else if ( ( vector.array() ==  Axis2Vector(BP::yPositive).array() ).all() ) {
            return BP::yPositive;

        } else if ( ( vector.array() ==  Axis2Vector(BP::yNegative).array() ).all() ) {
            return BP::yNegative;

        } else if ( ( vector.array() ==  Axis2Vector(BP::zPositive).array() ).all() ) {
            return BP::zPositive;

        } else if ( ( vector.array() ==  Axis2Vector(BP::zNegative).array() ).all() ) {
            return BP::zNegative;

        } else {
            // throw ERROR - invalid unit vector
            return BP::xPositive;
        }
    }



    // Transform user FieldData struct to code coordinates. Only transforms the momentum equations part.
    template< typename T >
    void TransformFieldDataToCode( FieldData<T> &fieldData,
                                   const AxisTransformationMap& axisTransformation )
    {
        // Create temporary copy to move data from 
        FieldData<T> userFieldData = fieldData;

        EnumFor<Axis>([&] (Axis::ENUMDATA codeAxis) { // Code axis

            fieldData.U[ codeAxis ] = userFieldData.U[ axisTransformation.UserAxis(codeAxis) ];

        } );
    }



    // Remaps the users boundary conditions
    void TransformBoundaryConditions( InputData &inputData, 
                                      const AxisTransformationMap &axisTransformation )
    {
        auto &boundaryConditions = inputData.boundaryConditions;

        // Temporary for boundary conditions as user specifies them
        const auto boundaryConditionsUser = boundaryConditions;
        fArray3 domainLowerBoundsUser = inputData.domainLowerBounds,
                domainUpperBoundsUser = inputData.domainUpperBounds;
        
        EnumFor<Axis>( [&] (Axis::ENUMDATA codeAxis) {

            inputData.domainLowerBounds( codeAxis ) = domainLowerBoundsUser( axisTransformation.UserAxis( codeAxis ) ); 
            inputData.domainUpperBounds( codeAxis ) = domainUpperBoundsUser( axisTransformation.UserAxis( codeAxis ) ); 

        } );

            
        // Transform boundary condition directions
        EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA patch) { 

            ForAllFieldData( [&] ( intType f) {

                boundaryConditions[f][patch] = boundaryConditionsUser[f][ axisTransformation.UserPatch( patch ) ];

                // Boundary profiles need to be transformed
                if ( boundaryConditions[f][patch].hasProfile1D ) {
                    InputData::Profile1D &profile1d = boundaryConditions[f][patch].profile1D;
                    const InputData::Profile1D profile1dUser = profile1d;

                    profile1d.axis = axisTransformation.CodeAxis( profile1dUser.axis );

                    if ( axisTransformation.UserAxisReversed( profile1dUser.axis ) ) {
                        profile1d.coordinates = - ReversedArray1D( profile1d.coordinates );
                        profile1d.values      = - ReversedArray1D( profile1d.values );
                    }

                }

            } );

        } );


        // Now the momentum directions
        TransformFieldDataToCode( boundaryConditions, axisTransformation );
    }



    // Reverse a mesh in a given axis
    void ReverseMesh(std::vector< InputData::MeshSegment > &meshSegments)
    {
        // Reverse the order of the segments
        std::reverse(meshSegments.begin(), meshSegments.end());

        // Now flip each segment
        for (auto &segment : meshSegments) {
            std::swap(segment.endCoordinate, segment.startCoordinate);
            segment.endCoordinate = - segment.endCoordinate;
            segment.startCoordinate = - segment.startCoordinate;
            segment.biasFactor = - segment.biasFactor;
        }
    }



    // Remaps the user mesh
    void TransformMesh(InputData &inputData, 
                       const AxisTransformationMap &axisTransformation )
    {
        using enum Axis::ENUMDATA;
        using enum BoundaryPatches::ENUMDATA;


        // Create temporary copy of mesh data to take data from
        EnumVector<Axis, std::vector<InputData::MeshSegment> > userMeshSegments = inputData.meshSegments;
        fArray3 userDomainLowerBounds = inputData.domainLowerBounds,
                userDomainUpperBounds = inputData.domainUpperBounds;

        EnumFor<Axis>( [&] (Axis::ENUMDATA codeAxis) {

            inputData.meshSegments[ codeAxis ] = userMeshSegments[ axisTransformation.UserAxis( codeAxis ) ];
            inputData.domainLowerBounds( codeAxis )   = userDomainLowerBounds[ axisTransformation.UserAxis( codeAxis ) ];
            inputData.domainUpperBounds( codeAxis )   = userDomainUpperBounds[ axisTransformation.UserAxis( codeAxis ) ];

            if ( axisTransformation.CodeAxisReversed( codeAxis ) ) {
                ReverseMesh(inputData.meshSegments[ codeAxis ]);
            }

        } );

    }



    template< typename AxisArray3 >
    AxisArray3 TransformAxisArray3ToCode( const AxisArray3 &userArray, 
                                          const AxisTransformationMap &axisTransformation )
    {
        static_assert( std::is_same< AxisArray3, fArray3 >::value ||
                       std::is_same< AxisArray3, iArray3 >::value );

        AxisArray3 codeArray = userArray;    
        EnumFor<Axis>([&] (Axis::ENUMDATA codeAxis) { 

            codeArray[ codeAxis ] = userArray[ axisTransformation.UserAxis( codeAxis ) ];

        } );

        return codeArray;
    }


    // Used for arrays that represent coordinates in the mesh
    template< typename AxisArray3 >
    AxisArray3 TransformPositionArray3ToCode( const AxisArray3 &userArray, 
                                              const AxisTransformationMap &axisTransformation )
    {
        static_assert( std::is_same< AxisArray3, fArray3 >::value ||
                       std::is_same< AxisArray3, fVector3 >::value );

        AxisArray3 codeArray = userArray;    
        EnumFor<Axis>([&] (Axis::ENUMDATA codeAxis) { 

            codeArray[ codeAxis ] = userArray[ axisTransformation.UserAxis( codeAxis ) ];

            bool axisReversed = axisTransformation.CodeAxisReversed( codeAxis );
            if ( axisReversed ) {
                codeArray[ codeAxis ] *= -1.0f;
            }

        } );

        return codeArray;
    }



    // Remaps the geometry data
    void TransformGeometry( InputData &inputData, 
                            const AxisTransformationMap &axisTransformation )
    {
        using enum Axis::ENUMDATA;
        using enum BoundaryPatches::ENUMDATA;

        // Blocks
        for ( InputData::SolidBlockData &solidBlock : inputData.solidBlocks ) {
            solidBlock.centerPosition = TransformPositionArray3ToCode( solidBlock.centerPosition, axisTransformation );
            solidBlock.dimensions     = TransformAxisArray3ToCode( solidBlock.dimensions    , axisTransformation );
            solidBlock.rotation       = TransformAxisArray3ToCode( solidBlock.rotation      , axisTransformation );
        }

        // Spheres
        for ( InputData::SolidSphereData &solidSphere : inputData.solidSpheres ) {
            solidSphere.centerPosition = TransformPositionArray3ToCode( solidSphere.centerPosition, axisTransformation );
        }
    }



    // Remaps the initial conditions
    void TransformInitialConditions( InputData &inputData,
                                     const AxisTransformationMap &axisTransformation )
    {
        TransformFieldDataToCode( inputData.constantInitialConditions, axisTransformation );
    }


    // Remaps any solver settings that have direction dependence
    void TransformSolver( InputData &inputData,
                          const AxisTransformationMap &axisTransformation )
    {
        TransformFieldDataToCode( inputData.smootherSettings.relaxation , axisTransformation );
    }


    // Remaps any output settings that have direction dependence
    void TransformOutput( InputData &inputData,
                          const AxisTransformationMap &axisTransformation )
    {
        // Probe locations
        for ( auto &probe : inputData.probes ) {
            probe.location = TransformPositionArray3ToCode( probe.location, axisTransformation );
        }
    }

}   // end anonymous namespace




// Sets the axis transformation map based on the sweeping directions
AxisTransformationMap CreateAxisTransformation( BoundaryPatches::ENUMDATA planeSweepDirection,
                                                BoundaryPatches::ENUMDATA lineSweepDirection ) 
{
    using BP = CAMIRA::BoundaryPatches::ENUMDATA;
    using directionVector = Eigen::Matrix<CAMIRA::intType, 3, 1>;

    AxisTransformationMap axisTransformation;

    // User input for plane sweep direction
    directionVector planeSweepVector = Axis2Vector( planeSweepDirection );

    // User input for line sweep direction
    directionVector lineSweepVector = Axis2Vector( lineSweepDirection );

    // Plane and line sweep direction must be orthogonal
    if ( abs( planeSweepVector.dot( lineSweepVector ) ) == 1 ) {
        // throw ERROR - plane and line sweep directions cannot be the same
    }

    // Point sweep direction, chosen to make right handed coordinate system
    directionVector pointSweepVector = lineSweepVector.cross( planeSweepVector );
    BP pointSweepDirection = Vector2Axis( pointSweepVector );

    // Update the axisTransformation map
    axisTransformation.Set( BP::xPositive, pointSweepDirection);
    axisTransformation.Set( BP::xNegative, Vector2Axis( -pointSweepVector ));

    axisTransformation.Set( BP::yPositive, lineSweepDirection);
    axisTransformation.Set( BP::yNegative, Vector2Axis( -lineSweepVector ));

    axisTransformation.Set( BP::zPositive, planeSweepDirection);
    axisTransformation.Set( BP::zNegative, Vector2Axis( -planeSweepVector ));

    return axisTransformation;
}




void TransformUserInputData( InputData &inputData,
                             const AxisTransformationMap &axisTransformation )
{
    inputData.smootherSettings.planeSweepDirection = BoundaryPatches::zPositive;
    inputData.smootherSettings.lineSweepDirection  = BoundaryPatches::yPositive;

    TransformBoundaryConditions( inputData, axisTransformation );
    TransformInitialConditions( inputData, axisTransformation );
    TransformGeometry( inputData, axisTransformation );
    TransformMesh( inputData, axisTransformation );
    TransformSolver( inputData, axisTransformation );
    TransformOutput( inputData, axisTransformation );
}





/*-------------------------------------------------------------------------------------*\
                            Solver Structures Transformations
\*-------------------------------------------------------------------------------------*/


namespace
{

    // Copies the mesh data from the axis of one mesh to the specified axis of another mesh.
    void CopyMeshAxis( Mesh &targetMesh,
                        const Mesh &sourceMesh,
                        const Axis::ENUMDATA &targetAxis,
                        const Axis::ENUMDATA &sourceAxis )
    {
        targetMesh.nCells( targetAxis )            = sourceMesh.nCells( sourceAxis );
        targetMesh.cellCenters[ targetAxis ]       = sourceMesh.cellCenters[ sourceAxis ];
        targetMesh.cellFaces[ targetAxis ]         = sourceMesh.cellFaces[ sourceAxis ];
        targetMesh.cellLengths[ targetAxis ]       = sourceMesh.cellLengths[ sourceAxis ];
        targetMesh.cellLengthsInv[ targetAxis ]    = sourceMesh.cellLengthsInv[ sourceAxis ];
        targetMesh.cellCenterDiffInv[ targetAxis ] = sourceMesh.cellCenterDiffInv[ sourceAxis ];
        targetMesh.interpFactors[ targetAxis ]     = sourceMesh.interpFactors[ sourceAxis ];

        targetMesh.extrapFactors[ LUT::PositivePatch[ targetAxis ] ] = sourceMesh.extrapFactors[ LUT::PositivePatch[ sourceAxis ] ];
        targetMesh.extrapFactors[ LUT::NegativePatch[ targetAxis ] ] = sourceMesh.extrapFactors[ LUT::NegativePatch[ sourceAxis ] ];
    }


    // Reverses a mesh along a coordinate direction
    void ReverseMeshAxis(Mesh &mesh, const Axis::ENUMDATA axis)
    {
        // Eigen .reverse is done in place, so need to made a copy
        mesh.cellCenters[axis]       = - ReversedArray1D( mesh.cellCenters[axis] );
        mesh.cellFaces[axis]         = - ReversedArray1D( mesh.cellFaces[axis] );
        mesh.cellLengths[axis]       = ReversedArray1D( mesh.cellLengths[axis] );
        mesh.cellLengthsInv[axis]    = ReversedArray1D( mesh.cellLengthsInv[axis] );
        mesh.cellCenterDiffInv[axis] = ReversedArray1D( mesh.cellCenterDiffInv[axis] );
        mesh.interpFactors[axis]     = 1 - ReversedArray1D( mesh.interpFactors[axis] );

        std::swap( mesh.extrapFactors[ LUT::PositivePatch[axis] ], mesh.extrapFactors[ LUT::NegativePatch[axis] ] );
    }



    template< typename T >
    void TransformAxisEnumVectorToUser( EnumVector<Axis, T> &axisVector,
                                        const AxisTransformationMap& axisTransformation )
    {
        // Create temporary copy to move data from 
        EnumVector<Axis, T> codeAxisVector = axisVector;

        EnumFor<Axis>([&] (Axis::ENUMDATA userAxis) { 

            axisVector[ userAxis ] = codeAxisVector[ axisTransformation.CodeAxis( userAxis ) ];

        } );
    }



    template< typename T >
    void TransformAxisEnumVectorToCode( EnumVector<Axis, T> &axisVector,
                                        const AxisTransformationMap& axisTransformation )
    {
        // Create temporary copy to move data from 
        EnumVector<Axis, T> userAxisVector = axisVector;

        EnumFor<Axis>([&] (Axis::ENUMDATA codeAxis) { 

            axisVector[ codeAxis ] = userAxisVector[ axisTransformation.UserAxis( codeAxis ) ];

        } );
    }



    template< typename T >
    void TransformBoundaryPatchEnumVectorToUser( EnumVector<BoundaryPatches, T> &boundaryPatchVector,
                                                 const AxisTransformationMap& axisTransformation )
    {
        // Create temporary copy to move data from 
        EnumVector<BoundaryPatches, T> codeBoundaryPatchVector = boundaryPatchVector;

        EnumFor<BoundaryPatches>([&] (BoundaryPatches::ENUMDATA userBP) { 

            boundaryPatchVector[ userBP ] = codeBoundaryPatchVector[ axisTransformation.CodePatch( userBP ) ];

        } );

    }

}   // end anonymous namespace



void TransformMeshToUserCoordinates( Mesh &mesh,
                                     const AxisTransformationMap &axisTransformation)
{
    // Temporary copy of the mesh to take data from 
    Mesh codeMesh = mesh;

    EnumFor<Axis>( [&] (Axis::ENUMDATA userAxis) {

        Axis::ENUMDATA codeAxis = axisTransformation.CodeAxis( userAxis );
        bool reverseAxis = axisTransformation.UserAxisReversed( userAxis );

        CopyMeshAxis( mesh, codeMesh, userAxis, codeAxis );

        if ( reverseAxis ) {
            ReverseMeshAxis(mesh, userAxis);
        }

    } );

}

void TransformScalarFieldToUserCoordinates( Tensor3D &scalarField,
                                            const AxisTransformationMap &axisTransformation )
{
    Eigen::array<intType , Axis::count> shuffleArray;
    Eigen::array<bool, Axis::count> reverseArray;

    EnumFor<Axis>( [&] (Axis::ENUMDATA userAxis) {

        Axis::ENUMDATA codeAxis = axisTransformation.CodeAxis( userAxis );
        bool reverseAxis = axisTransformation.UserAxisReversed( userAxis );

        // Shuffle and reverse arrays used for 3D arrays
        shuffleArray[ userAxis ] = codeAxis;
        reverseArray[ userAxis ] = reverseAxis;

    } );

    scalarField = Tensor3D( scalarField ).shuffle(shuffleArray).reverse(reverseArray);   // Have to make a copy
}


void TransformVectorFieldToUserCoordinates( EnumVector<Axis, Tensor3D> &vectorField,
                                            const AxisTransformationMap &axisTransformation )
{
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        TransformScalarFieldToUserCoordinates( vectorField[axis], axisTransformation );
    } );
    TransformAxisEnumVectorToUser( vectorField, axisTransformation );
}



void TransformBCDataToUserCoordinates( BoundaryConditionData &bcData,
                                       const AxisTransformationMap &axisTransformation )
{
    
    // Transform the EnumVector
    TransformAxisEnumVectorToUser( bcData.fields.U, axisTransformation );
    ForAllFieldData( [&] (intType f) {
        TransformBoundaryPatchEnumVectorToUser( bcData.fields[f], axisTransformation );
    } );


    // Transform the fixed boundary array
    EnumFor<Axis>( [&] (Axis::ENUMDATA userAxis) {

        Axis::ENUMDATA userAxis1 = LUT::LoOrthogonalAxis[ userAxis ],
                       userAxis2 = LUT::HiOrthogonalAxis[ userAxis ];

        Axis::ENUMDATA userTransformedAxis1 = axisTransformation.CodeAxis( userAxis1 ),
                       userTransformedAxis2 = axisTransformation.CodeAxis( userAxis2 );

        Axis::ENUMDATA codeAxis  = axisTransformation.CodeAxis( userAxis ),
                       codeAxis1 = LUT::LoOrthogonalAxis[ codeAxis ],
                       codeAxis2 = LUT::HiOrthogonalAxis[ codeAxis ];

        bool reverseAxis1 = axisTransformation.UserAxisReversed( userAxis1 ),
             reverseAxis2 = axisTransformation.UserAxisReversed( userAxis2 );

        Eigen::array<int , Axis::count-1> shuffleArray = {0, 1};
        Eigen::array<bool, Axis::count-1> reverseArray = {false, false};

        bool shouldTranspose =  ( codeAxis1 != userTransformedAxis1 ) || ( codeAxis2 != userTransformedAxis2 );
        if ( shouldTranspose )
            shuffleArray = {1, 0};

        if ( reverseAxis1 )
            reverseArray[ 0 ] = true;

        if ( reverseAxis2 )
            reverseArray[ 1 ] = true;

        ForAllFieldData( [&] (intType f) {
            bcData.fields[f][ LUT::PositivePatch[userAxis] ].value = Tensor2D( bcData.fields[f][ LUT::PositivePatch[userAxis] ].value.shuffle(shuffleArray).reverse(reverseArray) );
            bcData.fields[f][ LUT::NegativePatch[userAxis] ].value = Tensor2D( bcData.fields[f][ LUT::NegativePatch[userAxis] ].value.shuffle(shuffleArray).reverse(reverseArray) );
        } );

    } );

}


void TransformMeshToCodeCoordinates( Mesh &mesh,
                                     const AxisTransformationMap &axisTransformation )
{

    // Temporary copy of the mesh to take data from 
    Mesh userMesh = mesh;

    EnumFor<Axis>( [&] (Axis::ENUMDATA codeAxis) {

        Axis::ENUMDATA userAxis = axisTransformation.UserAxis( codeAxis );
        bool reverseAxis = axisTransformation.CodeAxisReversed( codeAxis );

        CopyMeshAxis( mesh, userMesh, codeAxis, userAxis );

        if ( reverseAxis ) {
            ReverseMeshAxis(mesh, codeAxis);
        }

    } );

}


void TransformFieldToCodeCoordinates( FieldData<Tensor3D> &fieldData, 
                                      const AxisTransformationMap &axisTransformation)
{

    Eigen::array<intType , Axis::count> shuffleArray;
    Eigen::array<bool, Axis::count> reverseArray;

    EnumFor<Axis>( [&] (Axis::ENUMDATA codeAxis) {

        Axis::ENUMDATA userAxis = axisTransformation.UserAxis( codeAxis );
        bool reverseAxis = axisTransformation.CodeAxisReversed( codeAxis );

        // Shuffle and reverse arrays used for 3D arrays
        shuffleArray[ codeAxis ] = userAxis;
        reverseArray[ codeAxis ] = reverseAxis;

    } );

    // 3D arrays
    ForAllFieldData( [&] (intType f) {
        fieldData[f] = Tensor3D( fieldData[f] ).shuffle(shuffleArray).reverse(reverseArray);   // Have to make a copy
    } );
    TransformAxisEnumVectorToCode( fieldData.U, axisTransformation );

}


}   // end namespace CAMIRA