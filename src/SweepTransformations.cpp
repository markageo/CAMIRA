#include "SweepTransformations.h"

namespace CFD
{

namespace
{

    // Transform an EnumVector of fields data to code coordinates. Only transforms the momentum equations part.
    template< typename T >
    void TransformFieldVectorToCode( EnumVector<Fields, T> &fieldsVector,
                                        const InputData::AxisTransformationMap& axisTransformation )
    {
        using F = Fields::ENUMDATA;
        EnumVector<Axis, F> axisField({ F::U, F::V, F::W });

        // Create temporary copy to move data from 
        EnumVector<Fields, T> userFieldsVector = fieldsVector;

        BoundaryPatches::ENUMDATA userPatch;
        Axis::ENUMDATA userAxis;

        EnumFor<Axis>([&] (Axis::ENUMDATA axis) { // Code axis

            userPatch = axisTransformation.UserPatch( PositivePatch[ axis ] );
            userAxis = BoundaryPatchAxis[ userPatch ];

            fieldsVector[ axisField[axis] ] = userFieldsVector[ axisField[userAxis] ];

        } );
    }



    // Remaps the users boundary conditions
    void TransformBoundaryConditions(InputData &inputData)
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;
        using F = CFD::Fields::ENUMDATA;
        using A = CFD::Axis::ENUMDATA;

        const InputData::AxisTransformationMap &axisTransformation = inputData.axisTransformation;
        InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;

        // Temporary for boundary conditions as user specifies them
        InputData::BoundaryConditionData boundaryConditionsUser = boundaryConditions;
        floatVector3 domainSizeUser = inputData.domainSize;
        

        BoundaryPatches::ENUMDATA userPatch;
        Axis::ENUMDATA userAxis;
        EnumFor<Axis>( [&] (A axis) {

            userPatch = axisTransformation.UserPatch( PositivePatch[ axis ] );
            userAxis = BoundaryPatchAxis[ userPatch ];

            inputData.domainSize( axis ) = domainSizeUser( userAxis ); 

        } );


        // Transform just the boundary conditions
        EnumFor<Fields>( [&] (F field) {
            
            EnumFor<BoundaryPatches>( [&] (BP patch) { 

                boundaryConditions[field][patch] = boundaryConditionsUser[field][ axisTransformation.UserPatch( patch ) ];

            } );

        } );

        // Now the momentum directions
        TransformFieldVectorToCode( boundaryConditions, axisTransformation );
    }



    // Reverse a mesh in a given axis
    void ReverseMesh(std::vector< InputData::MeshSegment > &meshSegments)
    {
        // Reverse the order of the segments
        std::reverse(meshSegments.begin(), meshSegments.end());

        // Now flip each segment
        for (auto &segment : meshSegments) {
            std::swap(segment.upperBound, segment.lowerBound);
            segment.upperBound = - segment.upperBound;
            segment.lowerBound = - segment.lowerBound;
            segment.biasFactor = - segment.biasFactor;
        }
    }



    // Remaps the user mesh
    void TransformMesh(InputData &inputData)
    {
        using enum Axis::ENUMDATA;
        using enum BoundaryPatches::ENUMDATA;


        // Create temporary copy of mesh data to take data from
        EnumVector<Axis, std::vector<InputData::MeshSegment> > userMeshSegments = inputData.meshSegments;
        floatVector3 userDomainSize;

        BoundaryPatches::ENUMDATA userPatch;
        Axis::ENUMDATA userAxis;

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

            userPatch = inputData.axisTransformation.UserPatch( PositivePatch[ axis ] );
            userAxis = BoundaryPatchAxis[ userPatch ];

            inputData.meshSegments[ axis ] = userMeshSegments[ userAxis ];
            inputData.domainSize( axis )   = userDomainSize[ userAxis ];

            if ( userPatch == NegativePatch[ userAxis ] ) {
                ReverseMesh(inputData.meshSegments[ axis ]);
            }

        } );

    }

    // Remaps the initial conditions
    void TransformInitialConditions( InputData &inputData )
    {
        TransformFieldVectorToCode( inputData.initialConditions, inputData.axisTransformation );
    }


    // Remaps any solver settings that have direction dependence
    void TransformSolver( InputData &inputData )
    {
        TransformFieldVectorToCode( inputData.schemes.implicitRelaxation      , inputData.axisTransformation );
        TransformFieldVectorToCode( inputData.linearSolverSettings.relaxation , inputData.axisTransformation );
        TransformFieldVectorToCode( inputData.planeSolverSettings.relaxation  , inputData.axisTransformation );
        TransformFieldVectorToCode( inputData.lineSolverSettings.relaxation   , inputData.axisTransformation );
    }

}   // end anonymous namespace



InputData TransformUserInputData(const InputData &userInputData )
{
    InputData codeInputData = userInputData;

    TransformBoundaryConditions(codeInputData);
    TransformInitialConditions(codeInputData);
    TransformMesh(codeInputData);
    TransformSolver(codeInputData);

    return codeInputData;
}



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

        targetMesh.extrapFactors[ PositivePatch[ targetAxis ] ] = sourceMesh.extrapFactors[ PositivePatch[ sourceAxis ] ];
        targetMesh.extrapFactors[ NegativePatch[ targetAxis ] ] = sourceMesh.extrapFactors[ NegativePatch[ sourceAxis ] ];
    }


    // Returns a copy of a 1D array that has been reversed 
    array1D ReversedArray1D( array1D &array )
    {
        Eigen::array<bool, 1> rev({true});
        return array.reverse( rev );
    };


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

        std::swap( mesh.extrapFactors[ PositivePatch[axis] ], mesh.extrapFactors[ NegativePatch[axis] ] );
    }


    // Transform an EnumVector of fields data to user coordinates. Only transforms the momentum equations part.
    template< typename T >
    void TransformFieldVectorToUser( EnumVector<Fields, T> &fieldsVector,
                                    const InputData::AxisTransformationMap& axisTransformation )
    {
        using F = Fields::ENUMDATA;
        EnumVector<Axis, F> axisField({ F::U, F::V, F::W });

        // Create temporary copy to move data from 
        EnumVector<Fields, T> codeFieldsVector = fieldsVector;

        BoundaryPatches::ENUMDATA codepatch;
        Axis::ENUMDATA codeAxis;

        EnumFor<Axis>([&] (Axis::ENUMDATA axis) { // User axis

            codepatch = axisTransformation.CodePatch( PositivePatch[ axis ] );
            codeAxis = BoundaryPatchAxis[ codepatch ];

            fieldsVector[ axisField[axis] ] = codeFieldsVector[ axisField[codeAxis] ];

        } );
    }

    template< typename T >
    void TransformFieldVectorToUser( ArrayAllocator<Fields, T> &fieldsVector,
                                    const InputData::AxisTransformationMap& axisTransformation )
    {
        using F = Fields::ENUMDATA;
        EnumVector<Axis, F> axisField({ F::U, F::V, F::W });

        // Create temporary copy to move data from 
        ArrayAllocator<Fields, T> codeFieldsVector = fieldsVector;

        BoundaryPatches::ENUMDATA codepatch;
        Axis::ENUMDATA codeAxis;

        EnumFor<Axis>([&] (Axis::ENUMDATA axis) { // User axis

            codepatch = axisTransformation.CodePatch( PositivePatch[ axis ] );
            codeAxis = BoundaryPatchAxis[ codepatch ];

            fieldsVector[ axisField[axis] ] = codeFieldsVector[ axisField[codeAxis] ];

        } );
    }

}   // end anonymous namespace



void TransformToUserCoordinates( Mesh &mesh, 
                                 ArrayAllocator<Fields, array3D> &fields, 
                                 const InputData::AxisTransformationMap &axisTransformation )                                  
{
    Eigen::array<intType , Axis::count> shuffleArray;
    Eigen::array<bool, Axis::count> reverseArray;

    // Temporary copy of the mesh to take data from 
    Mesh codeMesh = mesh;

    BoundaryPatches::ENUMDATA codePatch;
    Axis::ENUMDATA codeAxis;
    bool reverseAxis;

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        codePatch = axisTransformation.CodePatch( PositivePatch[ axis ] );
        codeAxis = BoundaryPatchAxis[ codePatch ];
        reverseAxis = ( codePatch == NegativePatch[ codeAxis ] );

        CopyMeshAxis( mesh, codeMesh, axis, codeAxis );

        if ( reverseAxis ) {
            ReverseMeshAxis(mesh, axis);
        }

        // Shuffle and reverse arrays used for 3D arrays
        shuffleArray[ axis ] = codeAxis;
        reverseArray[ axis ] = reverseAxis;

    } );

    // 3D arrays
    EnumFor<Fields>( [&] (Fields::ENUMDATA field) {
        fields[field] = array3D( fields[field] ).shuffle(shuffleArray).reverse(reverseArray);   // Have to make a copy
    } );
    TransformFieldVectorToUser( fields, axisTransformation );

}

}   // end namespace CFD