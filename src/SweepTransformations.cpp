#include "SweepTransformations.h"
#include <Eigen/Geometry>

namespace CFD
{

/*-------------------------------------------------------------------------------------*\
                          AxisTransformationMap Definitions
\*-------------------------------------------------------------------------------------*/

AxisTransformationMap::AxisTransformationMap() :
    m_codeBoundaryPatches( { BoundaryPatches::xPositive,
                             BoundaryPatches::xNegative,
                             BoundaryPatches::yPositive,
                             BoundaryPatches::yNegative,
                             BoundaryPatches::zPositive,
                             BoundaryPatches::zNegative } ),
    m_userBoundaryPatches( m_codeBoundaryPatches ),

    m_codeAxis( { Axis::X, Axis::Y, Axis::Z } ),
    m_userAxis( m_codeAxis )
    {};


void AxisTransformationMap::Set( const BoundaryPatches::ENUMDATA codePatch,
                                 const BoundaryPatches::ENUMDATA userPatch )
{
    A codeAxis = BoundaryPatchAxis[ codePatch ];
    A userAxis = BoundaryPatchAxis[ userPatch ];

    // Set the given patch
    m_codeBoundaryPatches[ codePatch ] = userPatch;
    m_userBoundaryPatches[ userPatch ] = codePatch;

    // Set the opposite patch
    BP codePatchOpposite;
    
    if ( codePatch == PositivePatch[ codeAxis ] ) {
        codePatchOpposite = NegativePatch[ codeAxis ];
    } else {
        codePatchOpposite = PositivePatch[ codeAxis ];
    }

    BP userPatchOpposite;
    if ( userPatch == PositivePatch[ userAxis ] ) {
        userPatchOpposite = NegativePatch[ userAxis ];
    } else {
        userPatchOpposite = PositivePatch[ userAxis ];
    }

    m_codeBoundaryPatches[ codePatchOpposite ] = userPatchOpposite;
    m_userBoundaryPatches[ userPatchOpposite ] = codePatchOpposite;


    // Update the axis
    m_codeAxis[ codeAxis ] = userAxis;
    m_userAxis[ userAxis ] = codeAxis;
}


// Code patch from user patch
const BoundaryPatches::ENUMDATA &AxisTransformationMap::CodePatch(const BP userPatch) const 
{   return m_userBoundaryPatches[ userPatch ]; }

// Code axis from user axis
const Axis::ENUMDATA &AxisTransformationMap::CodeAxis(const A userAxis) const
{ 
    BP userPatch = PositivePatch[ userAxis ];
    BP codePatch = CodePatch( userPatch );
    return BoundaryPatchAxis[ codePatch ];
}

// If code axis is mapped to the negative direction of a user axis
bool AxisTransformationMap::IsCodeAxisReversed( const A codeAxis ) const
{
    BP positiveCodePatch = PositivePatch[ codeAxis ];
    BP mappedUserPatch   = UserPatch( positiveCodePatch );
    A  mappedUserAxis    = BoundaryPatchAxis[ mappedUserPatch ];
    if ( mappedUserPatch == NegativePatch[ mappedUserAxis ] ) {
        return true;
    }
    return false;
}


// User patch from code patch
const BoundaryPatches::ENUMDATA &AxisTransformationMap::UserPatch(const BP codePatch) const
{   return m_codeBoundaryPatches[ codePatch ]; }

// User axis from code axis
const Axis::ENUMDATA &AxisTransformationMap::UserAxis(const A codeAxis) const
{
    BP codePatch = PositivePatch[ codeAxis ];
    BP userPatch = UserPatch( codePatch );
    return BoundaryPatchAxis[ userPatch ];
}

// If user axis is mapped to the negative direction of a code axis
bool AxisTransformationMap::IsUserAxisReversed( const A userAxis ) const
{
    BP positiveUserPatch = PositivePatch[ userAxis ];
    BP mappedCodePatch   = CodePatch( positiveUserPatch );
    A  mappedCodeAxis    = BoundaryPatchAxis[ mappedCodePatch ];
    if ( mappedCodePatch == NegativePatch[ mappedCodeAxis ] ) {
        return true;
    }
    return false;
}



/*-------------------------------------------------------------------------------------*\
                            User -> Code Axis Transformations
\*-------------------------------------------------------------------------------------*/


namespace 
{

    Eigen::Matrix<CFD::intType, 3, 1> Axis2Vector(const CFD::BoundaryPatches::ENUMDATA axis)
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;
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



    CFD::BoundaryPatches::ENUMDATA Vector2Axis(const Eigen::Matrix<CFD::intType, 3, 1> &vector)
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;
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



    // Sets the axis transformation map based on the sweeping directions
    AxisTransformationMap SetAxisTransformation( InputData &inputData) 
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;
        using directionVector = Eigen::Matrix<CFD::intType, 3, 1>;

        AxisTransformationMap axisTransformation;

        // User input for plane sweep direction
        BP planeSweepDirection = inputData.linearSolverSettings.planeSolverSettings.sweepDirection;
        directionVector planeSweepVector = Axis2Vector( planeSweepDirection );

        // User input for line sweep direction
        BP lineSweepDirection = inputData.linearSolverSettings.planeSolverSettings.lineSolverSettings.sweepDirection;
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

        // Change the sweep direction
        inputData.linearSolverSettings.planeSolverSettings.sweepDirection = BP::zPositive;
        inputData.linearSolverSettings.planeSolverSettings.lineSolverSettings.sweepDirection  = BP::yPositive;

        return axisTransformation;
    }



    // Transform an EnumVector of fields data to code coordinates. Only transforms the momentum equations part.
    template< typename T >
    void TransformFieldVectorToCode( EnumVector<Fields, T> &fieldsVector,
                                     const AxisTransformationMap& axisTransformation )
    {

        // Create temporary copy to move data from 
        EnumVector<Fields, T> userFieldsVector = fieldsVector;

        EnumFor<Axis>([&] (Axis::ENUMDATA codeAxis) { // Code axis

            fieldsVector[ AxisVelocity[codeAxis] ] = userFieldsVector[ AxisVelocity[ axisTransformation.UserAxis(codeAxis) ] ];

        } );
    }



    // Remaps the users boundary conditions
    void TransformBoundaryConditions( InputData &inputData, 
                                      const AxisTransformationMap &axisTransformation )
    {
        using BP = CFD::BoundaryPatches::ENUMDATA;
        using F = CFD::Fields::ENUMDATA;
        using A = CFD::Axis::ENUMDATA;

        InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;

        // Temporary for boundary conditions as user specifies them
        InputData::BoundaryConditionData boundaryConditionsUser = boundaryConditions;
        fVector3 domainSizeUser = inputData.domainSize;
        
        EnumFor<Axis>( [&] (A codeAxis) {

            inputData.domainSize( codeAxis ) = domainSizeUser( axisTransformation.UserAxis( codeAxis ) ); 

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
    void TransformMesh(InputData &inputData, 
                       const AxisTransformationMap &axisTransformation )
    {
        using enum Axis::ENUMDATA;
        using enum BoundaryPatches::ENUMDATA;


        // Create temporary copy of mesh data to take data from
        EnumVector<Axis, std::vector<InputData::MeshSegment> > userMeshSegments = inputData.meshSegments;
        fVector3 userDomainSize;

        EnumFor<Axis>( [&] (Axis::ENUMDATA codeAxis) {

            inputData.meshSegments[ codeAxis ] = userMeshSegments[ axisTransformation.UserAxis( codeAxis ) ];
            inputData.domainSize( codeAxis )   = userDomainSize[ axisTransformation.UserAxis( codeAxis ) ];

            if ( axisTransformation.IsCodeAxisReversed( codeAxis ) ) {
                ReverseMesh(inputData.meshSegments[ codeAxis ]);
            }

        } );

    }

    // Remaps the initial conditions
    void TransformInitialConditions( InputData &inputData,
                                     const AxisTransformationMap &axisTransformation )
    {
        TransformFieldVectorToCode( inputData.initialConditions, axisTransformation );
    }


    // Remaps any solver settings that have direction dependence
    void TransformSolver( InputData &inputData,
                          const AxisTransformationMap &axisTransformation )
    {
        TransformFieldVectorToCode( inputData.schemes.implicitRelaxation      , axisTransformation );
        TransformFieldVectorToCode( inputData.linearSolverSettings.relaxation , axisTransformation );
        TransformFieldVectorToCode( inputData.linearSolverSettings.planeSolverSettings.relaxation, axisTransformation );
        TransformFieldVectorToCode( inputData.linearSolverSettings.planeSolverSettings.lineSolverSettings.relaxation, axisTransformation );
    }

}   // end anonymous namespace



AxisTransformationMap TransformUserInputData(InputData &inputData )
{
    
    AxisTransformationMap axisTransformation = SetAxisTransformation( inputData );

    TransformBoundaryConditions( inputData, axisTransformation );
    TransformInitialConditions( inputData, axisTransformation );
    TransformMesh( inputData, axisTransformation );
    TransformSolver( inputData, axisTransformation );

    return axisTransformation;
}





/*-------------------------------------------------------------------------------------*\
                            Code -> User Axis Transformations
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


    // Template template parameter so both EnumVectors and ArrayAllocators can be used
    template< template<typename, typename> typename Container, typename T >
    void TransformFieldVectorToUser( Container<Fields, T> &fieldsVector,
                                     const AxisTransformationMap& axisTransformation )
    {
        static_assert( std::is_same< Container<Fields, T>, EnumVector<Fields, T> >::value    ||
                       std::is_same< Container<Fields, T>, ArrayAllocator<Fields, T> >::value );

        // Create temporary copy to move data from 
        Container<Fields, T> codeFieldsVector = fieldsVector;

        EnumFor<Axis>([&] (Axis::ENUMDATA userAxis) { // User axis

            fieldsVector[ AxisVelocity[userAxis] ] = codeFieldsVector[ AxisVelocity[ axisTransformation.CodeAxis( userAxis ) ] ];

        } );
    }

}   // end anonymous namespace



void TransformToUserCoordinates( Mesh &mesh, 
                                 ArrayAllocator<Fields, array3D> &fields, 
                                 const AxisTransformationMap &axisTransformation )                                  
{
    Eigen::array<intType , Axis::count> shuffleArray;
    Eigen::array<bool, Axis::count> reverseArray;

    // Temporary copy of the mesh to take data from 
    Mesh codeMesh = mesh;

    Axis::ENUMDATA codeAxis;
    bool reverseAxis;

    EnumFor<Axis>( [&] (Axis::ENUMDATA userAxis) {

        codeAxis = axisTransformation.CodeAxis( userAxis );
        reverseAxis = axisTransformation.IsUserAxisReversed( userAxis );

        CopyMeshAxis( mesh, codeMesh, userAxis, codeAxis );

        if ( reverseAxis ) {
            ReverseMeshAxis(mesh, userAxis);
        }

        // Shuffle and reverse arrays used for 3D arrays
        shuffleArray[ userAxis ] = codeAxis;
        reverseArray[ userAxis ] = reverseAxis;

    } );

    // 3D arrays
    EnumFor<Fields>( [&] (Fields::ENUMDATA field) {
        fields[field] = array3D( fields[field] ).shuffle(shuffleArray).reverse(reverseArray);   // Have to make a copy
    } );
    TransformFieldVectorToUser( fields, axisTransformation );

}

}   // end namespace CFD