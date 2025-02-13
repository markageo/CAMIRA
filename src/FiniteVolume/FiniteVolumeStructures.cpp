#include "FiniteVolume.h"
#include "../IO/InputProcessing.h"
#include "../Tools/FVTools.h"
#include "../Tools/FVLookups.h"
#include "../Tools/SweepTransformations.h"
#include "../IO/VTKReader.h"

#include <cmath>
#include <stdexcept>

namespace CFD
{


/*-------------------------------------------------------------------------------------*\
                                BoundaryConditionData
\*-------------------------------------------------------------------------------------*/


namespace
{

    // Interpolates user defined profile onto the mesh in 1D. Points that are outside the user defined region are set 
    // to the value of the nearest point.
    Tensor1D InterpProfile1D( const Tensor1D x,
                             const Tensor1D v, 
                             const Tensor1D xquery )
    {
        const intType nValuePoints = x.size(); 
        const intType nQueryPoints = xquery.size(); 
        Tensor1D vquery = Tensor1D( nQueryPoints );

        intType i = 0;
        for ( intType iq = 0; iq != nQueryPoints; iq++ ) {

            bool isBelowMinPoint = ( xquery(iq) <= x(0) );
            if ( isBelowMinPoint ) {
                vquery(iq) = v(0);
                continue;
            }

            bool isAboveMaxPoint = ( xquery(iq) >= x(nValuePoints-1) );
            if ( isAboveMaxPoint ) {
                vquery(iq) = v(nValuePoints-1);
                continue;
            }

            // Internal points
            while ( i != nValuePoints ) {

                bool isBetweenPoints = ( xquery(iq) >= x(i) )  &&  ( xquery(iq) <= x(i+1) );
                if ( isBetweenPoints ) {

                    // Interpolate
                    floatType v0 = v( i ),
                              v1 = v( i + 1 ),
                              x0 = x( i ),
                              x1 = x( i + 1);
                    vquery(iq) = v0 + ( xquery(iq) - x0 ) * ( v1 - v0 ) / ( x1 - x0 );
                    break;
                }
                i++;
            }

        }

        return vquery;
    }



    // Interpolates a 1D profile onto a meshed boundary face
    Tensor2D SetBoundaryProfile1D( const InputData::Profile1D &profile1D,   
                                  const Mesh& mesh,
                                  const BoundaryPatches::ENUMDATA boundaryPatch )
    {
        Axis::ENUMDATA profileAxis  = profile1D.axis,
                    normalAxis   = LUT::BoundaryPatchAxis[ boundaryPatch ],
                    constantAxis = static_cast< Axis::ENUMDATA >( 3 - profileAxis - normalAxis );

        // 1D profile on the mesh points
        Tensor1D interpolatedProfile1D = InterpProfile1D( profile1D.coordinates, profile1D.values, mesh.cellCenters[profileAxis] );

        // Copy out in 2D
        intType nCellsLo = mesh.nCells( LUT::LoOrthogonalAxis[normalAxis] ),
                nCellsHi = mesh.nCells( LUT::HiOrthogonalAxis[normalAxis] );
        Tensor2D boundaryPatchValues( nCellsLo, nCellsHi );
        int constantAxis2D = ( constantAxis == LUT::LoOrthogonalAxis[normalAxis] ) ? 0 : 1; // 3D axis enums cannot be used on the 2D plane
        for ( intType i = 0; i != mesh.nCells(constantAxis); i++ ) {
            boundaryPatchValues.chip(i, constantAxis2D) = interpolatedProfile1D;
        }

        return boundaryPatchValues;
    }



    // Sets a constant value for a boundary face patch
    Tensor2D SetBoundaryProfileConstant( const floatType value,   
                                        const Mesh& mesh,
                                        const BoundaryPatches::ENUMDATA boundaryPatch )
    {
        Axis::ENUMDATA normalAxis   = LUT::BoundaryPatchAxis[ boundaryPatch ];

        intType nCellsLo = mesh.nCells( LUT::LoOrthogonalAxis[normalAxis] ),
                nCellsHi = mesh.nCells( LUT::HiOrthogonalAxis[normalAxis] );

        return Tensor2D( nCellsLo, nCellsHi ).setConstant( value );
    }


}   // end anonymous namespace




BoundaryConditionData SetBoundaryConditionData( const InputData &inputData,
                                                const Mesh &mesh )
{

    BoundaryConditionData bcData;

    // Set boundary condition data for use in solver
    ForAllFieldData( [&] (intType f) {

        EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {

            const auto &bcConfig = inputData.boundaryConditions[f][bp];

            bcData.fields[f][bp].type = bcConfig.type;

            if ( bcConfig.hasUniformValue ) 
                bcData.fields[f][bp].value = SetBoundaryProfileConstant( bcConfig.uniformValue, mesh, bp );

            if ( bcConfig.hasProfile1D ) 
                bcData.fields[f][bp].value = SetBoundaryProfile1D( bcConfig.profile1D, mesh, bp );

        } );

    } );


    // Pressure field will be floating if none of the pathces are fixed
    bcData.pressureFieldIsFloating = true;
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        if ( bcData.fields.P[bp].type == BoundaryConditions::fixed ) {
            bcData.pressureFieldIsFloating = false;
        }
    } );

    return bcData;
}





/*-------------------------------------------------------------------------------------*\
                                    FVCoefficients
\*-------------------------------------------------------------------------------------*/

namespace
{

std::vector< TransportCoefficients::ENUMDATA > PicardEnums( const Axis::ENUMDATA momentumEquation, 
                                                                      const Axis::ENUMDATA velocity)
{
    using enum Axis::ENUMDATA;
    using C = TransportCoefficients::ENUMDATA;

    switch ( momentumEquation ) {
        case X: 

            switch ( velocity ) {
                case X: 
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
                case Y:
                    return {};
                    break;
                case Z:
                    return {};
                    break;
            }
            break;

        case Y: 

            switch ( velocity ) {
                case X: 
                    return {};
                    break;
                case Y:
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
                case Z:
                    return {};
                    break;
            }
            break;

        case Z:

            switch ( velocity ) {
                case X: 
                    return {};
                    break;

                case Y:
                    return {};
                    break;

                case Z:
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
            }
            break;
    }
    return {};
}


std::vector< TransportCoefficients::ENUMDATA > NewtonEnums( const Axis::ENUMDATA momentumEquation, 
                                                            const Axis::ENUMDATA velocity)
{
    using enum Axis::ENUMDATA;
    using C = TransportCoefficients::ENUMDATA;

    switch ( momentumEquation ) {
        case X: 

            switch ( velocity ) {
                case X: 
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
                case Y:
                    return {C::p, C::n, C::s};
                    break;
                case Z:
                    return {C::p, C::t, C::b};
                    break;
            }
            break;

        case Y: 

            switch ( velocity ) {
                case X: 
                    return {C::p, C::e, C::w};
                    break;
                case Y:
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
                case Z:
                    return {C::p, C::t, C::b};;
                    break;
            }
            break;

        case Z:

            switch ( velocity ) {
                case X: 
                    return {C::p, C::e, C::w};
                    break;

                case Y:
                    return {C::p, C::n, C::s};
                    break;

                case Z:
                    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
                    break;
            }
            break;
    }
    return {};
}



std::vector< TransportCoefficients::ENUMDATA > MomentumVelocityEnums( const Axis::ENUMDATA momentumEquation, 
                                                                      const Axis::ENUMDATA velocity,
                                                                      Linearisation li )
{
    switch ( li ) {
        case Linearisation::Picard:
            return PicardEnums( momentumEquation, velocity );

        case Linearisation::Newton:
            return NewtonEnums( momentumEquation, velocity );
    }
    return {};
}                                                                      



std::vector< TransportCoefficients::ENUMDATA > MomentumPressureEnums( const Axis::ENUMDATA momentumEquation )
{
    using enum Axis::ENUMDATA;
    using C = TransportCoefficients::ENUMDATA;

    switch ( momentumEquation ) {
        case X: 
            return { C::p, C::e, C::w };
            break;

        case Y: 
            return { C::p, C::n, C::s };
            break;

        case Z:
            return { C::p, C::t, C::b };
            break;
    }
    return {};
}



std::vector< TransportCoefficients::ENUMDATA > ContinuityPressureEnums( MomentumInterpolation mi ) 
{
    using C = TransportCoefficients::ENUMDATA;
    switch ( mi ) {
        case MomentumInterpolation::Implicit:
            return {C::p, C::n, C::e, C::s, C::w, C::t, C::b, C::nn, C::ee, C::ss, C::ww, C::tt, C::bb};  

        case MomentumInterpolation::SemiExplicit:
            return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
    }
    return {};
}




}   // end anonymous namespace


using C = TransportCoefficients::ENUMDATA;
using enum Axis::ENUMDATA;

// Momentum equations constructor
MomentumEquation::MomentumEquation( const Axis::ENUMDATA axis, 
                                    const iArray3 &dims,
                                    Linearisation li ) :
    AU( { EnumVector<TransportCoefficients, Tensor3D>( MomentumVelocityEnums(axis, X, li), dims + 2*nGhost),
          EnumVector<TransportCoefficients, Tensor3D>( MomentumVelocityEnums(axis, Y, li), dims + 2*nGhost),
          EnumVector<TransportCoefficients, Tensor3D>( MomentumVelocityEnums(axis, Z, li), dims + 2*nGhost) }  ),
    AP( MomentumPressureEnums( axis ), dims( axis ) + 2*nGhost ),
    B( CFD::Tensor3D( dims(X) + 2*nGhost, dims(Y) + 2*nGhost, dims(Z) + 2*nGhost ).setZero() ),
    F( CFD::Tensor3D( dims(X) + 2*nGhost, dims(Y) + 2*nGhost, dims(Z) + 2*nGhost ).setZero() ),
    diagCoeffInv( CFD::Tensor3D( dims(X) + 2*nGhost, dims(Y) + 2*nGhost, dims(Z) + 2*nGhost ).setZero() ),
    diff({ EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::e, C::w}, dims(X) ),
           EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::n, C::s}, dims(Y) ),
           EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::t, C::b}, dims(Z) ) }),
    diffBoundary( 0.0f ),
    BUBoundary(),
    BPBoundary(),   // These should be dimensioned only if needed
    component( axis ),
    linearisation( li )
{};


// Continuity equations constructor
ContinuityEquation::ContinuityEquation( const iArray3 &dims,
                                        MomentumInterpolation mi ) :
    AU( { EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::e, C::w}, dims( X ) + 2*nGhost),
          EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::n, C::s}, dims( Y ) + 2*nGhost),
          EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::t, C::b}, dims( Z ) + 2*nGhost) } ),
    AP( ContinuityPressureEnums( mi ), dims + 2*nGhost ),
    B( Tensor3D( dims(X) + 2*nGhost, dims(Y) + 2*nGhost, dims(Z) + 2*nGhost ).setZero() ),
    F( Tensor3D( dims(X) + 2*nGhost, dims(Y) + 2*nGhost, dims(Z) + 2*nGhost ).setZero() ),
    mwiSparseCoeffs( { std::array<Tensor1D, 4>{ Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero() } ,
                       std::array<Tensor1D, 4>{ Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero() } ,
                       std::array<Tensor1D, 4>{ Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero() } } ),
    mwiCompactCoeffs( { std::array<Tensor1D, 2>{ Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero() } ,
                        std::array<Tensor1D, 2>{ Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero() } ,
                        std::array<Tensor1D, 2>{ Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero() } } ),
    BUBoundary(),
    BPBoundary(),   // These should be dimensioned only if needed
    momentumInterpolation( mi )
{};


// Coefficients class constructor
FVCoefficients::FVCoefficients( const iArray3 &dims, 
                                Linearisation li,
                                MomentumInterpolation mi ) :
    Mom( { MomentumEquation(X, dims, li),  MomentumEquation(Y, dims, li),  MomentumEquation(Z, dims, li) } ),
    Cont( dims, mi ),
    nCells( dims )
{};



/*-------------------------------------------------------------------------------------*\
                                      InitialConditions
\*-------------------------------------------------------------------------------------*/

#ifdef CFD_HAS_VTK_LIB
FieldData<Tensor3D> SetInitialConditionFromVTKFile( const std::string &filename,
                                                    const Mesh &mesh,
                                                    const AxisTransformationMap &axisTransformation,
                                                    const InputData &inputData )
{
    VTK::FieldFileData fieldFileData = VTK::ReadVTKFields( filename );
    
    Mesh inputMesh( fieldFileData.cellFaces, inputData );
    TransformMeshToCodeCoordinates( inputMesh, axisTransformation);
    TransformFieldToCodeCoordinates( fieldFileData.cellFields, axisTransformation);

    TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                  extents = {mesh.nCells(0), mesh.nCells(1), mesh.nCells(2)};

    FieldData<Tensor3D> fields( Tensor3D( mesh.nCells(0) + 2*CFD::nGhost, 
                                          mesh.nCells(1) + 2*CFD::nGhost, 
                                          mesh.nCells(2) + 2*CFD::nGhost).setZero() );

    // Careful! Just checking the mesh is the same size, however it is possible that cell centers are at different locations.
    // Ideally should add the ability to do an interpolation.
    if ( ( inputMesh.nCells(0) != mesh.nCells(0) ) ||
         ( inputMesh.nCells(1) != mesh.nCells(1) ) ||
         ( inputMesh.nCells(2) != mesh.nCells(2) )  ) {
            throw std::runtime_error( "Mesh dimensions for initial condition do not match!" );
    }
    
    ForAllFieldData( [&] (intType f) {
        fields[f].slice( offsets, extents ) = fieldFileData.cellFields[f];
    } );
    
    return fields;
}
#endif



FieldData<Tensor3D> SetInitialConditionUniform( const FieldData<floatType> &constantInitialConditions,
                                                const Mesh &mesh )
{
    TensorIndex3D offsets = {nGhost, nGhost, nGhost},
                  extents = {mesh.nCells(0), mesh.nCells(1), mesh.nCells(2)};

    FieldData<Tensor3D> fields( Tensor3D( mesh.nCells(0) + 2*CFD::nGhost, 
                                          mesh.nCells(1) + 2*CFD::nGhost, 
                                          mesh.nCells(2) + 2*CFD::nGhost).setZero() );

    ForAllFieldData( [&] (intType i) { 
        fields[i].slice( offsets, extents ).setConstant( constantInitialConditions[i] );  
    } );
    
    return fields;
}



FieldData<Tensor3D> InitialiseFields( const Mesh &mesh, 
                                      const InputData &inputData,
                     [[maybe_unused]] const AxisTransformationMap &axisTransformation )
{
    FieldData<Tensor3D> fields;

    #if defined( CFD_HAS_VTK_LIB )
        switch ( inputData.initialConditionType ) {
            case InputData::InitialConditionTypes::uniform:
                fields = SetInitialConditionUniform( inputData.constantInitialConditions, mesh );
                break;

            case InputData::InitialConditionTypes::vtkFile:
                fields = SetInitialConditionFromVTKFile( inputData.initialConditionsFieldFilename, mesh, axisTransformation, inputData );
                break;
        }
    #else
        fields = SetInitialConditionUniform( inputData.constantInitialConditions, mesh );
    #endif

    return fields;
}


}   // end namespace CFD