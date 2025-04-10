#include "FiniteVolumeStructures.h"
#include "FiniteVolume.h"
#include "../IO/InputProcessing.h"
#include "../Core/FVTools.h"
#include "../Core/FVLookups.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"
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

std::vector< TransportCoefficients::ENUMDATA > MomentumVelocityEnums()
{
    using enum Axis::ENUMDATA;
    using C = TransportCoefficients::ENUMDATA;
    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
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

// Default constructor must be defined explicity since class holds references
FVCoefficients::FVCoefficients() :
    Cont( { m_gradientCoeffs,
            CFD::EnumVector<TransportCoefficients, Tensor3D>(),
            CFD::Tensor3D(),
            CFD::Tensor3D() } ),

    Mom( { MomentumEquation({ X,
                            m_momentumVelocityCoeffs, 
                            m_gradientCoeffs[X],
                            CFD::Tensor3D(),
                            CFD::Tensor3D() }),
                            
         MomentumEquation({ Y,
                            m_momentumVelocityCoeffs, 
                            m_gradientCoeffs[Y],
                            CFD::Tensor3D(),
                            CFD::Tensor3D() }), 
                            
         MomentumEquation({ Z,
                            m_momentumVelocityCoeffs, 
                            m_gradientCoeffs[Z],
                            CFD::Tensor3D(),
                            CFD::Tensor3D() }) 
        } )
{};



FVCoefficients::FVCoefficients( const iArray3 &dims, 
                                MomentumInterpolation mi) :
    m_momentumVelocityCoeffs( MomentumVelocityEnums() , dims + 2*nGhost ),
    m_gradientCoeffs( { CFD::EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::e, C::w}, dims( X ) + 2*nGhost),
                        CFD::EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::n, C::s}, dims( Y ) + 2*nGhost),
                        CFD::EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::t, C::b}, dims( Z ) + 2*nGhost) } ),

    Cont( { m_gradientCoeffs,
            EnumVector<TransportCoefficients, Tensor3D>( ContinuityPressureEnums( mi ), dims + 2*nGhost ),
            Tensor3D( dims(X) + 2*nGhost, dims(Y) + 2*nGhost, dims(Z) + 2*nGhost ).setZero(),
            Tensor3D( dims(X) + 2*nGhost, dims(Y) + 2*nGhost, dims(Z) + 2*nGhost ).setZero()
        } ),
    
    mwiSparseCoeffs( { std::array<Tensor1D, 4>{ Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero() } ,
                       std::array<Tensor1D, 4>{ Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero() } ,
                       std::array<Tensor1D, 4>{ Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero() } } ),
    mwiCompactCoeffs( { std::array<Tensor1D, 2>{ Tensor1D(dims(X)+1).setZero(), Tensor1D(dims(X)+1).setZero() } ,
                        std::array<Tensor1D, 2>{ Tensor1D(dims(Y)+1).setZero(), Tensor1D(dims(Y)+1).setZero() } ,
                        std::array<Tensor1D, 2>{ Tensor1D(dims(Z)+1).setZero(), Tensor1D(dims(Z)+1).setZero() } } ),
    momentumInterpolation( mi ),

    Mom( { MomentumEquation({ X,
                              m_momentumVelocityCoeffs, 
                              m_gradientCoeffs[X],
                              CFD::Tensor3D( dims(X) + 2*nGhost,  dims(Y) + 2*nGhost,  dims(Z) + 2*nGhost ).setZero(),
                              CFD::Tensor3D( dims(X) + 2*nGhost,  dims(Y) + 2*nGhost,  dims(Z) + 2*nGhost ).setZero() }),
              
              MomentumEquation({ Y,
                                 m_momentumVelocityCoeffs, 
                                 m_gradientCoeffs[Y],
                                 CFD::Tensor3D( dims(X) + 2*nGhost,  dims(Y) + 2*nGhost,  dims(Z) + 2*nGhost ).setZero(),
                                 CFD::Tensor3D( dims(X) + 2*nGhost,  dims(Y) + 2*nGhost,  dims(Z) + 2*nGhost ).setZero() }),

              MomentumEquation({ Z,
                                 m_momentumVelocityCoeffs, 
                                 m_gradientCoeffs[Z],
                                 CFD::Tensor3D( dims(X) + 2*nGhost,  dims(Y) + 2*nGhost,  dims(Z) + 2*nGhost ).setZero(),
                                 CFD::Tensor3D( dims(X) + 2*nGhost,  dims(Y) + 2*nGhost,  dims(Z) + 2*nGhost ).setZero() })
    } ),

    diff({ EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::e, C::w}, dims(X) ),
           EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::n, C::s}, dims(Y) ),
           EnumVector<TransportCoefficients, Tensor1D>( {C::p, C::t, C::b}, dims(Z) ) }),
       
    nCells( dims )
{};



// Copy constructor
FVCoefficients::FVCoefficients( const FVCoefficients &that ) :
    m_momentumVelocityCoeffs( that.m_momentumVelocityCoeffs ),
    m_gradientCoeffs( that.m_gradientCoeffs ),

    Cont( { m_gradientCoeffs,
            that.Cont.AP,
            that.Cont.B,
            that.Cont.F
        } ),
    
    mwiSparseCoeffs( that.mwiSparseCoeffs ),
    mwiCompactCoeffs( that.mwiCompactCoeffs),
    momentumInterpolation( that.momentumInterpolation ),


    Mom( { MomentumEquation({ X,
                              m_momentumVelocityCoeffs, 
                              m_gradientCoeffs[X],
                              that.Mom[X].B,
                              that.Mom[X].F }),
              
           MomentumEquation({ Y,
                              m_momentumVelocityCoeffs, 
                              m_gradientCoeffs[Y],
                              that.Mom[Y].B,
                              that.Mom[Y].F }),

           MomentumEquation({ Z,
                              m_momentumVelocityCoeffs, 
                              m_gradientCoeffs[Z],
                              that.Mom[Z].B,
                              that.Mom[Z].F })
        } ),
    diff( that.diff ),
    positiveFluxHiOrderAdvectionCoeffs( that.positiveFluxHiOrderAdvectionCoeffs ),
    negativeFluxHiOrderAdvectionCoeffs( that.negativeFluxHiOrderAdvectionCoeffs ),
    advectionScheme( that.advectionScheme ),
    timeScheme( that.timeScheme ),
    timeStep( that.timeStep ),
    advectionBlendingFactor( that.advectionBlendingFactor ),

    nu( that.nu ),
    rho( that.rho ),
    nCells( that.nCells )
{};



// Copy assignment
FVCoefficients& FVCoefficients::operator=( FVCoefficients that )
{
    // References to shared coefficients already refer to the correct object
    std::swap( this->m_momentumVelocityCoeffs , that.m_momentumVelocityCoeffs );
    std::swap( this->m_gradientCoeffs         , that.m_gradientCoeffs );

    // Continuity equation
    std::swap( this->Cont.AP, that.Cont.AP );
    std::swap( this->Cont.B , that.Cont.B );
    std::swap( this->Cont.F , that.Cont.F );

    std::swap( this->mwiSparseCoeffs      , that.mwiSparseCoeffs );
    std::swap( this->mwiCompactCoeffs     , that.mwiCompactCoeffs );
    std::swap( this->momentumInterpolation, that.momentumInterpolation );


    // Momentum equations
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        std::swap( this->Mom[axis].component, that.Mom[axis].component );
        std::swap( this->Mom[axis].B        , that.Mom[axis].B );
        std::swap( this->Mom[axis].F        , that.Mom[axis].F );
    } );
    std::swap( this->diff                              , that.diff );
    std::swap( this->positiveFluxHiOrderAdvectionCoeffs, that.positiveFluxHiOrderAdvectionCoeffs );
    std::swap( this->negativeFluxHiOrderAdvectionCoeffs, that.negativeFluxHiOrderAdvectionCoeffs );
    std::swap( this->advectionScheme                   , that.advectionScheme );
    std::swap( this->timeScheme                        , that.timeScheme );
    std::swap( this->timeStep                          , that.timeStep );
    std::swap( this->advectionBlendingFactor           , that.advectionBlendingFactor );

    std::swap( this->nu    , that.nu );
    std::swap( this->rho   , that.rho );
    std::swap( this->nCells, that.nCells );
    return *this;
}



// Move constructor, copy and swap
FVCoefficients::FVCoefficients( FVCoefficients &&that ) noexcept :
    m_momentumVelocityCoeffs( std::move( that.m_momentumVelocityCoeffs ) ),
    m_gradientCoeffs( std::move( that.m_gradientCoeffs ) ),

    Cont( { m_gradientCoeffs,
            std::move( that.Cont.AP ),
            std::move( that.Cont.B ),
            std::move( that.Cont.F )
        } ),
    
    mwiSparseCoeffs( std::move( that.mwiSparseCoeffs ) ),
    mwiCompactCoeffs( std::move( that.mwiCompactCoeffs) ),
    momentumInterpolation( std::move( that.momentumInterpolation ) ),

    Mom( { MomentumEquation({ X,
                              m_momentumVelocityCoeffs, 
                              m_gradientCoeffs[X],
                              std::move( that.Mom[X].B ),
                              std::move( that.Mom[X].F ) }),
              
           MomentumEquation({ Y,
                              m_momentumVelocityCoeffs, 
                              m_gradientCoeffs[Y],
                              std::move( that.Mom[Y].B ),
                              std::move( that.Mom[Y].F ) }),

           MomentumEquation({ Z,
                              m_momentumVelocityCoeffs, 
                              m_gradientCoeffs[Z],
                              std::move( that.Mom[Z].B ),
                              std::move( that.Mom[Z].F ) })
        } ),
    diff( std::move( that.diff ) ),
    positiveFluxHiOrderAdvectionCoeffs( std::move( that.positiveFluxHiOrderAdvectionCoeffs ) ),
    negativeFluxHiOrderAdvectionCoeffs( std::move( that.negativeFluxHiOrderAdvectionCoeffs ) ),
    advectionScheme( std::move( that.advectionScheme ) ),
    timeScheme( std::move( that.timeScheme ) ),
    timeStep( std::move( that.timeStep ) ),
    advectionBlendingFactor( std::move( that.advectionBlendingFactor ) ),

    nu( std::move( that.nu ) ),
    rho( std::move( that.rho ) ),
    nCells( std::move( that.nCells ) )
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