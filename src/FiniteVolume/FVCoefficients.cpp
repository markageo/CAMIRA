#include "FVCoefficients.h"
#include "../Core/FVTools.h"
#include "../Core/FVLookups.h"


namespace CAMIRA
{


namespace
{

std::vector< TransportCoefficients::ENUMDATA > MomentumVelocityEnums()
{
    using C = TransportCoefficients::ENUMDATA;
    return {C::p, C::n, C::e, C::s, C::w, C::t, C::b};
}                                                                      


std::vector< TransportCoefficients::ENUMDATA > GradientCoefficientEnums( Axis::ENUMDATA axis )
{
    using C = TransportCoefficients::ENUMDATA;
    using enum Axis::ENUMDATA;

    switch ( axis ) {
        case X:
            return {C::p, C::e, C::w};

        case Y:
            return {C::p, C::n, C::s};

        case Z:
            return {C::p, C::t, C::b};
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

// Default constructor must be defined explicity since class holds references
FVCoefficients::FVCoefficients() :
    Cont( { m_gradientCoeffs,
            CAMIRA::EnumVector<TransportCoefficients, Tensor3D>(),
            CAMIRA::Tensor3D(),
            CAMIRA::Tensor3D() } ),

    Mom( { MomentumEquation({ X,
                            m_momentumVelocityCoeffs, 
                            m_gradientCoeffs[X],
                            CAMIRA::Tensor3D(),
                            CAMIRA::Tensor3D() }),
                            
         MomentumEquation({ Y,
                            m_momentumVelocityCoeffs, 
                            m_gradientCoeffs[Y],
                            CAMIRA::Tensor3D(),
                            CAMIRA::Tensor3D() }), 
                            
         MomentumEquation({ Z,
                            m_momentumVelocityCoeffs, 
                            m_gradientCoeffs[Z],
                            CAMIRA::Tensor3D(),
                            CAMIRA::Tensor3D() }) 
        } )
{};



FVCoefficients::FVCoefficients( const iArray3 &dims, 
                                MomentumInterpolation mi) :
    FVCoefficients()   // Call the default constructor to set references
{
    using enum Axis::ENUMDATA;
    nCells = dims;
    momentumInterpolation = mi;

    for ( auto tc : MomentumVelocityEnums() ) {
        m_momentumVelocityCoeffs[tc] = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost ).setZero();
    }

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        for ( auto tc : GradientCoefficientEnums(axis) ) {
            m_gradientCoeffs[axis][tc] = Tensor1D( nCells(axis) + 2*nGhost ).setZero();;
        }
    } );

    // Cont.AU set as reference to m_gradientCoeffs
    for ( auto tc : ContinuityPressureEnums( momentumInterpolation ) ) {
        Cont.AP[tc] = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost ).setZero();
    }
    Cont.B = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost ).setZero();
    Cont.F = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost ).setZero();

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        mwiSparseCoeffs[axis]  = { Tensor1D(nCells(axis)+1).setZero(), Tensor1D(nCells(axis)+1).setZero(), Tensor1D(nCells(axis)+1).setZero(), Tensor1D(nCells(axis)+1).setZero() };
        mwiCompactCoeffs[axis] = { Tensor1D(nCells(axis)+1).setZero(), Tensor1D(nCells(axis)+1).setZero() };
    } );


    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        Mom[axis].component = axis;
        // Mom[axis].AU set as reference to m_momentumVelocityCoeffs
        // Mom[axis].AP set as reference to m_gradientCoeffs[axis]
        Mom[axis].B = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost ).setZero();
        Mom[axis].F = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost ).setZero();
    } );

    // Do not allocate nuTurb here as may not be used 

};



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
    positiveFluxHiOrderAdvectionCoeffs( that.positiveFluxHiOrderAdvectionCoeffs ),
    negativeFluxHiOrderAdvectionCoeffs( that.negativeFluxHiOrderAdvectionCoeffs ),
    advectionScheme( that.advectionScheme ),
    timeScheme( that.timeScheme ),
    timeStep( that.timeStep ),
    advectionBlendingFactor( that.advectionBlendingFactor ),

    nuTurb( that.nuTurb ),

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
    std::swap( this->positiveFluxHiOrderAdvectionCoeffs, that.positiveFluxHiOrderAdvectionCoeffs );
    std::swap( this->negativeFluxHiOrderAdvectionCoeffs, that.negativeFluxHiOrderAdvectionCoeffs );
    std::swap( this->advectionScheme                   , that.advectionScheme );
    std::swap( this->timeScheme                        , that.timeScheme );
    std::swap( this->timeStep                          , that.timeStep );
    std::swap( this->advectionBlendingFactor           , that.advectionBlendingFactor );

    std::swap( this->nuTurb, that.nuTurb );

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
    positiveFluxHiOrderAdvectionCoeffs( std::move( that.positiveFluxHiOrderAdvectionCoeffs ) ),
    negativeFluxHiOrderAdvectionCoeffs( std::move( that.negativeFluxHiOrderAdvectionCoeffs ) ),
    advectionScheme( std::move( that.advectionScheme ) ),
    timeScheme( std::move( that.timeScheme ) ),
    timeStep( std::move( that.timeStep ) ),
    advectionBlendingFactor( std::move( that.advectionBlendingFactor ) ),

    nuTurb( std::move( that.nuTurb ) ),

    nu( std::move( that.nu ) ),
    rho( std::move( that.rho ) ),
    nCells( std::move( that.nCells ) )
{};


}   // end namespace CAMIRA