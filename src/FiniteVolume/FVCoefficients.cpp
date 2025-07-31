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
        m_momentumVelocityCoeffs[tc] = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost );
        SetTensorZeroParallel( m_momentumVelocityCoeffs[tc] );
    }

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        for ( auto tc : GradientCoefficientEnums(axis) ) {
            m_gradientCoeffs[axis][tc] = Tensor1D( nCells(axis) + 2*nGhost );
            SetTensorZeroParallel( m_gradientCoeffs[axis][tc] );
        }
    } );

    // Cont.AU set as reference to m_gradientCoeffs
    for ( auto tc : ContinuityPressureEnums( momentumInterpolation ) ) {
        Cont.AP[tc] = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost );
        SetTensorZeroParallel( Cont.AP[tc] );
    }

    Cont.B = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost );
    SetTensorZeroParallel( Cont.B  );

    Cont.F = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost );
    SetTensorZeroParallel( Cont.F );

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        mwiSparseCoeffs[axis]  = { Tensor1D(nCells(axis)+1).setZero(), Tensor1D(nCells(axis)+1).setZero(), Tensor1D(nCells(axis)+1).setZero(), Tensor1D(nCells(axis)+1).setZero() };
        mwiCompactCoeffs[axis] = { Tensor1D(nCells(axis)+1).setZero(), Tensor1D(nCells(axis)+1).setZero() };
    } );


    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        Mom[axis].component = axis;
        // Mom[axis].AU set as reference to m_momentumVelocityCoeffs
        // Mom[axis].AP set as reference to m_gradientCoeffs[axis]
        Mom[axis].B = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost );
        SetTensorZeroParallel( Mom[axis].B );

        Mom[axis].F = Tensor3D( nCells(X) + 2*nGhost, nCells(Y) + 2*nGhost, nCells(Z) + 2*nGhost );
        SetTensorZeroParallel( Mom[axis].F );

    } );

    nuTurb = Tensor3D( nCells(X) + 2*nGhost, nCells(X) + 2*nGhost, nCells(X) + 2*nGhost );
    SetTensorZeroParallel( nuTurb );

};



}   // end namespace CAMIRA