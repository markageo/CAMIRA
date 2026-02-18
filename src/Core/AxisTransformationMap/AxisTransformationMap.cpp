#include "AxisTransformationMap.h"
#include "../FVLookups.h"

namespace CAMIRA
{

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
    A codeAxis = LUT::BoundaryPatchAxis[ codePatch ];
    A userAxis = LUT::BoundaryPatchAxis[ userPatch ];

    // Set the given patch
    m_codeBoundaryPatches[ codePatch ] = userPatch;
    m_userBoundaryPatches[ userPatch ] = codePatch;

    // Set the opposite patch
    BP codePatchOpposite;
    
    if ( codePatch == LUT::PositivePatch[ codeAxis ] ) {
        codePatchOpposite =LUT:: NegativePatch[ codeAxis ];
    } else {
        codePatchOpposite = LUT::PositivePatch[ codeAxis ];
    }

    BP userPatchOpposite;
    if ( userPatch == LUT::PositivePatch[ userAxis ] ) {
        userPatchOpposite = LUT::NegativePatch[ userAxis ];
    } else {
        userPatchOpposite = LUT::PositivePatch[ userAxis ];
    }

    m_codeBoundaryPatches[ codePatchOpposite ] = userPatchOpposite;
    m_userBoundaryPatches[ userPatchOpposite ] = codePatchOpposite;


    // Update the axis
    m_codeAxis[ codeAxis ] = userAxis;
    m_userAxis[ userAxis ] = codeAxis;
}


// Code patch from user patch
BoundaryPatches::ENUMDATA AxisTransformationMap::CodePatch(const BP userPatch) const 
{   return m_userBoundaryPatches[ userPatch ]; }

// Code axis from user axis
Axis::ENUMDATA AxisTransformationMap::CodeAxis(const A userAxis) const
{ 
    BP userPatch = LUT::PositivePatch[ userAxis ];
    BP codePatch = CodePatch( userPatch );
    return LUT::BoundaryPatchAxis[ codePatch ];
}

// If code axis is mapped to the negative direction of a user axis
bool AxisTransformationMap::CodeAxisReversed( const A codeAxis ) const
{
    BP positiveCodePatch = LUT::PositivePatch[ codeAxis ];
    BP mappedUserPatch   = UserPatch( positiveCodePatch );
    A  mappedUserAxis    = LUT::BoundaryPatchAxis[ mappedUserPatch ];
    if ( mappedUserPatch == LUT::NegativePatch[ mappedUserAxis ] ) {
        return true;
    }
    return false;
}

floatType AxisTransformationMap::CodeAxisReverseSign( const A codeAxis ) const
{
    if ( CodeAxisReversed( codeAxis ) ) {
        return -1;
    } else {
        return +1;
    }
}



// User patch from code patch
BoundaryPatches::ENUMDATA AxisTransformationMap::UserPatch(const BP codePatch) const
{   return m_codeBoundaryPatches[ codePatch ]; }

// User axis from code axis
Axis::ENUMDATA AxisTransformationMap::UserAxis(const A codeAxis) const
{
    BP codePatch = LUT::PositivePatch[ codeAxis ];
    BP userPatch = UserPatch( codePatch );
    return LUT::BoundaryPatchAxis[ userPatch ];
}

// If user axis is mapped to the negative direction of a code axis
bool AxisTransformationMap::UserAxisReversed( const A userAxis ) const
{
    BP positiveUserPatch = LUT::PositivePatch[ userAxis ];
    BP mappedCodePatch   = CodePatch( positiveUserPatch );
    A  mappedCodeAxis    = LUT::BoundaryPatchAxis[ mappedCodePatch ];
    if ( mappedCodePatch == LUT::NegativePatch[ mappedCodeAxis ] ) {
        return true;
    }
    return false;
}


floatType AxisTransformationMap::UserAxisReverseSign( const A userAxis ) const
{
    if ( UserAxisReversed( userAxis ) ) {
        return -1;
    } else {
        return +1;
    }
}


AxisTransformationMap AxisTransformationMap::Inverse() const
{
    AxisTransformationMap inverseAxisTrasformation;

    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        inverseAxisTrasformation.Set( bp, this->CodePatch( bp ) );  // This is opposite to the definition of the Set function
    } );

    return inverseAxisTrasformation;
}



AxisTransformationMap AxisTransformationMap::Identity() const
{
    AxisTransformationMap identityAxisTransformation;

    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA bp) {
        identityAxisTransformation.Set( bp, bp );
    } );

    return identityAxisTransformation;
}


}   // end namespace CAMIRA