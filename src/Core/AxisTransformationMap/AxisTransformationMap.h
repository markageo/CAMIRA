#ifndef CAMIRA_AXIS_TRANSFORMATION_MAP
#define CAMIRA_AXIS_TRANSFORMATION_MAP

#include "../Types.h"

namespace CAMIRA
{

// Structure for storing axis transformation, is just a one-to-one map
class AxisTransformationMap
{
    using BP = BoundaryPatches::ENUMDATA;
    using A = Axis::ENUMDATA;

    public:
        AxisTransformationMap();

        // Setting values
        void Set(const BP, const BP);


        // Code patch from user patch
        BP CodePatch(const BP) const;

        // Code axis from user axis
        A CodeAxis(const A) const;

        // If code axis is mapped to the negative direction of a user axis
        bool CodeAxisReversed( const A ) const;

        // Same as CodeAxisReversed but returns +1 if no reversed and -1 if reversed
        floatType CodeAxisReverseSign( const A ) const;
       

        // User patch from code patch
        BP UserPatch(const BP) const;

        // User axis from code axis
        A UserAxis(const A) const;

        // If user axis is mapped to the negative direction of a code axis
        bool UserAxisReversed( const A ) const;

        // Same as UserAxisReversed but returns +1 if no reversed and -1 if reversed
        floatType UserAxisReverseSign( const A ) const;


        // Returns an axis transformation map that is the inverse of the current one (i.e. code and user axis swapped around)
        AxisTransformationMap Inverse() const;

        // Returns an axis transformation map that is the identity (i.e. no transformation is applied)
        AxisTransformationMap Identity() const;


    private:

        // Lookups for patches
        EnumVector<BoundaryPatches, BP> m_codeBoundaryPatches;   // Code patch -> user patch
        EnumVector<BoundaryPatches, BP> m_userBoundaryPatches;   // User patch -> code patch

        // Lookups for axis
        EnumVector<Axis, A> m_codeAxis;    // Code axis -> user axis
        EnumVector<Axis, A> m_userAxis;    // User axis -> code axis
};

} // end namespace CAMIRA
 
#endif  // CAMIRA_AXIS_TRANSFORMATION_MAP