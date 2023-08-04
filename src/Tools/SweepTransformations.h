#ifndef SWEEP_TRANSFORMATIONS
#define SWEEP_TRANSFORMATIONS

#include "../Types.h"
#include "../IO/InputProcessing.h"
#include "../FiniteVolume/FiniteVolume.h"

// The problem is remapped so that the plane sweeping direction is always in the z direction and
// the line sweeping direction is always in the y direction (in the code). This is more memory
// efficient and is simpler to implement.

namespace CFD
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
       

        // User patch from code patch
        BP UserPatch(const BP) const;

        // User axis from code axis
        A UserAxis(const A) const;

        // If user axis is mapped to the negative direction of a code axis
        bool UserAxisReversed( const A ) const;


    private:

        // Lookups for patches
        EnumVector<BoundaryPatches, BP> m_codeBoundaryPatches;   // Code patch -> user patch
        EnumVector<BoundaryPatches, BP> m_userBoundaryPatches;   // User patch -> code patch

        // Lookups for axis
        EnumVector<Axis, A> m_codeAxis;    // Code axis -> user axis
        EnumVector<Axis, A> m_userAxis;    // User axis -> code axis
};


// Return a transformation map from a plane and line sweeping direction
AxisTransformationMap CreateAxisTransformation( BoundaryPatches::ENUMDATA, BoundaryPatches::ENUMDATA ); 

// Transform user input data so that sweeping is consistent with solver
AxisTransformationMap TransformUserInputData( InputData & );

// Transform back to the coordinates consistentent with the input file 
void TransformToUserCoordinates(Mesh &, FieldData<array3D> &, FieldData<array3D> &, const AxisTransformationMap &);


} // end namespace CFD
 
#endif  // SWEEP_TRANSFORMATIONS