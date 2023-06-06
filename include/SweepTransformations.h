#ifndef SWEEP_TRANSFORMATIONS
#define SWEEP_TRANSFORMATIONS

#include "Types.h"
#include "InputProcessing.h"
#include "FiniteVolume.h"

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
        void Set(const BP codePatch, const BP userPatch)
        {
            m_codeMap[codePatch] = userPatch;
            m_userMap[userPatch] = codePatch; 
        }


        // Code patch from user patch
        const BP &CodePatch(const BP userPatch) const 
        { return m_userMap.at( userPatch ); }

        // Code axis from user axis
        const A &CodeAxis(const A userAxis) const
        { 
            BP userPatch = PositivePatch[ userAxis ];
            BP codePatch = CodePatch( userPatch );
            return BoundaryPatchAxis[ codePatch ];
        }


        // User patch from code patch
        const BP &UserPatch(const BP codePatch) const
        { return m_codeMap.at( codePatch ); }

        // User axis from code axis
        const A &UserAxis(const A codeAxis) const
        {
            BP codePatch = PositivePatch[ codeAxis ];
            BP userPatch = UserPatch( codePatch );
            return BoundaryPatchAxis[ userPatch ];
        }

    private:
        std::map< BP, BP> m_codeMap;    // Code patch -> user patch
        std::map< BP, BP> m_userMap;    // User patch -> code patch
};


// Transform user input data so that sweeping is consistent with solver
AxisTransformationMap TransformUserInputData( InputData & );

// Transform back to the coordinates consistentent with the input file 
void TransformToUserCoordinates(Mesh &, ArrayAllocator<Fields, array3D> &, const AxisTransformationMap &);


} // end namespace CFD

#endif  // SWEEP_TRANSFORMATIONS