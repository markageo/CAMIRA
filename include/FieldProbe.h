#ifndef FIELD_PROBE
#define FIELD_PROBE

#include "Types.h"
#include "FiniteVolume.h"

namespace CFD
{

// https://en.wikipedia.org/wiki/Trilinear_interpolation
class FieldProbe
{
    public:

        FieldProbe( const Mesh &mesh,
                    const fVector3 &probePoint )
        {
            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                for ( intType i = 1; i != mesh.cellCenters[axis].size(); i++ ) {

                    floatType x0 = mesh.cellCenters[axis](i-1);
                    floatType x1 = mesh.cellCenters[axis](i);

                    bool pointHigher = probePoint[axis] >= x0;
                    bool pointLower  = probePoint[axis] <= x1;

                    if ( pointHigher && pointLower ) {
                        m_latticeIndex[axis](0) = i-1;
                        m_latticeIndex[axis](1) = i;
                        m_latticeCoord[axis] = ( probePoint[axis] - x0 ) / ( x1 - x0 );
                        break;
                    }
                }

            } );

        }

        // For arrays with ghost cells
        floatType GetFieldValue( const array3D &field )
        {
            using enum Axis::ENUMDATA;
            const floatType c000 = field( G( m_latticeIndex[X](0), m_latticeIndex[Y](0), m_latticeIndex[Z](0) ) ),
                            c100 = field( G( m_latticeIndex[X](1), m_latticeIndex[Y](0), m_latticeIndex[Z](0) ) ),
                            c010 = field( G( m_latticeIndex[X](0), m_latticeIndex[Y](1), m_latticeIndex[Z](0) ) ),
                            c001 = field( G( m_latticeIndex[X](0), m_latticeIndex[Y](0), m_latticeIndex[Z](1) ) ),
                            c101 = field( G( m_latticeIndex[X](1), m_latticeIndex[Y](0), m_latticeIndex[Z](1) ) ),
                            c011 = field( G( m_latticeIndex[X](0), m_latticeIndex[Y](1), m_latticeIndex[Z](1) ) ),
                            c110 = field( G( m_latticeIndex[X](1), m_latticeIndex[Y](1), m_latticeIndex[Z](0) ) ),
                            c111 = field( G( m_latticeIndex[X](1), m_latticeIndex[Y](1), m_latticeIndex[Z](1) ) );

            // Linear interpolation in x direction
            const floatType c00 = c000 * ( 1 - m_latticeCoord[X] ) + c100 *  m_latticeCoord[X],
                            c01 = c001 * ( 1 - m_latticeCoord[X] ) + c101 *  m_latticeCoord[X],
                            c10 = c010 * ( 1 - m_latticeCoord[X] ) + c110 *  m_latticeCoord[X],
                            c11 = c011 * ( 1 - m_latticeCoord[X] ) + c111 *  m_latticeCoord[X];

            // Linear interpolation in y direction
            const floatType c0 = c00 * ( 1 - m_latticeCoord[Y] ) + c10 * m_latticeCoord[Y],
                            c1 = c01 * ( 1 - m_latticeCoord[Y] ) + c11 * m_latticeCoord[Y];

            // Linear interpolation in z direction
            const floatType c = c0 * ( 1 - m_latticeCoord[Z] ) + c1 * m_latticeCoord[Z];

            return c;
        }

    private:
        EnumVector<Axis, iVector2> m_latticeIndex;  // Index of lattice points to itnerpolate from
        EnumVector<Axis, floatType> m_latticeCoord;   // Normalised coordinates in the lattice
    
};

}   // end namespace CFD

#endif // FIELD_PROBE