#include "DerivedQuantities.h"
#include "../../Core/FVTools.h"
#include "../../Core/Types.h"
#include "../FiniteVolume/FiniteVolume.h"


namespace CAMIRA
{


FieldProbe::FieldProbe( const Mesh &mesh,
            const fArray3 &probePoint ) :
    m_probePoint( probePoint )
{
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        
        // These values are valid if the mesh is only 1 cell thick
        m_latticeIndex[axis](0) = 0;
        m_latticeIndex[axis](1) = 0;
        m_latticeCoord[axis]    = 0.0f;

        for ( intType i = 1; i != mesh.cellCenters[axis].size(); i++ ) {

            floatType x0 = mesh.cellCenters[axis](i-1);
            floatType x1 = mesh.cellCenters[axis](i);

            bool pointHigher = probePoint[axis] >= x0;
            bool pointLower  = probePoint[axis] <= x1;

            if ( pointHigher && pointLower ) {
                m_latticeIndex[axis](0)  = i-1;
                m_latticeIndex[axis](1)  = i;
                m_latticeDims(axis) = x1 - x0;
                m_latticeCoord[axis]     = ( probePoint[axis] - x0 ) / ( x1 - x0 );
                break;
            }
        }

    } );

}


// For arrays with ghost cells
// https://en.wikipedia.org/wiki/Trilinear_interpolation
floatType FieldProbe::GetFieldValue( const Tensor3D &field ) const
{
    using namespace FVT;
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
    const floatType c00 = c000 * ( 1 - m_latticeCoord[X] )  +  c100 * m_latticeCoord[X],
                    c01 = c001 * ( 1 - m_latticeCoord[X] )  +  c101 * m_latticeCoord[X],
                    c10 = c010 * ( 1 - m_latticeCoord[X] )  +  c110 * m_latticeCoord[X],
                    c11 = c011 * ( 1 - m_latticeCoord[X] )  +  c111 * m_latticeCoord[X];

    // Linear interpolation in y direction
    const floatType c0 = c00 * ( 1 - m_latticeCoord[Y] )  +  c10 * m_latticeCoord[Y],
                    c1 = c01 * ( 1 - m_latticeCoord[Y] )  +  c11 * m_latticeCoord[Y];

    // Linear interpolation in z direction
    const floatType c = c0 * ( 1 - m_latticeCoord[Z] )  +  c1 * m_latticeCoord[Z];

    return c;
}


// For arrays with ghost cells
fVector3 FieldProbe::GetFieldGradient( const Tensor3D &field ) const
{
    using namespace FVT;
    using enum Axis::ENUMDATA;

    fVector3 gradientVector;

    const floatType c000 = field( G( m_latticeIndex[X](0), m_latticeIndex[Y](0), m_latticeIndex[Z](0) ) ),
                    c100 = field( G( m_latticeIndex[X](1), m_latticeIndex[Y](0), m_latticeIndex[Z](0) ) ),
                    c010 = field( G( m_latticeIndex[X](0), m_latticeIndex[Y](1), m_latticeIndex[Z](0) ) ),
                    c001 = field( G( m_latticeIndex[X](0), m_latticeIndex[Y](0), m_latticeIndex[Z](1) ) ),
                    c101 = field( G( m_latticeIndex[X](1), m_latticeIndex[Y](0), m_latticeIndex[Z](1) ) ),
                    c011 = field( G( m_latticeIndex[X](0), m_latticeIndex[Y](1), m_latticeIndex[Z](1) ) ),
                    c110 = field( G( m_latticeIndex[X](1), m_latticeIndex[Y](1), m_latticeIndex[Z](0) ) ),
                    c111 = field( G( m_latticeIndex[X](1), m_latticeIndex[Y](1), m_latticeIndex[Z](1) ) );
        
    // x gradient -----------------------------------------------------------------------------------------
    const floatType gx00 = ( c100 - c000 ) / m_latticeDims[X],
                    gx10 = ( c110 - c010 ) / m_latticeDims[X],
                    gx01 = ( c101 - c001 ) / m_latticeDims[X],
                    gx11 = ( c111 - c011 ) / m_latticeDims[X];

    // Interpolate y direction
    const floatType gx0 = gx00 * ( 1 - m_latticeCoord[Y] )  +  gx10 * m_latticeCoord[Y],
                    gx1 = gx01 * ( 1 - m_latticeCoord[Y] )  +  gx11 * m_latticeCoord[Y];

    // Interpolate z direction
    const floatType gx = gx0 * ( 1 - m_latticeCoord[Z] )  +  gx1 * m_latticeCoord[Z];
    gradientVector(X) = gx;


    // y gradient -----------------------------------------------------------------------------------------
    const floatType gy00 = ( c010 - c000 ) / m_latticeDims[Y],
                    gy10 = ( c110 - c100 ) / m_latticeDims[Y],
                    gy01 = ( c011 - c001 ) / m_latticeDims[Y],
                    gy11 = ( c111 - c101 ) / m_latticeDims[Y];

    // Interpolate x direction
    const floatType gy0 = gy00 * ( 1 - m_latticeCoord[X] )  +  gy10 * m_latticeCoord[X],
                    gy1 = gy01 * ( 1 - m_latticeCoord[X] )  +  gy11 * m_latticeCoord[X];

    // Interpolate z direction
    const floatType gy = gy0 * ( 1 - m_latticeCoord[Z] )  +  gy1 * m_latticeCoord[Z];
    gradientVector(Y) = gy;


    // z gradient -----------------------------------------------------------------------------------------
    const floatType gz00 = ( c001 - c000 ) / m_latticeDims[Z],
                    gz10 = ( c101 - c100 ) / m_latticeDims[Z],
                    gz01 = ( c011 - c010 ) / m_latticeDims[Z],
                    gz11 = ( c111 - c110 ) / m_latticeDims[Z];

    // Interpolate x direction
    const floatType gz0 = gz00 * ( 1 - m_latticeCoord[X] )  +  gz10 * m_latticeCoord[X],
                    gz1 = gz01 * ( 1 - m_latticeCoord[X] )  +  gz11 * m_latticeCoord[X];

    // Interpolate y direction
    const floatType gz = gz0 * ( 1 - m_latticeCoord[Y] )  +  gz1 * m_latticeCoord[Y];
    gradientVector(Z) = gz;


    return gradientVector;
}


const fArray3& FieldProbe::Coordinates() const
{ return m_probePoint; }


floatType FieldProbe::Coordinate( const intType axis ) const
{ return m_probePoint( axis ); }


}   // end namespace CAMIRA

