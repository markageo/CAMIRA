#ifndef CAMIRA_FIELD_PROBE
#define CAMIRA_FIELD_PROBE

#include "../Core/Types.h"
#include "../FiniteVolume/Mesh.h"

namespace CAMIRA
{

// Calculate values and gradients between points using linear interpolation
class FieldProbe
{
    public:
        FieldProbe( const Mesh &, const fArray3 & );
        floatType GetFieldValue( const Tensor3D & ) const;
        fVector3 GetFieldGradient( const Tensor3D & ) const;
        const fArray3& Coordinates() const;
        floatType Coordinate( const intType ) const;

    private:
        fArray3 m_probePoint;
        fArray3 m_latticeDims;
        EnumVector<Axis, iArray2> m_latticeIndex;     // Index of lattice points to interpolate from
        EnumVector<Axis, floatType> m_latticeCoord;   // Normalised coordinates in the lattice
};


// Get probe values for all FieldData
inline FieldData<floatType> ProbeAllFieldValues( const FieldData<Tensor3D> &fields,
                                                 const FieldProbe &fieldProbe )
{
    FieldData<floatType> probeValues;

    ForAllFieldData( [&] (intType f) { 
        probeValues[f] =  fieldProbe.GetFieldValue( fields[f] );
    } );

    return probeValues;
}



}   // end namespace CAMIRA

#endif // CAMIRA_FIELD_PROBE