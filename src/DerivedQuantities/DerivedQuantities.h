#ifndef CFD_DERIVED_QUANTITIES
#define CFD_DERIVED_QUANTITIES

#include "../Types.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../FiniteVolume/Mesh.h"

namespace CFD
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


// Calculate forces exerted by fluid on immersed boundary geometry
class ForceCalculator
{
    public: 
        ForceCalculator( const IBData &, const Mesh &, const FieldData<Tensor3D> &, const floatType rho, const floatType nu );
        fVector3 GetForce() const;

    private:
        floatType GetPressureIntegralForce( const Axis::ENUMDATA, const IBCell::SourceTermData & ) const;
        floatType GetViscousIntegralForce( const Axis::ENUMDATA, const IBCell::SourceTermData &, const TensorIndex3D & ) const;
        floatType GetAdvectiveIntegralForce( const Axis::ENUMDATA, const IBCell::SourceTermData & ) const;

        const IBData &m_ibData;
        const Mesh &m_mesh;
        const FieldData<Tensor3D> &m_fields;
        const floatType m_rho, m_nu;
};



}   // end namespace CFD

#endif // CFD_DERIVED_QUANTITIES