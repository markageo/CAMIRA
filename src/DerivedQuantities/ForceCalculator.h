#ifndef CAMIRA_FORCE_CALCULATOR
#define CAMIRA_FORCE_CALCULATOR

#include "../Core/Types.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../FiniteVolume/Mesh.h"

namespace CAMIRA
{


// Calculate forces exerted by fluid on immersed boundary geometry
class ForceCalculator
{
    public: 
        ForceCalculator( const IBData &, const Mesh &, const FieldData<Tensor3D> &, const floatType, const floatType );
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



}   // end namespace CAMIRA

#endif // CAMIRA_FORCE_CALCULATOR