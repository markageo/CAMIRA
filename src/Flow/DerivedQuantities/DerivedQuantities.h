#ifndef CAMIRA_DERIVED_QUANTITIES
#define CAMIRA_DERIVED_QUANTITIES

#include "../../Core/Types.h"
#include "../../Core/Mesh/Mesh.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../FiniteVolume/BoundaryConditionData.h"

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



// Calculates some y+ statistics
class YPlusCalculator
{
    public: 
        YPlusCalculator( const InputData &,
                         const AxisTransformationMap &,
                         const BoundaryConditionData &, 
                         const IBData &, 
                         const Mesh &, 
                         const FieldData<Tensor3D> & );

        void Update();
        floatType minYPlus, maxYPlus, averageYPlus;

    private:
        const FieldData<Tensor3D> &m_fields;
        const floatType m_rho, m_nu;

        struct WallCellData {
            floatType yPlus;
            TensorIndex3D cellIndex;
            floatType wallDistance;
            fVector3 normalVector;
            fVector3 wallTangentialVelocity;
        };

        std::vector< WallCellData > m_wallCells;
};



}   // end namespace CAMIRA

#endif // CAMIRA_DERIVED_QUANTITIES