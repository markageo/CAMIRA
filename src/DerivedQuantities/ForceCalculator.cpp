#include "DerivedQuantities.h"
#include "../Core/FVTools.h"

namespace CAMIRA
{

ForceCalculator::ForceCalculator( const IBData &ibData, 
                                  const Mesh &mesh, 
                                  const FieldData<Tensor3D> &fields,
                                  const floatType rho,
                                  const floatType nu ) : 
                                    m_ibData( ibData ), 
                                    m_mesh(mesh), 
                                    m_fields(fields), 
                                    m_rho(rho), 
                                    m_nu(nu) 
                                {};


fVector3 ForceCalculator::GetForce() const
{
    using enum Axis::ENUMDATA;

    fVector3 pressureForce{ 0.0f, 0.0f, 0.0f },
             viscousForce{ 0.0f, 0.0f, 0.0f },
             advectiveForce{ 0.0f, 0.0f, 0.0f };

    for ( const auto &ibCellComponent : m_ibData.ibCells ) {
        for ( const auto &ibCell : ibCellComponent ) { 
            for ( const auto &sourceTermData : ibCell.sourceTermsData ) {

                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                    pressureForce(axis)  += GetPressureIntegralForce(axis, sourceTermData);
                    viscousForce(axis)   += GetViscousIntegralForce(axis, sourceTermData, ibCell.cellIndex);
                    advectiveForce(axis) += GetAdvectiveIntegralForce(axis, sourceTermData);
                } );

            }
        }
    }

    fVector3 totalForce = pressureForce + viscousForce + advectiveForce;

    return totalForce;

}


floatType ForceCalculator::GetViscousIntegralForce( const Axis::ENUMDATA forceComponent,
                                                    const IBCell::SourceTermData &sourceTermData,
                                                    const TensorIndex3D &cellIndex ) const
{
    using FVT::G;

    // Gradient in face normal direction
    floatType u2 = sourceTermData.faceValues.U[forceComponent],
              u1 = m_fields.U[forceComponent]( G( cellIndex ) ),
              dx = m_mesh.cellLengths[sourceTermData.direction]( cellIndex[sourceTermData.direction] ) / 2.0f;

    floatType grad =  static_cast<floatType>( sourceTermData.directionIndex )  // Direction index to account for whether face is hi or lo side of cell
                   * ( u2 - u1 ) / dx;        

    return - m_rho * m_nu * grad * sourceTermData.faceAreaComponent;
}


floatType ForceCalculator::GetPressureIntegralForce( const Axis::ENUMDATA forceComponent,
                                                     const IBCell::SourceTermData &sourceTermData ) const
{
    if ( forceComponent == sourceTermData.direction ) {
        return sourceTermData.faceValues.P * sourceTermData.faceAreaComponent;
    } else {
        return 0.0f;
    } 
}


floatType ForceCalculator::GetAdvectiveIntegralForce( const Axis::ENUMDATA forceComponent,
                                                      const IBCell::SourceTermData &sourceTermData ) const
{
    return   m_rho 
           * sourceTermData.faceValues.U[ forceComponent ]
           * sourceTermData.faceValues.U[ sourceTermData.direction ]
           * sourceTermData.faceAreaComponent;
}




}   // end namespace CAMIRA

