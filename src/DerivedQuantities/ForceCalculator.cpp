#include "DerivedQuantities.h"
#include "../Core/FVTools.h"

namespace CFD
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

    fVector3 force{ 0.0f, 0.0f, 0.0f };

    for ( const auto &ibCell : m_ibData.ibCells ) { 
        for ( const auto &sourceTermData : ibCell.sourceTermsData ) {

            EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
                force( axis ) += GetPressureIntegralForce(axis, sourceTermData)
                               + GetViscousIntegralForce(axis, sourceTermData, ibCell.cellIndex)
                               + GetAdvectiveIntegralForce(axis, sourceTermData);
            } );

        }
    }

    return force;

}


namespace 
{
    // Calculates first or derivative accounting for ghost cell
    floatType CalculateVelocityGradient( const Axis::ENUMDATA gradDirection,
                                         const TensorIndex3D &cellIndex,
                                         const Axis::ENUMDATA &velocityComponent,
                                         const EnumVector<Axis, Tensor3D> &velocityFields,
                                         const IBData &ibData,
                                         const Mesh &mesh )
    {
        using FVT::G;

        const Tensor3D &U = velocityFields[velocityComponent];

        floatType up = U( G( cellIndex ) ),
                  uHi, uLo;
        bool uHiSet = false,
             uLoSet = false;

        // Check if current cell is an ib cell (has source terms added due to IB)
        for ( const auto &ibCell : ibData.ibCells ) {
            
            if ( ibCell.cellIndex != cellIndex )
                continue;
            
            for ( const auto &sourceTermData : ibCell.sourceTermsData ) {

                if ( sourceTermData.direction != gradDirection ) 
                    continue;

                if ( sourceTermData.directionIndex == +1 ) {
                    uHi    = sourceTermData.ghostCellValues.U[velocityComponent];
                    uHiSet = true;
                } else {
                    uLo    = sourceTermData.ghostCellValues.U[velocityComponent];
                    uLoSet = true;
                }
                
            }

            break;
        }

        if ( !uHiSet ) {
            TensorIndex3D cellIndexHi = cellIndex;
            cellIndexHi[gradDirection] += 1;
            uHi = U( G( cellIndexHi ) );
        }
        
        if ( !uLoSet ) {
            TensorIndex3D cellIndexLo = cellIndex;
            cellIndexLo[gradDirection] -= 1;
            uLo = U( G( cellIndexLo ) );
        }

        intType idx = cellIndex[gradDirection];
        floatType lambdaHi = mesh.interpFactors[gradDirection](idx+1),
                  lambdaLo = mesh.interpFactors[gradDirection](idx  ),
                  dx       = mesh.cellLengthsInv[gradDirection](idx);

        floatType gradient = (   lambdaHi                  * uHi
                               + (1 - lambdaHi - lambdaLo) * up
                               - (1 - lambdaLo)            * uLo
                             ) / dx;

        return gradient;
    }


    floatType ExtrapolateVelocityGradient( const floatType grad_p, 
                                           const floatType grad_a,
                                           const IBCell::SourceTermData &sourceTermData )
    {
        // First extrapolate to immersed boundary surface, we do it like this because we already have the coefficients stored
        floatType grad_ib = sourceTermData.ibExtrapFactor_p * grad_p    
                          + sourceTermData.ibExtrapFactor_a * grad_a;

        // Then to approximated boundary face
        floatType grad_f  = sourceTermData.faceExtrapCoeff_p  * grad_p
                          + sourceTermData.faceExtrapCoeff_a  * grad_a
                          + sourceTermData.faceExtrapCoeff_ib * grad_ib;

        return grad_f;
    }


} // end anonymous namespace


floatType ForceCalculator::GetViscousIntegralForce( const Axis::ENUMDATA forceComponent,
                                                    const IBCell::SourceTermData &sourceTermData,
                                                    const TensorIndex3D &cellIndex ) const
{
    using FVT::G;

    // Gradient in face normal direction
    floatType u2 = sourceTermData.faceValues.U[forceComponent],
              u1 = m_fields.U[forceComponent]( G( cellIndex ) ),
              dx = m_mesh.cellLengths[sourceTermData.direction]( cellIndex[sourceTermData.direction] ) / 2.0f;

    floatType grad1 =  static_cast<floatType>( sourceTermData.directionIndex )  // Direction index to account for whether face is hi or lo side of cell
                    * ( u2 - u1 ) / dx;        

    // Gradient in force component direction
    floatType grad2{0.0f};
    if ( forceComponent == sourceTermData.direction ) {
        grad2 = grad1;
    } else {
        floatType grad2_p = CalculateVelocityGradient( forceComponent, cellIndex                 , sourceTermData.direction, m_fields.U, m_ibData, m_mesh );
        floatType grad2_a = CalculateVelocityGradient( forceComponent, sourceTermData.cellIndex_a, sourceTermData.direction, m_fields.U, m_ibData, m_mesh );
        grad2 = ExtrapolateVelocityGradient( grad2_p, grad2_a, sourceTermData );
    }

    return - m_rho * m_nu * ( grad1 + grad2 ) * sourceTermData.faceAreaComponent;
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




}   // end namespace CFD

