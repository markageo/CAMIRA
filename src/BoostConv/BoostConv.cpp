#include "BoostConv.h"

#include "../Core/FVTools.h"

#include <utility>
#include <iostream>

namespace CAMIRA
{


namespace 
{


floatType FieldDataDotProduct( const FieldData<Tensor3D> &f1,
                               const FieldData<Tensor3D> &f2 )
{   
    using FVT::G;
    floatType result = 0.0f;

    // Exclude ghost cells
    const intType ni = f1[0].dimension(0) - 2*nGhost,
                  nj = f1[0].dimension(1) - 2*nGhost,
                  nk = f1[0].dimension(2) - 2*nGhost;

    // ForAllFieldData( [&] (intType f) {
    //     for ( intType k = 0; k != nk; k++ ) {
    //         for ( intType j = 0; j != nj; j++ ) {
    //             for ( intType i = 0; i != ni; i++ ) {
    //                 result += f1[f](G(i, j, k)) * f2[f](G(i, j, k));
    //             }
    //         }
    //     }
    // } );

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        for ( intType k = 0; k != nk; k++ ) {
            for ( intType j = 0; j != nj; j++ ) {
                for ( intType i = 0; i != ni; i++ ) {
                    result += f1.U[axis](G(i, j, k)) * f2.U[axis](G(i, j, k));
                }
            }
        }
    } );

    return result;
}


}   // end anonymous namespace



BoostConv::BoostConv( const intType N, 
                      const floatType relaxation,
                      const FieldData<Tensor3D> &initialResidual ) :
                        m_N( N ),
                        m_relaxation( relaxation ),
                        m_previousResidual( initialResidual ),
                        m_isFirstIteration(true),
                        stateModification( Tensor3D( initialResidual[0].dimensions() ).setZero() )
{
    m_v.reserve( m_N );
    m_w.reserve( m_N );
}



void BoostConv::UpdateStateModification( const FieldData<Tensor3D> &currentResidual )
{

    // No correction after the first iteration
    if ( m_isFirstIteration ) {
        m_previousResidual = currentResidual;
        ForAllFieldData( [&] (intType f) {
            stateModification[f].setZero();
        } );
        m_isFirstIteration = false;
        return;
    }


    // Grow the basis dimension until we reach N
    intType dim = m_D.rows();
    if ( dim < m_N ) {

        // Increase the dimension of everything 
        dim += 1;
        
        // Add new set of vectors
        m_v.emplace_back();
        m_w.emplace_back();
        
        // Resize least squares matrix and known term
        m_D.conservativeResize( dim, dim );
        m_t.conservativeResize( dim );

    } else {

        // Discard the oldest vectors
        for ( size_t m = 0; m != static_cast<size_t>(m_N)-1; m++ ) {
            m_v[m] = m_v[m+1];
            m_w[m] = m_w[m+1];
        }
        
        // Discard old elements of the least squares matrix
        for ( intType i = 0; i != dim-1; i++ ) {
            for ( intType j = 0; j != dim-1; j++ ) {
                m_D(i, j) = m_D(i+1, j+1);
            }
        }

    }

    // Update the newest vectors
    ForAllFieldData( [&] (intType f) {
        // m_v[dim-1][f] = m_previousResidual[f] - currentResidual[f]; // Other way around??
        m_v[dim-1][f] = currentResidual[f] - m_previousResidual[f];
        m_w[dim-1][f] = stateModification[f] - m_v[dim-1][f];
    } );

    // For next iteration
    m_previousResidual = currentResidual;

    // Update the least squares matrix
    for ( intType i = 0; i != dim; i++ ) {

        // Matrix element
        m_D(i, dim-1) = FieldDataDotProduct( m_v[dim-1], m_v[i] );

        // Matrix is symmetric
        m_D(dim-1, i) = m_D(i, dim-1);

    }


    // Build the known term
    for ( intType i = 0; i != dim; i++ ) {
        m_t(i) = FieldDataDotProduct( m_v[i], currentResidual );
    }

    // Solve the least squares problem
    m_c = m_D.fullPivLu().solve( m_t );

    // std::cout << "t: " << std::endl << m_t << std::endl;

    // std::cout << "D: " << std::endl << m_D << std::endl;

    // std::cout << "D*c: " << std::endl << m_D * m_c << std::endl;

    // std::cout << "c: " << std::endl << m_c << std::endl;

    // // Compute the modification to the state
    // ForAllFieldData( [&] (intType f) {

    //     stateModification[f]  = m_w[0][f]            // First one is assigned
    //                           * m_w[0][f].constant( m_relaxation * m_c(0) );

    //     for ( intType i = 1; i != dim; i++ ) {
    //         stateModification[f] += m_w[i][f]
    //                               * m_w[i][f].constant( m_relaxation * m_c(i) );
    //     }

    // } );

    // Compute the modification to the state (velocity only)
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        stateModification.U[axis]  = m_w[0].U[axis]            // First one is assigned
                                   * m_w[0].U[axis].constant( m_relaxation * m_c(0) );

        for ( intType i = 1; i != dim; i++ ) {
            stateModification.U[axis] += m_w[i].U[axis]
                                       * m_w[i].U[axis].constant( m_relaxation * m_c(i) );
        }

    } );

}



} // end namespace CAMIRA