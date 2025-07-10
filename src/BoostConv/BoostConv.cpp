#include "BoostConv.h"

#include "../Core/FVTools.h"

#include <utility>
#include <iostream>

namespace CAMIRA
{


namespace 
{


floatType FieldDotProduct( const Tensor3D &f1,
                           const Tensor3D &f2 )
{   
    using FVT::G;
    floatType result = 0.0f;

    // Exclude ghost cells
    const intType ni = f1.dimension(0) - 2*nGhost,
                  nj = f1.dimension(1) - 2*nGhost,
                  nk = f1.dimension(2) - 2*nGhost;

    for ( intType k = 0; k != nk; k++ ) {
        for ( intType j = 0; j != nj; j++ ) {
            for ( intType i = 0; i != ni; i++ ) {
                result += f1(G(i, j, k)) * f2(G(i, j, k));
            }
        }
    }

    return result;
}


}   // end anonymous namespace



BoostConv::BoostConv( const intType basisSize, 
                      const intType startIteration,
                      const floatType relaxation,
                      const FieldData<Tensor3D> &initialResidual ) :
                        m_N( basisSize ),
                        m_startIteration( startIteration ),
                        m_relaxation( relaxation ),
                        m_previousResidual( initialResidual ),
                        m_iteration(0),
                        stateModification( Tensor3D( initialResidual[0].dimensions() ).setZero() )
{
    m_v.reserve( m_N );
    m_w.reserve( m_N );
}



void BoostConv::UpdateStateModification( const FieldData<Tensor3D> &currentResidual )
{

    // No correction for first few iterations
    m_iteration++;
    if ( m_iteration < m_startIteration ) {
         m_previousResidual = currentResidual;
         return;
    }

    // Grow the basis dimension until we reach N
    intType dim = m_D[0].rows();
    if ( dim < m_N ) {

        // Increase the dimension of everything 
        dim += 1;
        
        // Add new set of vectors
        m_v.emplace_back();
        m_w.emplace_back();
            
        ForAllFieldData( [&] (intType f) {
            // Resize least squares matrix and known term
            m_D[f].conservativeResize( dim, dim );
            m_t[f].conservativeResize( dim );
        } );

    } else {

        // Discard the oldest vectors
        for ( size_t m = 0; m != static_cast<size_t>(m_N)-1; m++ ) {
            m_v[m] = m_v[m+1];
            m_w[m] = m_w[m+1];
        }
        
        // Discard old elements of the least squares matrix
        ForAllFieldData( [&] (intType f) {
            for ( intType i = 0; i != dim-1; i++ ) {
                for ( intType j = 0; j != dim-1; j++ ) {
                    m_D[f](i, j) = m_D[f](i+1, j+1);
                }
            }
        } );

    }

    // Update the newest vectors
    ForAllFieldData( [&] (intType f) {
        m_v[dim-1][f] = currentResidual[f] - m_previousResidual[f];
        m_w[dim-1][f] = stateModification[f] - m_v[dim-1][f];
    } );

    // For next iteration
    m_previousResidual = currentResidual;

    
    ForAllFieldData( [&] (intType f) {

        // Update the least squares matrix
        for ( intType i = 0; i != dim; i++ ) {

            // Matrix element
            m_D[f](i, dim-1) = FieldDotProduct( m_v[dim-1][f], m_v[i][f] );

            // Matrix is symmetric
            m_D[f](dim-1, i) = m_D[f](i, dim-1);

        }

        // Build the known term
        for ( intType i = 0; i != dim; i++ ) {
            m_t[f](i) = FieldDotProduct( m_v[i][f], currentResidual[f] );
        }

        // Solve the least squares problem
        m_c[f] = m_D[f].fullPivLu().solve( m_t[f] );


        // Set the state modification
        stateModification[f]  = m_w[0][f]            // First one is assigned
                              * m_w[0][f].constant( m_relaxation * m_c[f](0) );

        for ( intType i = 1; i != dim; i++ ) {
            stateModification[f] += m_w[i][f]
                                  * m_w[i][f].constant( m_relaxation * m_c[f](i) );
        }

    } );

}



} // end namespace CAMIRA