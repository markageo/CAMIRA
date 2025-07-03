#ifndef CAMIRA_BOOST_CONV
#define CAMIRA_BOOST_CONV

#include "../Core/Types.h"
#include <Eigen/Dense>

#include <vector>


namespace CAMIRA
{


class BoostConv
{

private:

    intType m_N;                                                    // Number of solutions to store
    floatType m_relaxation;                                         // Relaxation on modification to residual  
    FieldData<Tensor3D> m_previousResidual;                         // Residual from previous iteration                                                 
    std::vector<FieldData<Tensor3D>> m_v, m_w;                      // Vector pairs
    Eigen::Matrix<floatType, Eigen::Dynamic, Eigen::Dynamic> m_D;   // m_N x m_N
    Eigen::Matrix<floatType, Eigen::Dynamic, 1>  m_t, m_c;          // m_N x 1   
    bool m_isFirstIteration;                    

public:

    FieldData<Tensor3D> stateModification;                          // xi in the mathematical notation

    BoostConv( const intType, const floatType, const FieldData<Tensor3D> & );

    void UpdateStateModification( const FieldData<Tensor3D> & );

};



} // end namespace CAMIRA

#endif // CAMIRA_BOOST_CONV