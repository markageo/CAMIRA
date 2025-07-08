#include "TurbulenceModels.h"

#include "../Core/Types.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Geometry/Geometry.h"
#include "TurbulenceModelTools.h"

#include <algorithm>
#include <cmath>


// Li, C., Li, X., Su, Y., & Zhu, Y. (2012). A new zero-equation turbulence model for micro-scale climate simulation. Building and Environment, 47, 243-255.

namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::ZEQ4>::SetTurbulenceModelData( const InputData &inputData,
                                                                      const Mesh &mesh,
                                                   [[ maybe_unused ]] const IBData &ibData,
                                                                      const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    m_eddyViscosityRelaxation = inputData.eddyViscosityRelaxation;

    m_heightAxis              = inputData.zeq4ModelData.heightAxis;
    m_averageBuildingHeight   = inputData.zeq4ModelData.averageBuildingHeight;
    m_averageBuildingWidth    = inputData.zeq4ModelData.averageBuildingWidth;
    m_referenceHeight         = inputData.zeq4ModelData.referenceHeight;
    m_Cmu                     = 0.09f;
    m_Ig                      = 0.1f;
    m_alpha                   = 0.22;

    // Length scale, distance to nearest wall
    Polyhedron geometry = MakeGeometry( inputData );
    Tree tree = MakeAABBTree( geometry );
    m_wallDistance = NearestWallDistance( mesh, tree, bcData );

    m_velocityDeformationRate = Tensor3D( mesh.nCells[X] + 2*nGhost, mesh.nCells[Y] + 2*nGhost, mesh.nCells[Z] + 2*nGhost ).setZero();
}



void TurbulenceModel<TurbulenceModels::ZEQ4>::SetTurbulenceViscosityField( Tensor3D &nuTurbulent,
                                                                           const FieldData<Tensor3D> &fields,
                                                                           const IBData &ibData,
                                                                           const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    CalculateVelocityDeformationRate( m_velocityDeformationRate, fields, ibData, mesh );

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                const TensorIndex3D cellIndex  = {i, j, k};
                const TensorIndex3D cellIndexG = G(i, j, k);

                const floatType velocityMagnitude = sqrt( std::pow(fields.U[X](cellIndexG), 2.0f) 
                                                        + std::pow(fields.U[Y](cellIndexG), 2.0f) 
                                                        + std::pow(fields.U[Z](cellIndexG), 2.0f) );

                const floatType z = mesh.cellCenters[m_heightAxis]( cellIndex[m_heightAxis] );

                const floatType nuIn  = m_velocityDeformationRate(cellIndexG) 
                                      * std::pow( 
                                                1.8f
                                            * ( 1.0f - exp( -0.645 * std::pow( m_averageBuildingWidth / m_averageBuildingHeight, 0.8 ) ) )
                                            * ( -2.0f * std::min( z / m_averageBuildingHeight , static_cast<floatType>( 1.0 ) ) )
                                            * std::pow( m_wallDistance(cellIndexG), 2.0 )
                                        , 2.0f );

                const floatType nuOut = sqrt( m_Cmu )
                                      * std::pow( m_Ig, 2.0f )
                                      * std::pow( m_referenceHeight, 2.0f * m_alpha + 0.1 )
                                      * std::pow( z, 0.9 - 2 * m_alpha ) / m_alpha
                                      * velocityMagnitude
                                      * m_wallDistance(cellIndexG);

                const floatType nuTurbulentNew = std::max( nuIn, nuOut );


                nuTurbulent(cellIndexG) = (1.0f - m_eddyViscosityRelaxation ) * nuTurbulent(cellIndexG)
                                        + m_eddyViscosityRelaxation * nuTurbulentNew;
                
            }
        }
    }


}

}   // end namespace CAMIRA