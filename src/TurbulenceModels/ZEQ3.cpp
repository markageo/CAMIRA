#include "TurbulenceModels.h"

#include "../Core/Types.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Geometry/Geometry.h"
#include "TurbulenceModelTools.h"

#include <algorithm>
#include <cmath>

#include <iostream>

// Liu, J., Srebric, J., & Yu, N. (2013, July). A rapid and reliable numerical method for predictions of outdoor thermal environment in actual 
// urban areas. In Heat Transfer Summer Conference (Vol. 55492, p. V003T21A008). American Society of Mechanical Engineers.

namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::ZEQ3>::SetTurbulenceModelData( const InputData &inputData,
                                                                      const Mesh &mesh,
                                                   [[ maybe_unused ]] const IBData &ibData,
                                                                      const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    m_eddyViscosityRelaxation      = inputData.eddyViscosityRelaxation;

    m_heightAxis                   = inputData.zeq3ModelData.heightAxis;
    m_averageBuildingHeight        = inputData.zeq3ModelData.averageBuildingHeight;
    m_inflowVelocityBuildingHeight = inputData.zeq3ModelData.inflowVelocityBuildingHeight;
    m_roughnessLength              = inputData.zeq3ModelData.roughnessLength;
    m_alpha                        = 3.15f;
    m_zh                           = 0.1f * m_averageBuildingHeight;
    m_b                            = 1.75f;

    // Length scale, distance to nearest wall
    Polyhedron geometry = MakeGeometry( inputData );
    Tree tree = MakeAABBTree( geometry );
    m_wallDistance = NearestWallDistance( mesh, tree, bcData );
}



void TurbulenceModel<TurbulenceModels::ZEQ3>::SetTurbulenceViscosityField( Tensor3D &nuTurbulent,
                                                                           const FieldData<Tensor3D> &fields,
                                                                           const IBData &ibData,
                                                                           const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                const TensorIndex3D cellIndex  = {i, j, k};
                const TensorIndex3D cellIndexG = G(i, j, k);

                const floatType velocityMagnitude = sqrt( std::pow(fields.U[X](cellIndexG), 2.0f) 
                                                        + std::pow(fields.U[Y](cellIndexG), 2.0f) 
                                                        + std::pow(fields.U[Z](cellIndexG), 2.0f) );

                const floatType nuIn  = m_alpha
                                      * m_zh
                                      * exp( - m_b * m_wallDistance(cellIndexG) / m_averageBuildingHeight )
                                      * velocityMagnitude
                                      * std::pow( m_wallDistance(cellIndexG) / m_averageBuildingHeight, 2.0f );

                const floatType z = mesh.cellCenters[m_heightAxis]( cellIndex[m_heightAxis] );

                const floatType nuOut = ( z == 0 ) ? 0.0f       //  Avoid division by zero
                                                   : 0.16f
                                                   * m_inflowVelocityBuildingHeight
                                                   * ( z + m_roughnessLength ) / std::log( ( z + m_roughnessLength ) / m_roughnessLength );  

                const floatType nuTurbulentNew = std::max( nuIn, nuOut );

                nuTurbulent(cellIndexG) = (1.0f - m_eddyViscosityRelaxation ) * nuTurbulent(cellIndexG)
                                        + m_eddyViscosityRelaxation * nuTurbulentNew;
                
            }
        }
    }

}

}   // end namespace CAMIRA