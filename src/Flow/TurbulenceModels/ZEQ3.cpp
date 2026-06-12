#include "TurbulenceModels.h"

#include "Core/Types.h"
#include "Core/FVTools.h"
#include "Core/Geometry/Geometry.h"
#include "Core/Mesh/Mesh.h"
#include "Flow/FiniteVolume/FiniteVolume.h"
#include "TurbulenceModelTools.h"

#include <algorithm>
#include <cmath>

#include <iostream>

// Liu, J., Srebric, J., & Yu, N. (2013, July). A rapid and reliable numerical method for predictions of outdoor thermal environment in actual 
// urban areas. In Heat Transfer Summer Conference (Vol. 55492, p. V003T21A008). American Society of Mechanical Engineers.

namespace CAMIRA
{

using namespace CORE;

namespace FLOW
{

void TurbulenceModel<TurbulenceModels::ZEQ3>::SetTurbulenceModelData( const InputData &inputData,
                                                                      const AxisTransformationMap &axisTransformation,
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
    Polyhedron P; 
    MakePolyhedron( P, inputData.geometryData, axisTransformation );
    Tree tree;
    MakeAABBTree( tree, P );
    
    NearestWallDistance( m_wallDistance, mesh, tree, bcData );
}



void TurbulenceModel<TurbulenceModels::ZEQ3>::SetTurbulenceViscosityField( Tensor3D &nuTurbulent,
                                                                           const FieldData<Tensor3D> &fields,
                                                        [[ maybe_unused ]] const IBData &ibData,
                                                                           const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    // Includes ghost cells
    #pragma omp parallel for collapse(3)
    for ( intType k = -1; k != mesh.nCells[Z] + 1; k++ ) {
        for ( intType j = -1; j != mesh.nCells[Y] + 1; j++ ) {
            for ( intType i = -1; i != mesh.nCells[X] + 1; i++ ) {

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
  
                const floatType z = GetVerticalHeight( mesh, cellIndex, m_heightAxis );

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

}   // end namespace FLOW

}   // end namespace CAMIRA