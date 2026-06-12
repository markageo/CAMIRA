#include "TurbulenceModels.h"

#include "Core/Types.h"
#include "Core/FVTools.h"
#include "Core/Geometry/Geometry.h"
#include "Core/Mesh/Mesh.h"
#include "Flow/FiniteVolume/FiniteVolume.h"
#include "TurbulenceModelTools.h"

#include <algorithm>
#include <cmath>


// Li, C., Li, X., Su, Y., & Zhu, Y. (2012). A new zero-equation turbulence model for micro-scale climate simulation. Building and Environment, 47, 243-255.

namespace CAMIRA
{

using namespace CORE;

namespace FLOW
{

void TurbulenceModel<TurbulenceModels::ZEQ4>::SetTurbulenceModelData( const InputData &inputData,
                                                                      const AxisTransformationMap &axisTransformation,
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
    Polyhedron P; 
    MakePolyhedron( P, inputData.geometryData, axisTransformation );
    Tree tree;
    MakeAABBTree( tree, P );
    
    NearestWallDistance( m_wallDistance, mesh, tree, bcData );

    m_velocityDeformationRate = Tensor3D( mesh.nCells[X] + 2*nGhost, mesh.nCells[Y] + 2*nGhost, mesh.nCells[Z] + 2*nGhost );
    SetTensorZeroParallel( m_velocityDeformationRate );
}



void TurbulenceModel<TurbulenceModels::ZEQ4>::SetTurbulenceViscosityField( Tensor3D &nuTurbulent,
                                                                           const FieldData<Tensor3D> &fields,
                                                        [[ maybe_unused ]] const IBData &ibData,
                                                                           const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    CalculateVelocityDeformationRate( m_velocityDeformationRate, fields, ibData, mesh );

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

                const floatType z = GetVerticalHeight( mesh, cellIndex, m_heightAxis );

                const floatType nuIn  = m_velocityDeformationRate(cellIndexG) 
                                      * std::pow( 
                                                1.8f
                                            * ( 1.0f - exp( -0.645 * std::pow( m_averageBuildingWidth / m_averageBuildingHeight, 0.8 ) ) )
                                            * ( -2.0f * std::min( z / m_averageBuildingHeight , static_cast<floatType>( 1.0 ) ) )
                                            * m_wallDistance(cellIndexG)
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

}   // end namespace FLOW

}   // end namespace CAMIRA