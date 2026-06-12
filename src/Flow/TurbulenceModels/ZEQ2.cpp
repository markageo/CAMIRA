#include "TurbulenceModels.h"

#include "Core/Types.h"
#include "Core/FVTools.h"
#include "Core/Geometry/Geometry.h"
#include "Core/Mesh/Mesh.h"
#include "Flow/FiniteVolume/FiniteVolume.h"
#include "TurbulenceModelTools.h"

#include <cmath>

// Davidovic, D. (2010). Improvements in Numerical Airflow Modeling Around Multiple Buildings (Doctoral dissertation, Pennsylvania State University).

namespace CAMIRA
{

using namespace CORE;

namespace FLOW
{

void TurbulenceModel<TurbulenceModels::ZEQ2>::SetTurbulenceModelData( const InputData &inputData,
                                                                      const AxisTransformationMap &axisTransformation,
                                                                      const Mesh &mesh,
                                                   [[ maybe_unused ]] const IBData &ibData,
                                                                      const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    m_eddyViscosityRelaxation                = inputData.eddyViscosityRelaxation;

    m_heightAxis                             = inputData.zeq2ModelData.heightAxis;
    m_averageBuildingHeight                  = inputData.zeq2ModelData.averageBuildingHeight;
    m_inflowVelocityBuildingHeight           = inputData.zeq2ModelData.inflowVelocityBuildingHeight;
    m_inflowTKEBuildingHeight                = inputData.zeq2ModelData.inflowTKEBuildingHeight;
    m_inflowIntergralTimeScaleBuildingHeight = inputData.zeq2ModelData.inflowIntergralTimeScaleBuildingHeight;
    m_roughnessLength                        = inputData.zeq2ModelData.roughnessLength;
    m_wallDistanceLengthScale                = 1.0f;
    m_nu                                     = inputData.nu;

    // Length scale, distance to nearest wall
    Polyhedron P; 
    MakePolyhedron( P, inputData.geometryData, axisTransformation );
    Tree tree;
    MakeAABBTree( tree, P );
    
    NearestWallDistance( m_wallDistance, mesh, tree, bcData );
}



void TurbulenceModel<TurbulenceModels::ZEQ2>::SetTurbulenceViscosityField( Tensor3D &nuTurbulent,
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

                const floatType aReTurb = std::pow( ( m_inflowTKEBuildingHeight * m_inflowIntergralTimeScaleBuildingHeight ) / m_nu, 1.0f/3.0f );
                const floatType bReBulk = std::pow( ( m_inflowVelocityBuildingHeight * m_averageBuildingHeight ) / m_nu            , 1.0f/3.0f );

                const floatType wallDistanceNormalised = m_wallDistance(cellIndexG) / m_wallDistanceLengthScale;

                const floatType nuIn = aReTurb 
                                     * wallDistanceNormalised
                                     * exp( -bReBulk * wallDistanceNormalised )
                                     * velocityMagnitude
                                     * m_wallDistance(cellIndexG);

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