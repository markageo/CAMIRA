#include "TurbulenceModels.h"

#include "../Core/Types.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Geometry/Geometry.h"
#include "TurbulenceModelTools.h"

#include <cmath>

// Davidovic, D. (2010). Improvements in Numerical Airflow Modeling Around Multiple Buildings (Doctoral dissertation, Pennsylvania State University).

namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::ZEQ2>::SetTurbulenceModelData( const InputData &inputData,
                                                                      const Mesh &mesh,
                                                   [[ maybe_unused ]] const IBData &ibData,
                                                                      const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    m_eddyViscosityRelaxation                = inputData.eddyViscosityRelaxation;

    m_averageBuildingHeight                  = inputData.zeq2ModelData.averageBuildingHeight;
    m_inflowVelocityBuildingHeight           = inputData.zeq2ModelData.inflowVelocityBuildingHeight;
    m_inflowTKEBuildingHeight                = inputData.zeq2ModelData.inflowTKEBuildingHeight;
    m_inflowIntergralTimeScaleBuildingHeight = inputData.zeq2ModelData.inflowIntergralTimeScaleBuildingHeight;
    m_wallDistanceLengthScale                = 1.0f;
    m_nu                                     = inputData.nu;

    // Length scale, distance to nearest wall
    Polyhedron geometry = MakeGeometry( inputData );
    Tree tree = MakeAABBTree( geometry );
    m_wallDistance = NearestWallDistance( mesh, tree, bcData );
}



void TurbulenceModel<TurbulenceModels::ZEQ2>::SetTurbulenceViscosityField( Tensor3D &nuTurbulent,
                                                                           const FieldData<Tensor3D> &fields,
                                                        [[ maybe_unused ]] const IBData &ibData,
                                                                           const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    // Includes ghost cells
    for ( intType k = -1; k != mesh.nCells[Z] + 1; k++ ) {
        for ( intType j = -1; j != mesh.nCells[Y] + 1; j++ ) {
            for ( intType i = -1; i != mesh.nCells[X] + 1; i++ ) {

                const TensorIndex3D cellIndexG = G(i, j, k);

                const floatType velocityMagnitude = sqrt( std::pow(fields.U[X](cellIndexG), 2.0f) 
                                                        + std::pow(fields.U[Y](cellIndexG), 2.0f) 
                                                        + std::pow(fields.U[Z](cellIndexG), 2.0f) );

                const floatType aReTurb = std::pow( ( m_inflowTKEBuildingHeight * m_inflowIntergralTimeScaleBuildingHeight ) / m_nu, 1.0f/3.0f );
                const floatType bReBulk = std::pow( ( m_inflowVelocityBuildingHeight * m_averageBuildingHeight ) / m_nu            , 1.0f/3.0f );

                const floatType wallDistanceNormalised = m_wallDistance(cellIndexG) / m_wallDistanceLengthScale;

                const floatType nuTurbulentNew = aReTurb 
                                               * wallDistanceNormalised
                                               * exp( -bReBulk * wallDistanceNormalised )
                                               * velocityMagnitude
                                               * m_wallDistance(cellIndexG);
                                                       

                nuTurbulent(cellIndexG) = (1.0f - m_eddyViscosityRelaxation ) * nuTurbulent(cellIndexG)
                                        + m_eddyViscosityRelaxation * nuTurbulentNew;
                
            }
        }
    }

}

}   // end namespace CAMIRA