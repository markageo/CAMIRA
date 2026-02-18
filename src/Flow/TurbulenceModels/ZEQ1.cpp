#include "TurbulenceModels.h"

#include "Core/Types.h"
#include "Core/FVTools.h"
#include "Core/Geometry/Geometry.h"
#include "Core/Mesh/Mesh.h"
#include "Flow/FiniteVolume/FiniteVolume.h"
#include "TurbulenceModelTools.h"

// Qian, Y. (2004). Development of an Algebraic Turbulence Model for Airflow and Contaminant Simulations around a Building 
// (Doctoral dissertation, Pennsylvania State University).

namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::ZEQ1>::SetTurbulenceModelData( const InputData &inputData,
                                                                      const AxisTransformationMap &axisTransformation,
                                                                      const Mesh &mesh,
                                                   [[ maybe_unused ]] const IBData &ibData,
                                                                      const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    m_eddyViscosityRelaxation                 = inputData.eddyViscosityRelaxation;

    m_reynoldsNumberBuildingHeight            = inputData.zeq1ModelData.reynoldsNumberBuildingHeight;
    m_inflowTurbulenceIntensityBuildingHeight = inputData.zeq1ModelData.inflowTurbulenceIntensityBuildingHeight; 

    // Length scale, distance to nearest wall
    Polyhedron geometry = MakeGeometry( inputData.geometryData, axisTransformation );
    Tree tree = MakeAABBTree( geometry );
    NearestWallDistance( m_wallDistance, mesh, tree, bcData );
}



void TurbulenceModel<TurbulenceModels::ZEQ1>::SetTurbulenceViscosityField( Tensor3D &nuTurbulent,
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

                const TensorIndex3D cellIndexG = G(i, j, k);

                const floatType velocityMagnitude = sqrt( std::pow(fields.U[X](cellIndexG), 2.0f) 
                                                        + std::pow(fields.U[Y](cellIndexG), 2.0f) 
                                                        + std::pow(fields.U[Z](cellIndexG), 2.0f) );

                const floatType nuTurbulentNew = 0.2f 
                                               * ( 1.0e5 / m_reynoldsNumberBuildingHeight )
                                               * m_inflowTurbulenceIntensityBuildingHeight
                                               * velocityMagnitude
                                               * m_wallDistance(cellIndexG);

                nuTurbulent(cellIndexG) = (1.0f - m_eddyViscosityRelaxation ) * nuTurbulent(cellIndexG)
                                        + m_eddyViscosityRelaxation * nuTurbulentNew;
                
            }
        }
    }

}

}   // end namespace CAMIRA