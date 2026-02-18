#include "TurbulenceModels.h"

#include "Core/Types.h"
#include "Core/FVTools.h"
#include "Core/Geometry/Geometry.h"
#include "Core/Mesh/Mesh.h"
#include "Flow/FiniteVolume/FiniteVolume.h"
#include "TurbulenceModelTools.h"

// Chen, Q., & Xu, W. (1998). A zero-equation turbulence model for indoor airflow simulation. Energy and buildings, 28(2), 137-144.

namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::ZEQ0>::SetTurbulenceModelData( const InputData &inputData,
                                                                      const AxisTransformationMap &axisTransformation,
                                                                      const Mesh &mesh,
                                                   [[ maybe_unused ]] const IBData &ibData,
                                                                      const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    m_eddyViscosityRelaxation = inputData.eddyViscosityRelaxation;

    // Constant of proportionality
    m_proportionalityConstant = 0.03874;

    // Length scale, distance to nearest wall
    Polyhedron geometry = MakeGeometry( inputData.geometryData, axisTransformation );
    Tree tree = MakeAABBTree( geometry );
    NearestWallDistance( m_wallDistance, mesh, tree, bcData );
}



void TurbulenceModel<TurbulenceModels::ZEQ0>::SetTurbulenceViscosityField( Tensor3D &nuTurbulent,
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

                const floatType nuTurbulentNew = m_proportionalityConstant
                                               * m_wallDistance(cellIndexG)
                                               * velocityMagnitude;

                nuTurbulent(cellIndexG) = (1.0f - m_eddyViscosityRelaxation ) * nuTurbulent(cellIndexG)
                                        + m_eddyViscosityRelaxation * nuTurbulentNew;
                
            }
        }
    }


}

}   // end namespace CAMIRA