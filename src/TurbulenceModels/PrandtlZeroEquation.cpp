#include "TurbulenceModels.h"

#include "../Core/Types.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Geometry/Geometry.h"
#include "TurbulenceModelTools.h"


namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::PrandtlZeroEquation>::SetTurbulenceModelData( const InputData &inputData, 
                                                                                     const AxisTransformationMap &axisTransformation,
                                                                                     const Mesh &mesh,
                                                                  [[ maybe_unused ]] const IBData &ibData,
                                                                                     const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    m_eddyViscosityRelaxation = inputData.eddyViscosityRelaxation;

    // Constant of proportionality
    m_vonKarmanConstant = 0.4;

    // Length scale, distance to nearest wall
    Polyhedron geometry = MakeGeometry( inputData, axisTransformation );
    Tree tree = MakeAABBTree( geometry );
    NearestWallDistance( m_wallDistance, mesh, tree, bcData );
    
    m_velocityDeformationRate = Tensor3D( mesh.nCells[X] + 2*nGhost, mesh.nCells[Y] + 2*nGhost, mesh.nCells[Z] + 2*nGhost );
    SetTensorZeroParallel( m_velocityDeformationRate );
}



void TurbulenceModel<TurbulenceModels::PrandtlZeroEquation>::SetTurbulenceViscosityField( Tensor3D &nuTurbulent,
                                                                                          const FieldData<Tensor3D> &fields,
                                                                         [[maybe_unused]] const IBData &ibData,
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

                const TensorIndex3D cellIndexG = G(i, j, k);

                const floatType lmix = m_vonKarmanConstant * m_wallDistance(cellIndexG);

                const floatType nuTurbulentNew = lmix * lmix * m_velocityDeformationRate(cellIndexG);

                nuTurbulent(cellIndexG) = (1.0f - m_eddyViscosityRelaxation ) * nuTurbulent(cellIndexG)
                                        + m_eddyViscosityRelaxation * nuTurbulentNew;
                
            }
        }
    }

}

}   // end namespace CAMIRA