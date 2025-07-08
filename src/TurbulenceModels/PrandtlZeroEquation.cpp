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
                                                                                     const Mesh &mesh,
                                                                  [[ maybe_unused ]] const IBData &ibData,
                                                                                     const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    m_eddyViscosityRelaxation = inputData.eddyViscosityRelaxation;

    // Constant of proportionality
    m_vonKarmanConstant = 0.4;

    // Length scale, distance to nearest wall
    Polyhedron geometry = MakeGeometry( inputData );
    Tree tree = MakeAABBTree( geometry );
    m_wallDistance = NearestWallDistance( mesh, tree, bcData );

    m_velocityDeformationRate = Tensor3D( mesh.nCells[X] + 2*nGhost, mesh.nCells[Y] + 2*nGhost, mesh.nCells[Z] + 2*nGhost ).setZero();
}



void TurbulenceModel<TurbulenceModels::PrandtlZeroEquation>::SetTurbulenceViscosityField( Tensor3D &nuTurbulent,
                                                                                          const FieldData<Tensor3D> &fields,
                                                                         [[maybe_unused]] const IBData &ibData,
                                                                                          const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    CalculateVelocityDeformationRate( m_velocityDeformationRate, fields, ibData, mesh );

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

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