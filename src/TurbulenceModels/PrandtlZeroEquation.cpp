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

    // Constant of proportionality
    m_vonKarmanConstant = 0.4;

    // Length scale, distance to nearest wall
    Polyhedron geometry = MakeGeometry( inputData );
    Tree tree = MakeAABBTree( geometry );
    m_wallDistance = NearestWallDistance( mesh, tree, bcData );

    m_velocityDeformationRate = Tensor3D( mesh.nCells[X] + 2*nGhost, mesh.nCells[Y] + 2*nGhost, mesh.nCells[Z] + 2*nGhost ).setZero();
}



void TurbulenceModel<TurbulenceModels::PrandtlZeroEquation>::SetTurbulenceViscosityField( EnumVector<Axis, Tensor3D> &nuTurbulent,
                                                                                          const FieldData<Tensor3D> &fields,
                                                                         [[maybe_unused]] const IBData &ibData,
                                                                                          const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    CalculateVelocityDeformationRate( m_velocityDeformationRate, fields, ibData, mesh );

    EnumFor<Axis>( [&] (Axis::ENUMDATA faceNormal) {

        for ( intType k = 0; k != mesh.nFacesNormal[faceNormal][Z]; k++ ) {
            for ( intType j = 0; j != mesh.nFacesNormal[faceNormal][Y]; j++ ) {
                for ( intType i = 0; i != mesh.nFacesNormal[faceNormal][X]; i++ ) {

                    TensorIndex3D faceIndex    = { i, j, k },
                                  loCellIndexG = G( i, j, k ),
                                  hiCellIndexG = G( i, j, k );
                    loCellIndexG[faceNormal] -= 1;
                    
                    const floatType lambda = mesh.interpFactors[faceNormal]( faceIndex[faceNormal] );
                    const floatType Sface = (1.0f - lambda) * m_velocityDeformationRate(loCellIndexG)  +  lambda * m_velocityDeformationRate(hiCellIndexG);
                    
                    const floatType lmix = m_vonKarmanConstant * m_wallDistance[faceNormal](faceIndex);

                    nuTurbulent[faceNormal](faceIndex) = lmix * lmix * Sface;
                    
                }
            }
        }

    } );


     // Go back through and correct for the immersed boundary faces
    for ( const auto &ibCellComponent : ibData.ibCells ) {
        for ( const auto &ibCell : ibCellComponent ) { 

            const TensorIndex3D &cellIndex = ibCell.cellIndex;

            for ( const auto &sourceTermData : ibCell.sourceTermsData ) {

                const Axis::ENUMDATA faceNormal = sourceTermData.direction;

                // Get face index
                TensorIndex3D faceIndex = cellIndex;
                faceIndex[faceNormal] += sourceTermData.faceDirectionIndex;

                // Cell face velocity deformation rate. Just take internal cell value
                const floatType Sface = m_velocityDeformationRate( G(cellIndex) );

                const floatType lmix = m_vonKarmanConstant * m_wallDistance[faceNormal](faceIndex);

                nuTurbulent[faceNormal](faceIndex) = lmix * lmix * Sface;

            }
        }
    }

}

}   // end namespace CAMIRA