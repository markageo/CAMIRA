#include "TurbulenceModels.h"

#include "../Core/Types.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Geometry/Geometry.h"
#include "TurbulenceModelTools.h"


namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::PrandtlZeroEquation>::SetTurbulenceModelData( const Mesh &mesh,
                                                                                     const Polyhedron &geometry,
                                                                                     const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    // Constant of proportionality
    m_vonKarmanConstant = 0.4;

    // Length scale, distance to nearest wall
    Tree tree = MakeAABBTree( geometry );
    m_wallDistance = NearestWallDistance( mesh, tree, bcData );

    m_velocityDeformationRate = Tensor3D( mesh.nCells[X], mesh.nCells[Y], mesh.nCells[Z] ).setZero();
}



void TurbulenceModel<TurbulenceModels::PrandtlZeroEquation>::SetTurbulenceViscosityField( EnumVector<Axis, Tensor3D> &nuTurbulent,
                                                                                          const FieldData<Tensor3D> &fields,
                                                                                          const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    CalculateVelocityDeformationRate( m_velocityDeformationRate, fields, mesh );

    EnumFor<Axis>( [&] (Axis::ENUMDATA faceNormal) {

        for ( intType k = 0; k != mesh.nFacesNormal[faceNormal][Z]; k++ ) {
            for ( intType j = 0; j != mesh.nFacesNormal[faceNormal][Y]; j++ ) {
                for ( intType i = 0; i != mesh.nFacesNormal[faceNormal][X]; i++ ) {

                    TensorIndex3D faceIndex   = { i, j, k },
                                  loCellIndex = { i, j, k },
                                  hiCellIndex = { i, j, k };
                    loCellIndex[faceNormal] -= 1;
                    
                    floatType lambda = mesh.interpFactors[faceNormal]( faceIndex[faceNormal] );
                    floatType Sface = 0.0f;
                    if        ( faceIndex[faceNormal] == 0 ) {
                        Sface = m_velocityDeformationRate(hiCellIndex);
                    } else if ( faceIndex[faceNormal] == mesh.nFacesNormal[faceNormal][faceNormal]-1 ) {
                        Sface = m_velocityDeformationRate(loCellIndex);
                    } else {
                        Sface = (1.0f - lambda) * m_velocityDeformationRate(loCellIndex)  +  lambda * m_velocityDeformationRate(hiCellIndex);
                    }
                    
                    floatType lmix = m_vonKarmanConstant * m_wallDistance[faceNormal](faceIndex);

                    nuTurbulent[faceNormal](faceIndex) = lmix * lmix * Sface;
                    
                }
            }
        }

    } );
}

}   // end namespace CAMIRA