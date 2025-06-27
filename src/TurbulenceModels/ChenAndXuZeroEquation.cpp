#include "TurbulenceModels.h"

#include "../Core/Types.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Geometry/Geometry.h"
#include "TurbulenceModelTools.h"


// DEBUGGING
#include <memory>
#include "../IO/VTKWriter.h"
// END DEBUGGING

namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::ChenAndXuZeroEquation>::SetTurbulenceModelData( const Mesh &mesh,
                                                                                       const Polyhedron &geometry,
                                                                                       const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    // Constant of proportionality
    m_proportionalityConstant = 0.03874;

    // Length scale, distance to nearest wall
    Tree tree = MakeAABBTree( geometry );
    m_wallDistance = NearestWallDistance( mesh, tree, bcData );
}



void TurbulenceModel<TurbulenceModels::ChenAndXuZeroEquation>::SetTurbulenceViscosityField( EnumVector<Axis, Tensor3D> &nuTurbulent,
                                                                                            const FieldData<Tensor3D> &fields,
                                                                                            const IBData &ibData,
                                                                                            const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    EnumFor<Axis>( [&] (Axis::ENUMDATA faceNormal) {

        for ( intType k = 0; k != mesh.nFacesNormal[faceNormal][Z]; k++ ) {
            for ( intType j = 0; j != mesh.nFacesNormal[faceNormal][Y]; j++ ) {
                for ( intType i = 0; i != mesh.nFacesNormal[faceNormal][X]; i++ ) {

                    TensorIndex3D faceIndex    = {i, j, k},
                                  loCellIndexG = G( i, j, k ),
                                  hiCellIndexG = G( i, j, k );
                    loCellIndexG[faceNormal] -= 1;
                    
                    // Domain boundaries evaluated using ghost cells
                    floatType lambda = mesh.interpFactors[faceNormal]( faceIndex[faceNormal] );  
                    
                    // const floatType faceVelocityMagnitude = sqrt(
                    //                                               std::pow( fields.U[X](loCellIndexG) * ( 1.0f - lambda )  +  fields.U[X](hiCellIndexG) * lambda , 2.0f )
                    //                                             + std::pow( fields.U[Y](loCellIndexG) * ( 1.0f - lambda )  +  fields.U[Y](hiCellIndexG) * lambda , 2.0f )
                    //                                             + std::pow( fields.U[Z](loCellIndexG) * ( 1.0f - lambda )  +  fields.U[Z](hiCellIndexG) * lambda , 2.0f ) 
                    //                                         );


                    const floatType hiVelocityMagnitude = sqrt( std::pow(fields.U[X](hiCellIndexG), 2.0f) + std::pow(fields.U[Y](hiCellIndexG), 2.0f) + std::pow(fields.U[Z](hiCellIndexG), 2.0f) );
                    const floatType loVelocityMagnitude = sqrt( std::pow(fields.U[X](loCellIndexG), 2.0f) + std::pow(fields.U[Y](loCellIndexG), 2.0f) + std::pow(fields.U[Z](loCellIndexG), 2.0f) );
                    const floatType faceVelocityMagnitude = ( 1.0f - lambda ) * loVelocityMagnitude  +  lambda * hiVelocityMagnitude;

                    nuTurbulent[faceNormal](faceIndex) = m_proportionalityConstant
                                                       * m_wallDistance[faceNormal](faceIndex)
                                                       * faceVelocityMagnitude;
                    
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

                // Velocity magnitude on IB face
                const floatType faceVelocityMagnitude = sqrt( 
                                                            + std::pow( sourceTermData.faceValues.U[X], 2.0f )
                                                            + std::pow( sourceTermData.faceValues.U[Y], 2.0f )
                                                            + std::pow( sourceTermData.faceValues.U[Z], 2.0f ) 
                                                        );

                // Modify the turbulent viscosity
                nuTurbulent[faceNormal](faceIndex) = m_proportionalityConstant
                                                   * m_wallDistance[faceNormal](faceIndex)
                                                   * faceVelocityMagnitude;


            }
        }
    }


}

}   // end namespace CAMIRA