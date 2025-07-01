#include "TurbulenceModels.h"

#include "../Core/Types.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Geometry/Geometry.h"
#include "TurbulenceModelTools.h"

#include <algorithm>
#include <cmath>

#include <iostream>

// Liu, J., Srebric, J., & Yu, N. (2013, July). A rapid and reliable numerical method for predictions of outdoor thermal environment in actual 
// urban areas. In Heat Transfer Summer Conference (Vol. 55492, p. V003T21A008). American Society of Mechanical Engineers.

namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::ZEQ3>::SetTurbulenceModelData( const InputData &inputData,
                                                                      const Mesh &mesh,
                                                   [[ maybe_unused ]] const IBData &ibData,
                                                                      const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    m_heightAxis                   = inputData.zeq3ModelData.heightAxis;
    m_averageBuildingHeight        = inputData.zeq3ModelData.averageBuildingHeight;
    m_inflowVelocityBuildingHeight = inputData.zeq3ModelData.inflowVelocityBuildingHeight;
    m_roughnessLength              = inputData.zeq3ModelData.roughnessLength;
    m_alpha                        = 3.15f;
    m_zh                           = 0.1f * m_averageBuildingHeight;
    m_b                            = 1.75f;

    // Length scale, distance to nearest wall
    Polyhedron geometry = MakeGeometry( inputData );
    Tree tree = MakeAABBTree( geometry );
    m_wallDistance = NearestWallDistance( mesh, tree, bcData );
}



void TurbulenceModel<TurbulenceModels::ZEQ3>::SetTurbulenceViscosityField( EnumVector<Axis, Tensor3D> &nuTurbulent,
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

                    const floatType nuIn  = m_alpha
                                          * m_zh
                                          * exp( - m_b * m_wallDistance[faceNormal](faceIndex) / m_averageBuildingHeight )
                                          * faceVelocityMagnitude
                                          * std::pow( m_wallDistance[faceNormal](faceIndex) / m_averageBuildingHeight, 2.0f );

                    const floatType z = ( m_heightAxis == faceNormal ) ? mesh.cellFaces[m_heightAxis]( faceIndex[m_heightAxis] )
                                                                       : mesh.cellCenters[m_heightAxis]( faceIndex[m_heightAxis] );

                    const floatType nuOut = ( z == 0 ) ? 0.0f       //  Avoid division by zero
                                                       : 0.16f
                                                       * m_inflowVelocityBuildingHeight
                                                       * ( z + m_roughnessLength ) / std::log( ( z + m_roughnessLength ) / m_roughnessLength );  

                    nuTurbulent[faceNormal](faceIndex) = std::max( nuIn, nuOut );
                    
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
                const floatType nuIn  = m_alpha
                                      * m_zh
                                      * exp( - m_b * m_wallDistance[faceNormal](faceIndex) / m_averageBuildingHeight )
                                      * faceVelocityMagnitude
                                      * std::pow( m_wallDistance[faceNormal](faceIndex) / m_averageBuildingHeight, 2.0f );

                const floatType z = ( m_heightAxis == faceNormal ) ? mesh.cellFaces[m_heightAxis]( faceIndex[m_heightAxis] )
                                                                   : mesh.cellCenters[m_heightAxis]( faceIndex[m_heightAxis] );

                const floatType nuOut = ( z == 0 ) ? 0.0f       //  Avoid division by zero
                                                   : 0.16f
                                                   * m_inflowVelocityBuildingHeight
                                                   * ( z + m_roughnessLength ) / std::log( ( z + m_roughnessLength ) / m_roughnessLength );  

                nuTurbulent[faceNormal](faceIndex) = std::max( nuIn, nuOut );

            }
        }
    }


}

}   // end namespace CAMIRA