#include "TurbulenceModels.h"

#include "../Core/Types.h"
#include "../Core/FVTools.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../Geometry/Geometry.h"
#include "TurbulenceModelTools.h"

#include <algorithm>
#include <cmath>


// Li, C., Li, X., Su, Y., & Zhu, Y. (2012). A new zero-equation turbulence model for micro-scale climate simulation. Building and Environment, 47, 243-255.

namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::ZEQ4>::SetTurbulenceModelData( const InputData &inputData,
                                                                      const Mesh &mesh,
                                                   [[ maybe_unused ]] const IBData &ibData,
                                                                      const BoundaryConditionData &bcData )
{ 
    using enum Axis::ENUMDATA;

    m_heightAxis            = inputData.zeq4ModelData.heightAxis;
    m_averageBuildingHeight = inputData.zeq4ModelData.averageBuildingHeight;
    m_averageBuildingWidth  = inputData.zeq4ModelData.averageBuildingWidth;
    m_referenceHeight       = inputData.zeq4ModelData.referenceHeight;
    m_Cmu                   = 0.09f;
    m_Ig                    = 0.1f;
    m_alpha                 = 0.22;

    // Length scale, distance to nearest wall
    Polyhedron geometry = MakeGeometry( inputData );
    Tree tree = MakeAABBTree( geometry );
    m_wallDistance = NearestWallDistance( mesh, tree, bcData );

    m_velocityDeformationRate = Tensor3D( mesh.nCells[X] + 2*nGhost, mesh.nCells[Y] + 2*nGhost, mesh.nCells[Z] + 2*nGhost ).setZero();
}



void TurbulenceModel<TurbulenceModels::ZEQ4>::SetTurbulenceViscosityField( EnumVector<Axis, Tensor3D> &nuTurbulent,
                                                                           const FieldData<Tensor3D> &fields,
                                                                           const IBData &ibData,
                                                                           const Mesh &mesh )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    CalculateVelocityDeformationRate( m_velocityDeformationRate, fields, ibData, mesh );

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
                    
                    // Cell face velocity magntiude
                    // const floatType faceVelocityMagnitude = sqrt(
                    //                                               std::pow( fields.U[X](loCellIndexG) * ( 1.0f - lambda )  +  fields.U[X](hiCellIndexG) * lambda , 2.0f )
                    //                                             + std::pow( fields.U[Y](loCellIndexG) * ( 1.0f - lambda )  +  fields.U[Y](hiCellIndexG) * lambda , 2.0f )
                    //                                             + std::pow( fields.U[Z](loCellIndexG) * ( 1.0f - lambda )  +  fields.U[Z](hiCellIndexG) * lambda , 2.0f ) 
                    //                                         );
                    
                    const floatType hiVelocityMagnitude = sqrt( std::pow(fields.U[X](hiCellIndexG), 2.0f) + std::pow(fields.U[Y](hiCellIndexG), 2.0f) + std::pow(fields.U[Z](hiCellIndexG), 2.0f) );
                    const floatType loVelocityMagnitude = sqrt( std::pow(fields.U[X](loCellIndexG), 2.0f) + std::pow(fields.U[Y](loCellIndexG), 2.0f) + std::pow(fields.U[Z](loCellIndexG), 2.0f) );
                    const floatType faceVelocityMagnitude = ( 1.0f - lambda ) * loVelocityMagnitude  +  lambda * hiVelocityMagnitude;

                    // Cell face velocity deformation rate
                    const floatType Sface = (1.0f - lambda) * m_velocityDeformationRate(loCellIndexG)  +  lambda * m_velocityDeformationRate(hiCellIndexG);

                    // Vertical coordinate
                    const floatType z = ( m_heightAxis == faceNormal ) ? mesh.cellFaces[m_heightAxis]( faceIndex[m_heightAxis] )
                                                                       : mesh.cellCenters[m_heightAxis]( faceIndex[m_heightAxis] );

                    const floatType nuIn  = Sface 
                                          * std::pow( 
                                                  1.8f
                                                * ( 1.0f - exp( -0.645 * std::pow( m_averageBuildingWidth / m_averageBuildingHeight, 0.8 ) ) )
                                                * ( -2.0f * std::min( z / m_averageBuildingHeight , static_cast<floatType>( 1.0 ) ) )
                                                * std::pow( m_wallDistance[faceNormal](faceIndex), 2.0 )
                                            , 2.0f );

                    const floatType nuOut = sqrt( m_Cmu )
                                          * std::pow( m_Ig, 2.0f )
                                          * std::pow( m_referenceHeight, 2.0f * m_alpha + 0.1 )
                                          * std::pow( z, 0.9 - 2 * m_alpha ) / m_alpha
                                          * faceVelocityMagnitude
                                          * m_wallDistance[faceNormal](faceIndex);

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

                // Cell face velocity deformation rate. Just take internal cell value
                const floatType Sface = m_velocityDeformationRate( G(cellIndex) );

                // Vertical coordinate
                const floatType z = ( m_heightAxis == faceNormal ) ? mesh.cellFaces[m_heightAxis]( faceIndex[m_heightAxis] )
                                                                   : mesh.cellCenters[m_heightAxis]( faceIndex[m_heightAxis] );

                const floatType nuIn  = Sface 
                                        * std::pow( 
                                                1.8f
                                            * ( 1.0f - exp( -0.645 * std::pow( m_averageBuildingWidth / m_averageBuildingHeight, 0.8 ) ) )
                                            * ( -2.0f * std::min( z / m_averageBuildingHeight , static_cast<floatType>( 1.0 ) ) )
                                            * std::pow( m_wallDistance[faceNormal](faceIndex), 2.0 )
                                        , 2.0f );

                const floatType nuOut = sqrt( m_Cmu )
                                        * std::pow( m_Ig, 2.0f )
                                        * std::pow( m_referenceHeight, 2.0f * m_alpha + 0.1 )
                                        * std::pow( z, 0.9 - 2 * m_alpha ) / m_alpha
                                        * faceVelocityMagnitude
                                        * m_wallDistance[faceNormal](faceIndex);

                // Modify the turbulent viscosity
                nuTurbulent[faceNormal](faceIndex) = std::max( nuIn, nuOut );

            }
        }
    }


}

}   // end namespace CAMIRA