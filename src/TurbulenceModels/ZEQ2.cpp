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



void TurbulenceModel<TurbulenceModels::ZEQ2>::SetTurbulenceViscosityField( EnumVector<Axis, Tensor3D> &nuTurbulent,
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

                    const floatType aReTurb = std::pow( ( m_inflowTKEBuildingHeight * m_inflowIntergralTimeScaleBuildingHeight ) / m_nu, 1.0f/3.0f );
                    const floatType bReBulk = std::pow( ( m_inflowVelocityBuildingHeight * m_averageBuildingHeight ) / m_nu            , 1.0f/3.0f );

                    const floatType wallDistanceNormalised = m_wallDistance[faceNormal](faceIndex) / m_wallDistanceLengthScale;
                    nuTurbulent[faceNormal](faceIndex) = aReTurb 
                                                       * wallDistanceNormalised
                                                       * exp( -bReBulk * wallDistanceNormalised )
                                                       * faceVelocityMagnitude
                                                       * m_wallDistance[faceNormal](faceIndex);
                    
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
                const floatType aReTurb = std::pow( ( m_inflowTKEBuildingHeight * m_inflowIntergralTimeScaleBuildingHeight ) / m_nu, 1.0f/3.0f );
                const floatType bReBulk = std::pow( ( m_inflowVelocityBuildingHeight * m_averageBuildingHeight ) / m_nu            , 1.0f/3.0f );

                const floatType wallDistanceNormalised = m_wallDistance[faceNormal](faceIndex) / m_wallDistanceLengthScale;
                nuTurbulent[faceNormal](faceIndex) = aReTurb 
                                                    * wallDistanceNormalised
                                                    * exp( -bReBulk * wallDistanceNormalised )
                                                    * faceVelocityMagnitude
                                                    * m_wallDistance[faceNormal](faceIndex);

            }
        }
    }


}

}   // end namespace CAMIRA