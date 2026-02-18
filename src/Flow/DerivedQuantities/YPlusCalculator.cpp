#include "DerivedQuantities.h"
#include "../../Core/FVTools.h"
#include "../../Core/FVLookups.h"
#include "../../Core/Geometry/Geometry.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"

#include <algorithm>
#include <numeric>
#include <iostream>

namespace CAMIRA
{


namespace 
{


bool IsInternalDomainCell( const TensorIndex3D &cellIndex,
                           const Mesh &mesh )
{
    bool isInternalDomainCell = true;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        if ( cellIndex[axis] < 0 || cellIndex[axis] >= mesh.nCells[axis] )
            isInternalDomainCell = false;
    } );
    return isInternalDomainCell;
}


bool HasDomainWallNeighbourCell( const TensorIndex3D &cellIndex,
                                 const BoundaryConditionData &bcData,
                                 const Mesh &mesh )
{
    bool hasDomainWallNeighbourCell = false;

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        
        // Negative side
        if ( cellIndex[axis] == 0 ) {
            if ( bcData.isWall[ LUT::NegativePatch[axis] ] )
                hasDomainWallNeighbourCell = true;
        }

        // Positive side
        if ( cellIndex[axis] == mesh.nCells[axis]-1 ) {
            if ( bcData.isWall[ LUT::PositivePatch[axis] ] )
                hasDomainWallNeighbourCell = true;
        }

    } );

    return hasDomainWallNeighbourCell;

}


bool HasNeighbouringSolidCell( intType i,
                               intType j,
                               intType k,
                               const Tensor3D &mask,
                               const Mesh &mesh,
                               const BoundaryConditionData &bcData )
{
    using FVT::G;

    TensorIndex3D cellN = {i  , j+1, k  },
                  cellS = {i  , j-1, k  },
                  cellE = {i+1, j  , k  },
                  cellW = {i-1, j  , k  },
                  cellT = {i  , j  , k+1},
                  cellB = {i  , j  , k-1};

    // Need to be careful because ghost cells are masked at domain boundaries, but they are not neccesarily solid cells.
    // First check for IB cells
    bool hasNeighbourIBCell =  ( static_cast<intType>( mask( G(cellN) ) ) == CellType::Solid  &&  IsInternalDomainCell(cellN, mesh) )
                            || ( static_cast<intType>( mask( G(cellS) ) ) == CellType::Solid  &&  IsInternalDomainCell(cellS, mesh) ) 
                            || ( static_cast<intType>( mask( G(cellE) ) ) == CellType::Solid  &&  IsInternalDomainCell(cellE, mesh) )
                            || ( static_cast<intType>( mask( G(cellW) ) ) == CellType::Solid  &&  IsInternalDomainCell(cellW, mesh) )
                            || ( static_cast<intType>( mask( G(cellT) ) ) == CellType::Solid  &&  IsInternalDomainCell(cellT, mesh) )
                            || ( static_cast<intType>( mask( G(cellB) ) ) == CellType::Solid  &&  IsInternalDomainCell(cellB, mesh) );

    if ( hasNeighbourIBCell ) {
        return true;
    }
    
    // If there is no IB cells, it is possible there is a domain boundary adjacent to the current cell
    if ( HasDomainWallNeighbourCell( {i, j, k}, bcData, mesh ) ) {
        return true;
    }
    
    return false;

}


}   // end anonymous namespace



YPlusCalculator::YPlusCalculator( const InputData &inputData,
                                  const AxisTransformationMap &axisTransformation,
                                  const BoundaryConditionData &bcData,
                                  const IBData &ibData, 
                                  const Mesh &mesh, 
                                  const FieldData<Tensor3D> &fields) :
                                    m_fields( fields ), 
                                    m_rho( inputData.rho ), 
                                    m_nu( inputData.nu ) 
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    Polyhedron P = MakeGeometry( inputData.geometryData, axisTransformation );
    Tree tree = MakeAABBTree( P );

    floatType inf = std::numeric_limits<floatType>::infinity();

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                TensorIndex3D cellIndex = {i, j, k};

                // Ignore if it is inside the geometry
                if ( ibData.mask( G(cellIndex) ) == static_cast<floatType>( CellType::Solid ) )
                    continue;

                // Check if it is the nearest cell to the wall
                if ( !HasNeighbouringSolidCell(i, j, k, ibData.mask, mesh, bcData) )
                    continue;

                const floatType xq = mesh.cellCenters[X](i),
                                yq = mesh.cellCenters[Y](j),
                                zq = mesh.cellCenters[Z](k);
                const fVector3 cellCoords{xq, yq, zq };

                fVector3 nearestWallPoint = cellCoords;
                fVector3 wallTangentialVelocity = {0.0f, 0.0f, 0.0f};
                floatType nearestWallPointDistance = inf;
                
                // Distance to nearest geometry object
                if ( !tree.empty() ) {
                    nearestWallPoint = NearestPoint( tree, xq, yq, zq );
                    nearestWallPointDistance = sqrt( std::pow( xq - nearestWallPoint(X), 2.0f ) 
                                                   + std::pow( yq - nearestWallPoint(Y), 2.0f )
                                                   + std::pow( zq - nearestWallPoint(Z), 2.0f ) );

                    wallTangentialVelocity = {0.0f, 0.0f, 0.0f};
                }


                // Check if there is a closer point on the domain boundaries
                EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

                    if ( bcData.isWall[ LUT::NegativePatch[axis] ] ) {
                        floatType domainWallDistance = std::min( nearestWallPointDistance,
                                                                 abs( cellCoords[axis] - mesh.cellFaces[axis]( 0 ) ) );

                        if ( domainWallDistance < nearestWallPointDistance ) {
                            nearestWallPointDistance = domainWallDistance;

                            EnumFor<Axis>( [&] (Axis::ENUMDATA component) {
                                nearestWallPoint(component) = ( axis == component ) ? mesh.cellFaces[axis](0) 
                                                                                    : mesh.cellCenters[axis](0);
                            } );

                            floatType wallTangentLoVelocityComponent = bcData.fields.U[ LUT::LoOrthogonalAxis[axis] ][ LUT::NegativePatch[axis] ].value( cellIndex[ LUT::LoOrthogonalAxis[axis] ],
                                                                                                                                                         cellIndex[ LUT::HiOrthogonalAxis[axis] ] );
                            
                            floatType wallTangentHiVelocityComponent = bcData.fields.U[ LUT::HiOrthogonalAxis[axis] ][ LUT::NegativePatch[axis] ].value( cellIndex[ LUT::LoOrthogonalAxis[axis] ],
                                                                                                                                                         cellIndex[ LUT::HiOrthogonalAxis[axis] ] );

                            wallTangentialVelocity(axis) = 0.0f;
                            wallTangentialVelocity( LUT::LoOrthogonalAxis[axis] ) = wallTangentLoVelocityComponent;
                            wallTangentialVelocity( LUT::HiOrthogonalAxis[axis] ) = wallTangentHiVelocityComponent;

                        }
                    }


                    if ( bcData.isWall[ LUT::PositivePatch[axis] ] ) {
                        floatType domainWallDistance = std::min( nearestWallPointDistance,
                                                                 abs( mesh.cellFaces[axis]( mesh.nFacesNormal[axis](axis) - 1 ) - cellCoords[axis] ) );

                        if ( domainWallDistance < nearestWallPointDistance ) {
                            nearestWallPointDistance = domainWallDistance;
                            
                            EnumFor<Axis>( [&] (Axis::ENUMDATA component) {
                                nearestWallPoint(component) = ( axis == component ) ? mesh.cellFaces[axis]( mesh.nFacesNormal[axis](axis) - 1 ) 
                                                                                    : mesh.cellCenters[axis]( mesh.nCells[axis] - 1 );
                            } );

                            floatType wallTangentLoVelocityComponent = bcData.fields.U[ LUT::LoOrthogonalAxis[axis] ][ LUT::PositivePatch[axis] ].value( cellIndex[ LUT::LoOrthogonalAxis[axis] ],
                                                                                                                                                         cellIndex[ LUT::HiOrthogonalAxis[axis] ] );
                            
                            floatType wallTangentHiVelocityComponent = bcData.fields.U[ LUT::HiOrthogonalAxis[axis] ][ LUT::PositivePatch[axis] ].value( cellIndex[ LUT::LoOrthogonalAxis[axis] ],
                                                                                                                                                         cellIndex[ LUT::HiOrthogonalAxis[axis] ] );

                            wallTangentialVelocity(axis) = 0.0f;
                            wallTangentialVelocity( LUT::LoOrthogonalAxis[axis] ) = wallTangentLoVelocityComponent;
                            wallTangentialVelocity( LUT::HiOrthogonalAxis[axis] ) = wallTangentHiVelocityComponent;

                        }
                    }

                } );

                // Outward pointing (unit) normal vector from solid surface
                fVector3 normalVector = cellCoords - nearestWallPoint;
                normalVector /= normalVector.norm(); 


                m_wallCells.emplace_back( 0.0f,
                                          cellIndex,
                                          nearestWallPointDistance,
                                          normalVector, 
                                          wallTangentialVelocity );

            }
        }
    }

    // Update the y+ values
    Update();

}



void YPlusCalculator::Update()
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    // Possible there is not wall cells
    if ( m_wallCells.empty() ) {
        minYPlus = 0.0f;
        maxYPlus = 0.0f;
        averageYPlus = 0.0f;
        return;
    }

    for ( auto &wallCellData : m_wallCells ) {

        // Approximate the wall shear stress with one sided difference
        fVector3 fluidVelocity = { m_fields.U[X]( G( wallCellData.cellIndex ) ),
                                   m_fields.U[Y]( G( wallCellData.cellIndex ) ),
                                   m_fields.U[Z]( G( wallCellData.cellIndex ) ) };

        fVector3 relativeVelocityToWall = fluidVelocity - wallCellData.wallTangentialVelocity;

        fVector3 tangentialVelocityComponent = relativeVelocityToWall - relativeVelocityToWall.dot( wallCellData.normalVector ) * wallCellData.normalVector;

        floatType wallShearStress = m_nu * tangentialVelocityComponent.norm() / wallCellData.wallDistance;

        floatType frictionVelocity = sqrt( wallShearStress );

        wallCellData.yPlus = frictionVelocity * wallCellData.wallDistance / m_nu;

    }

    // Calculate the statistics
    auto compare = []( const WallCellData &a, const WallCellData &b ) {
        return a.yPlus < b.yPlus;
    };

    minYPlus = std::min_element( m_wallCells.begin(), m_wallCells.end(), compare )->yPlus;
    maxYPlus = std::max_element( m_wallCells.begin(), m_wallCells.end(), compare )->yPlus;

    auto add = [](floatType sum, const WallCellData &cell) {
        return sum + cell.yPlus;
    };
    averageYPlus = std::accumulate( m_wallCells.begin(), m_wallCells.end(), 0.0f, add ) / static_cast<floatType>( m_wallCells.size() );
}


}   // end namespace CAMIRA