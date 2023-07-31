#include "FiniteVolume.h"
#include "../Tools/FVTools.h"
#include "../Tools/FVLookups.h"

namespace CFD
{

using namespace FVT;

namespace
{

void LinearInterpolationInteriorFaceVelocity( EnumVector<Axis, array3D> &faceVelocities, 
                                              const EnumVector<Axis, array3D> &cellVelocities, 
                                              const Mesh &mesh, 
                                              const Axis::ENUMDATA axis,
                                              const Axis::ENUMDATA velocityComponent )
{
    using enum Axis::ENUMDATA;

    array3D &faceVel = faceVelocities[ velocityComponent ];
    const array3D &cellVel = cellVelocities[ velocityComponent ];

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++ ) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D idx = {i, j, k},
                             HiIndex = idx,
                             LoIndex = idx;
                LoIndex[axis] -= 1;
                floatType interpFactor = mesh.interpFactors[ axis ]( idx[axis] );
                faceVel( idx ) = (1 - interpFactor)*cellVel( G(LoIndex) ) + interpFactor*cellVel( G(HiIndex) );

            }
        }
    }
}




void UpwindInteriorFaceVelocity( EnumVector<Axis, array3D> &faceVelocities, 
                                 const EnumVector<Axis, array3D> &cellVelocities, 
                                 const EnumVector<Axis, array3D> &faceFluxes,
                                 const Mesh &mesh, 
                                 const Axis::ENUMDATA axis,
                                 const Axis::ENUMDATA velocityComponent )
{
    using enum Axis::ENUMDATA;

    array3D &faceVel = faceVelocities[ velocityComponent ];
    const array3D &cellVel = cellVelocities[ velocityComponent ];

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++ ) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                arrayIndex3D idx = {i, j, k},
                             HiIndex = idx,
                             LoIndex = idx;
                LoIndex[axis] -= 1;

                if ( faceFluxes[axis](idx) >= 0.0f ) {
                    faceVel( idx ) = cellVel( G(LoIndex) );
                } else {
                    faceVel( idx ) = cellVel( G(HiIndex) );
                }

            }
        }
    }
}





void BoundaryFaceVelocitiy( EnumVector<Axis, array3D> &faceVelocities, 
                            const EnumVector<Axis, array3D> &cellVelocities, 
                            const Mesh &mesh, 
                            const EnumVector< Axis, BoundaryConditionData >& boundaryConditions,
                            const BoundaryPatches::ENUMDATA boundaryPatch,
                            const Axis::ENUMDATA component )
{
    using BC = BoundaryConditions::ENUMDATA;
    
    Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[ boundaryPatch ];
    
    static constexpr arrayIndex3D offsets = {nGhost, nGhost, nGhost};
    arrayIndex3D extents = {mesh.nCells(Axis::X), mesh.nCells(Axis::Y), mesh.nCells(Axis::Z)};
    
    intType faceEndIndex, fieldEndIndex;
    if ( boundaryPatch == LUT::PositivePatch[ axis ] ) {
        faceEndIndex = mesh.nCells(axis);
        fieldEndIndex = mesh.nCells(axis)-1;
    } else {
        faceEndIndex = 0;
        fieldEndIndex = 0;
    }

    switch ( boundaryConditions[component][boundaryPatch].type ) 
    {    
        case BC::zeroGradient:
        {
            faceVelocities[component].chip(faceEndIndex, axis) = cellVelocities[component].slice(offsets, extents).chip(fieldEndIndex, axis);          
            break;
        }
            

        case BC::fixed:
        {
            faceVelocities[component].chip(faceEndIndex, axis) = boundaryConditions[component][boundaryPatch].value;
            break;
        }
            

        case BC::extrapolated:
        {
            floatType extrapFactor_p = mesh.extrapFactors[boundaryPatch].p;
            floatType extrapFactor_a = mesh.extrapFactors[boundaryPatch].a;
            faceVelocities[component].chip(faceEndIndex, axis) = cellVelocities[component].slice(offsets, extents).chip(fieldEndIndex  , axis) 
                                                                    * cellVelocities[component].slice(offsets, extents).chip(fieldEndIndex  , axis).constant( extrapFactor_p )
                                                               + cellVelocities[component].slice(offsets, extents).chip(fieldEndIndex+1, axis) 
                                                                    * cellVelocities[component].slice(offsets, extents).chip(fieldEndIndex+1, axis).constant( extrapFactor_a );
            break;
        }
            

        default:
            break;
    }


}

}   // end anonymous namespace


// ---------------------------------------- Face Advected Velocities ----------------------------------------

// Calculates face velocity fluxes. i.e. normal component of velocity on faces
void UpdateFaceAdvectedVelocities( EnumVector< Axis, EnumVector< Axis, array3D > > &faceAdvectedVelocities, 
                                   const Mesh &mesh, 
                                   const EnumVector< Axis, array3D > &cellVelocities, 
                                   const EnumVector< Axis, array3D > &faceFluxes,
                                   const FieldData< BoundaryConditionData > &bcData )
{
    // Internal faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        EnumFor<Axis>( [&] (Axis::ENUMDATA velocityComponent) {
            UpwindInteriorFaceVelocity( faceAdvectedVelocities[axis], cellVelocities, faceFluxes, mesh, axis, velocityComponent);
        } );
    } );
    
    // Boundary faces
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA boundaryPatch) {

        Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[ boundaryPatch ];
        EnumFor<Axis>( [&] (Axis::ENUMDATA velocityComponent) {
            BoundaryFaceVelocitiy( faceAdvectedVelocities[axis], cellVelocities, mesh, bcData.U, boundaryPatch, velocityComponent );
        } );

    } );
}



EnumVector< Axis, EnumVector<Axis, array3D> > InitialiseAdvectedFaceVelocities( const Mesh &mesh, 
                                                                                const EnumVector<Axis, array3D> &cellVelocities, 
                                                                                const EnumVector<Axis, array3D> &faceFluxes,
                                                                                const FieldData< BoundaryConditionData > &bcData)
{
    // First index is the face normal. Second index is the velocity components.
    // Faces are staggered in the negative direction:
    //   cellFaceVelocity[X][X](i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceVelocity[Y][X](i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceVelocity[Z][X](i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    EnumVector< Axis, EnumVector<Axis, array3D> > faceAdvectedVelocities( { EnumVector<Axis, array3D>( array3D(mesh.nCells(0) + 1, mesh.nCells(1)    , mesh.nCells(2)    ).setZero() ),
                                                                    EnumVector<Axis, array3D>( array3D(mesh.nCells(0)    , mesh.nCells(1) + 1, mesh.nCells(2)    ).setZero() ),
                                                                    EnumVector<Axis, array3D>( array3D(mesh.nCells(0)    , mesh.nCells(1)    , mesh.nCells(2) + 1).setZero() ) } );
                                                     
    UpdateFaceAdvectedVelocities(faceAdvectedVelocities, mesh, cellVelocities, faceFluxes, bcData);

    return faceAdvectedVelocities;
}





// ---------------------------------------- Face Fluxes ----------------------------------------

// Calculates face velocity fluxes. i.e. normal component of velocity on faces
void UpdateFaceFluxes( EnumVector< Axis, array3D > &faceFluxes, 
                       const Mesh &mesh, 
                       const EnumVector< Axis, array3D > &cellVelocities, 
                       const FieldData< BoundaryConditionData > &bcData )
{
    // Internal faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        LinearInterpolationInteriorFaceVelocity( faceFluxes, cellVelocities, mesh, axis, axis);
    } );
    
    // Boundary faces
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA boundaryPatch) {

        Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[ boundaryPatch ];
        BoundaryFaceVelocitiy( faceFluxes, cellVelocities, mesh, bcData.U, boundaryPatch, axis );

    } );
}



EnumVector<Axis, array3D> InitialiseFaceFluxes( const Mesh &mesh, 
                                                const EnumVector<Axis, array3D> &cellVelocities, 
                                                const FieldData< BoundaryConditionData > &bcData)
{
    // Faces are staggered in the negative direction:
    //   cellFaceFlux[X](i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceFlux[Y](i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceFlux[Z](i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    EnumVector<Axis, array3D> faceFluxes( {{Axis::X, {mesh.nCells(0) + 1, mesh.nCells(1)    , mesh.nCells(2)    }},
                                           {Axis::Y, {mesh.nCells(0)    , mesh.nCells(1) + 1, mesh.nCells(2)    }},
                                           {Axis::Z, {mesh.nCells(0)    , mesh.nCells(1)    , mesh.nCells(2) + 1}}} );
                                                     
    UpdateFaceFluxes(faceFluxes, mesh, cellVelocities, bcData);

    return faceFluxes;
}


} // end namespace CFD