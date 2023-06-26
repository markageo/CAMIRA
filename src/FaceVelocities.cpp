#include "FiniteVolume.h"

namespace CFD
{

namespace
{


void FaceVelocity( EnumVector<Axis, array3D> &faceFluxes, 
                   const EnumVector<Axis, array3D> &cellVelocities, 
                   const Mesh &mesh, 
                   const Axis::ENUMDATA axis,
                   const Axis::ENUMDATA velocityComponent )
{
    using enum Axis::ENUMDATA;

    array3D &faceVel = faceFluxes[ velocityComponent ];
    const array3D &cellVel = cellVelocities[ velocityComponent ];

    // Starting index and number of faces to iterate over
    iVector3 startIndex, nFaces;
    EnumFor<Axis>( [&] ( Axis::ENUMDATA a) {
        startIndex[a] = 0;
        nFaces[a] = faceFluxes[ velocityComponent ].dimension(a);
    } );
    startIndex[axis] += 1;
    nFaces[axis] -= 1;

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



void BoundaryFaceVelocitiy( EnumVector<Axis, array3D> &faceFluxes, 
                            const EnumVector<Axis, array3D> &cellVelocities, 
                            const Mesh &mesh, 
                            const EnumVector< Axis, EnumVector< BoundaryPatches, InputData::BoundaryConditionData > >& boundaryConditions,
                            const BoundaryPatches::ENUMDATA boundaryPatch,
                            const Axis::ENUMDATA component )
{
    using BC = BoundaryConditions::ENUMDATA;
    
    Axis::ENUMDATA axis = BoundaryPatchAxis[ boundaryPatch ];
    
    static constexpr arrayIndex3D offsets = {nGhost, nGhost, nGhost};
    arrayIndex3D extents = {mesh.nCells(Axis::X), mesh.nCells(Axis::Y), mesh.nCells(Axis::Z)};
    
    intType faceEndIndex, fieldEndIndex;
    if ( boundaryPatch == PositivePatch[ axis ] ) {
        faceEndIndex = mesh.nCells(axis);
        fieldEndIndex = mesh.nCells(axis)-1;
    } else {
        faceEndIndex = 0;
        fieldEndIndex = 0;
    }

    floatType extrapFactor_p, extrapFactor_a;
    switch ( boundaryConditions[component][boundaryPatch].type ) 
    {    
        case BC::zeroGradient:
            faceFluxes[component].chip(faceEndIndex, axis) = cellVelocities[component].slice(offsets, extents).chip(fieldEndIndex, axis);          
            break;

        case BC::uniform:
            faceFluxes[component].chip(faceEndIndex, axis) = faceFluxes[component].chip(faceEndIndex, axis).constant( boundaryConditions[component][boundaryPatch].value );
            break;

        case BC::extrapolated:
            extrapFactor_p = mesh.extrapFactors[boundaryPatch].p;
            extrapFactor_a = mesh.extrapFactors[boundaryPatch].a;
            faceFluxes[component].chip(faceEndIndex, axis) = cellVelocities[component].slice(offsets, extents).chip(fieldEndIndex  , axis) 
                                                            * cellVelocities[component].slice(offsets, extents).chip(fieldEndIndex  , axis).constant( extrapFactor_p )
                                                       + cellVelocities[component].slice(offsets, extents).chip(fieldEndIndex+1, axis) 
                                                            * cellVelocities[component].slice(offsets, extents).chip(fieldEndIndex+1, axis).constant( extrapFactor_a );
            break;

        default:
            break;
    }


}

}   // end anonymous namespace





EnumVector<Axis, array3D> InitialiseFaceFluxes( const Mesh &mesh, 
                                                const EnumVector<Axis, array3D> &cellVelocities, 
                                                const InputData &inputData)
{
    // Faces are staggered in the negative direction:
    //   cellFaceVelocity_x(i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceVelocity_y(i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceVelocity_z(i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    EnumVector<Axis, array3D> faceFluxes( {{Axis::X, {mesh.nCells(0) + 1, mesh.nCells(1)    , mesh.nCells(2)    }},
                                           {Axis::Y, {mesh.nCells(0)    , mesh.nCells(1) + 1, mesh.nCells(2)    }},
                                           {Axis::Z, {mesh.nCells(0)    , mesh.nCells(1)    , mesh.nCells(2) + 1}}} );
                                                     
    UpdateFaceFluxes(faceFluxes, mesh, cellVelocities, inputData);

    return faceFluxes;
}




// Calculates face velocity fluxes. i.e. normal component of velocity on faces
void UpdateFaceFluxes( EnumVector<Axis, array3D> &faceFluxes, 
                       const Mesh &mesh, 
                       const EnumVector<Axis, array3D> &cellVelocities, 
                       const InputData &inputData )
{
    // Internal faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        FaceVelocity( faceFluxes, cellVelocities, mesh, axis, axis);
    } );
    
    // Boundary faces
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA boundaryPatch) {

        Axis::ENUMDATA axis = BoundaryPatchAxis[ boundaryPatch ];
        BoundaryFaceVelocitiy( faceFluxes, cellVelocities, mesh, inputData.boundaryConditions.U, boundaryPatch, axis );

    } );
}

} // end namespace CFD