#include "FiniteVolume.h"

namespace CFD
{

namespace
{


void FaceVelocity( ArrayAllocator<Fields, array3D> &faceVelocities, 
                   const ArrayAllocator<Fields, array3D> &fields, 
                   const Mesh &mesh, 
                   const Axis::ENUMDATA axis,
                   const Fields::ENUMDATA field )
{
    using enum Axis::ENUMDATA;

    array3D &faceVel = faceVelocities[ field ];
    const array3D &cellVel = fields[ field ];

    // Starting index and number of faces to iterate over
    iVector3 startIndex, nFaces;
    EnumFor<Axis>( [&] ( Axis::ENUMDATA a) {
        startIndex[a] = 0;
        nFaces[a] = faceVelocities[ field ].dimension(a);
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


void BoundaryFaceVelocitiy( ArrayAllocator<Fields, array3D> &faceVelocities, 
                            const ArrayAllocator<Fields, array3D> &fields, 
                            const Mesh &mesh, 
                            const InputData::BoundaryConditionData& boundaryConditions,
                            const BoundaryPatches::ENUMDATA boundaryPatch,
                            const Fields::ENUMDATA field )
{
    using BC = BoundaryConditions::ENUMDATA;
    
    Axis::ENUMDATA axis = BoundaryPatchAxis[ boundaryPatch ];
    Fields::ENUMDATA axisVel = AxisVelocity[ axis ];
    
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
    switch ( boundaryConditions[field][boundaryPatch].type ) 
    {    
        case BC::zeroGradient:
            faceVelocities[field].chip(faceEndIndex, axis) = fields[field].slice(offsets, extents).chip(fieldEndIndex, axis);          
            break;

        case BC::uniform:
            faceVelocities[field].chip(faceEndIndex, axis) = faceVelocities[field].chip(faceEndIndex, axis).constant( boundaryConditions[field][boundaryPatch].value );
            break;

        case BC::extrapolated:
            extrapFactor_p = mesh.extrapFactors[boundaryPatch].p;
            extrapFactor_a = mesh.extrapFactors[boundaryPatch].a;
            faceVelocities[field].chip(faceEndIndex, axis) = fields[field].slice(offsets, extents).chip(fieldEndIndex  , axis) 
                                                                * fields[field].slice(offsets, extents).chip(fieldEndIndex  , axis).constant( extrapFactor_p )
                                                             + fields[field].slice(offsets, extents).chip(fieldEndIndex+1, axis) 
                                                                * fields[field].slice(offsets, extents).chip(fieldEndIndex+1, axis).constant( extrapFactor_a );
            break;

        default:
            break;
    }


}

}   // end anonymous namespace





ArrayAllocator<Fields, array3D> InitialiseFaceVelocities(const Mesh &mesh, 
                                                         const ArrayAllocator<Fields, array3D> &fields, 
                                                         const InputData &inputData)
{
    using F = Fields::ENUMDATA;

    // Faces are staggered in the negative direction:
    //   cellFaceVelocity_x(i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceVelocity_y(i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceVelocity_z(i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    ArrayAllocator<Fields, array3D> faceVelocities( {{F::U, {mesh.nCells(0) + 1, mesh.nCells(1)    , mesh.nCells(2)    }},
                                                     {F::V, {mesh.nCells(0)    , mesh.nCells(1) + 1, mesh.nCells(2)    }},
                                                     {F::W, {mesh.nCells(0)    , mesh.nCells(1)    , mesh.nCells(2) + 1}}} );
                                                     
    UpdateFaceVelocities(faceVelocities, mesh, fields, inputData);

    return faceVelocities;
}




// Calculates face velocity fluxes. i.e. normal component of velocity on faces
void UpdateFaceVelocities( ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                           const Mesh &mesh, 
                           const ArrayAllocator<Fields, CFD::array3D> &fields, 
                           const InputData &inputData)
{
    // Internal faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        FaceVelocity( faceVelocities, fields, mesh, axis, AxisVelocity[axis]);
    } );
    
    // Boundary faces
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA boundaryPatch) {
        Fields::ENUMDATA field = AxisVelocity[ BoundaryPatchAxis[ boundaryPatch ] ];
        BoundaryFaceVelocitiy( faceVelocities, fields, mesh, inputData.boundaryConditions, boundaryPatch, field );
    } );
}

} // end namespace CFD