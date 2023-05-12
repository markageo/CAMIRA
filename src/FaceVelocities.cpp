#include "FiniteVolume.h"

namespace CFD
{

namespace
{

void FaceVelocityXnormal( array3D &faceVel, 
                          const array3D &cellVel, 
                          const Mesh &mesh)
{
    using enum Axis::ENUMDATA;

    floatType interpFactor;
    for (intType k = 0; k != faceVel.dimension(2); k++ ) {
        for (intType j = 0; j != faceVel.dimension(1); j++) {
            for (intType i = 1; i != faceVel.dimension(0)-1; i++) {

                interpFactor = mesh.interpFactors[X](i);
                faceVel(i, j, k) = (1 - interpFactor)*cellVel( G(i-1, j, k) ) + interpFactor*cellVel( G(i, j, k) );

            }
        }
    }
}


void FaceVelocityYnormal( array3D &faceVel, 
                          const array3D &cellVel, 
                          const Mesh &mesh)
{
    using enum Axis::ENUMDATA;

    floatType interpFactor;
    for (intType k = 0; k != faceVel.dimension(2); k++ ) {
        for (intType j = 1; j != faceVel.dimension(1)-1; j++) {
            for (intType i = 0; i != faceVel.dimension(0); i++) {
                
                interpFactor = mesh.interpFactors[Y](j);
                faceVel(i, j, k) = (1 - interpFactor)*cellVel( G(i, j-1, k) ) + interpFactor*cellVel( G(i, j, k) );

            }
        }
    }
}


void FaceVelocityZnormal( array3D &faceVel, 
                          const array3D &cellVel, 
                          const Mesh &mesh)
{
    using enum Axis::ENUMDATA;

    floatType interpFactor;
    for (intType k = 1; k != faceVel.dimension(2)-1; k++ ) {
        for (intType j = 0; j != faceVel.dimension(1); j++) {
            for (intType i = 0; i != faceVel.dimension(0); i++) {

                interpFactor = mesh.interpFactors[Z](k);
                faceVel(i, j, k) = (1 - interpFactor)*cellVel( G(i, j, k-1) ) + interpFactor*cellVel( G(i, j, k) );

            }
        }
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


void UpdateFaceVelocities( ArrayAllocator<Fields, CFD::array3D> &faceVelocities, 
                           const Mesh &mesh, 
                           const ArrayAllocator<Fields, CFD::array3D> &fields, 
                           const InputData &inputData)
{
    using F = Fields::ENUMDATA;
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;

    const InputData::BoundaryConditionData &boundaryConditions = inputData.boundaryConditions;

    // Non-boundary cells
    FaceVelocityXnormal( faceVelocities[F::U], fields[F::U], mesh);
    FaceVelocityYnormal( faceVelocities[F::V], fields[F::V], mesh);
    FaceVelocityZnormal( faceVelocities[F::W], fields[F::W], mesh);

    // To allow each axis to be computed by looping
    constexpr std::array<Fields::ENUMDATA, 3> faceVelocityFields = {F::U, F::V, F::W};
    Fields::ENUMDATA axisVel;
    intType faceEndIndex, fieldEndIndex;
    floatType extrapFactor_p, extrapFactor_a;
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;

    // Ghost cells mean that a slice of the fields tensor must be used when setting BCs
    Eigen::array<intType, 3> offsets = {nGhost, nGhost, nGhost};
    Eigen::array<intType, 3> extents = {mesh.nCells(X), mesh.nCells(Y), mesh.nCells(Z)};

    for (int axis = 0; axis != Axis::count; axis++) {

        positivePatch = PositivePatch[axis];
        negativePatch = NegativePatch[axis];
        axisVel = faceVelocityFields[axis];

        // Axis positive boundary
        faceEndIndex = mesh.nCells(axis);
        fieldEndIndex = mesh.nCells(axis)-1;
        switch ( boundaryConditions[axisVel][positivePatch].type ) {
            
            case BC::zeroGradient:
                faceVelocities[axisVel].chip(faceEndIndex, axis) = fields[axisVel].slice(offsets, extents).chip(fieldEndIndex , axis);          
                break;

            case BC::uniform:
                faceVelocities[axisVel].chip(faceEndIndex, axis) = faceVelocities[axisVel].chip(faceEndIndex, axis).constant( boundaryConditions[axisVel][positivePatch].value );
                break;

            case BC::extrapolated:
                extrapFactor_p = mesh.extrapFactors[positivePatch].p;
                extrapFactor_a = mesh.extrapFactors[positivePatch].a;
                faceVelocities[axisVel].chip(faceEndIndex, axis) = fields[axisVel].slice(offsets, extents).chip(fieldEndIndex  , axis) 
                                                                 * fields[axisVel].slice(offsets, extents).chip(fieldEndIndex  , axis).constant( extrapFactor_p )
                                                                 + fields[axisVel].slice(offsets, extents).chip(fieldEndIndex-1, axis) 
                                                                 * fields[axisVel].slice(offsets, extents).chip(fieldEndIndex-1, axis).constant( extrapFactor_a );
                break;

            default:
                break;
        }


        // Axis negative boundary
        faceEndIndex = 0;
        fieldEndIndex = 0;
        switch ( boundaryConditions[axisVel][negativePatch].type ) {
            
            case BC::zeroGradient:
                faceVelocities[axisVel].chip(faceEndIndex, axis) = fields[axisVel].slice(offsets, extents).chip(fieldEndIndex, axis);          
                break;

            case BC::uniform:
                faceVelocities[axisVel].chip(faceEndIndex, axis) = faceVelocities[axisVel].chip(faceEndIndex, axis).constant( boundaryConditions[axisVel][negativePatch].value );
                break;

            case BC::extrapolated:
                extrapFactor_p = mesh.extrapFactors[negativePatch].p;
                extrapFactor_a = mesh.extrapFactors[negativePatch].a;
                faceVelocities[axisVel].chip(faceEndIndex, axis) = fields[axisVel].slice(offsets, extents).chip(fieldEndIndex  , axis) 
                                                                                * fields[axisVel].slice(offsets, extents).chip(fieldEndIndex, axis).constant( extrapFactor_p )
                                                                 + fields[axisVel].slice(offsets, extents).chip(fieldEndIndex+1, axis) 
                                                                                * fields[axisVel].slice(offsets, extents).chip(fieldEndIndex+1, axis).constant( extrapFactor_a );
                break;

            default:
                break;
        }
    }

}

} // end namespace CFD