#include "FiniteVolumeFunctions.h"

// Implementation file for face velocity update functions

namespace
{

void FaceVelocityXnormal(CFD::array3D &faceVel, const CFD::array3D &cellVel, const CFD::Mesh &mesh)
{
    using namespace CFD;
    using enum Axis::ENUMDATA;

    floatType interpFactor;
    for (iterType k = 0; k != faceVel.dimension(2); k++ ) {
        for (iterType j = 0; j != faceVel.dimension(1); j++) {
            for (iterType i = 1; i != faceVel.dimension(0)-1; i++) {
                interpFactor = mesh.interpFactors[X](i);
                
                faceVel(i, j, k) = (1 - interpFactor)*cellVel(i-1, j, k) + interpFactor*cellVel(i, j, k);

            }
        }
    }
}


void FaceVelocityYnormal(CFD::array3D &faceVel, const CFD::array3D &cellVel, const CFD::Mesh &mesh)
{
    using namespace CFD;
    using enum Axis::ENUMDATA;

    floatType interpFactor;
    for (iterType k = 0; k != faceVel.dimension(2); k++ ) {
        for (iterType j = 1; j != faceVel.dimension(1)-1; j++) {
            interpFactor = mesh.interpFactors[Y](j);
            for (iterType i = 0; i != faceVel.dimension(0); i++) {
                
                faceVel(i, j, k) = (1 - interpFactor)*cellVel(i, j-1, k) + interpFactor*cellVel(i, j, k);

            }
        }
    }
}


void FaceVelocityZnormal(CFD::array3D &faceVel, const CFD::array3D &cellVel, const CFD::Mesh &mesh)
{
    using namespace CFD;
    using enum Axis::ENUMDATA;

    floatType interpFactor;
    for (iterType k = 1; k != faceVel.dimension(2)-1; k++ ) {
        interpFactor = mesh.interpFactors[Z](k);
        for (iterType j = 0; j != faceVel.dimension(1); j++) {
            for (iterType i = 0; i != faceVel.dimension(0); i++) {

                faceVel(i, j, k) = (1 - interpFactor)*cellVel(i, j, k-1) + interpFactor*cellVel(i, j, k);

            }
        }
    }
}


}   // end anonymous namespace


namespace CFD
{

void UpdateFaceVelocities( ArrayAllocator<Fields> &faceVelocities, const Mesh &mesh, const ArrayAllocator<Fields> &fields, 
    const InputData::BoundaryConditionData &boundaryConditions)
{
    using F = Fields::ENUMDATA;
    using BC = BoundaryConditions::ENUMDATA;
    using enum Axis::ENUMDATA;

    // Non-boundary cells
    FaceVelocityXnormal( faceVelocities[F::U], fields[F::U], mesh);
    FaceVelocityYnormal( faceVelocities[F::V], fields[F::V], mesh);
    FaceVelocityZnormal( faceVelocities[F::W], fields[F::W], mesh);

    // Just to assign the boundary index of the face flux and field arrays
    constexpr std::array<Fields::ENUMDATA, 3> faceVelocityFields = {F::U, F::V, F::W};
    Fields::ENUMDATA axisVel;
    iterType faceEndIndex, fieldEndIndex;
    floatType extrapFactor_p, extrapFactor_a;
    BoundaryPatches::ENUMDATA positivePatch, negativePatch;

    for (int axis = 0; axis != Axis::count; axis++) {

        positivePatch = positivePatches[axis];
        negativePatch = negativePatches[axis];
        axisVel = faceVelocityFields[axis];

        // Axis positive boundary
        faceEndIndex = faceVelocities[axisVel].dimension(axis)-1;
        fieldEndIndex = fields[axisVel].dimension(axis)-1;
        switch ( boundaryConditions[axisVel][positivePatch].type ) {
            
            case BC::zeroGradient:
                faceVelocities[axisVel].chip(faceEndIndex, axis) = fields[axisVel].chip(fieldEndIndex, axis);          
                break;

            case BC::uniform:
                faceVelocities[axisVel].chip(faceEndIndex, axis) = faceVelocities[axisVel].chip(faceEndIndex, axis).constant( boundaryConditions[axisVel][positivePatch].value );
                break;

            case BC::extrapolated:
                extrapFactor_p = mesh.extrapFactors[positivePatch].p;
                extrapFactor_a = mesh.extrapFactors[positivePatch].a;
                faceVelocities[axisVel].chip(faceEndIndex, axis) = fields[axisVel].chip(fieldEndIndex, axis) * fields[axisVel].chip(fieldEndIndex, axis).constant( extrapFactor_p )
                                                                 + fields[axisVel].chip(fieldEndIndex-1, axis) * fields[axisVel].chip(fieldEndIndex-1, axis).constant( extrapFactor_a );
                break;

            default:
                break;
        }


        // Axis negative boundary
        faceEndIndex = 0;
        fieldEndIndex = 0;
        switch ( boundaryConditions[axisVel][negativePatch].type ) {
            
            case BC::zeroGradient:
                faceVelocities[axisVel].chip(faceEndIndex, axis) = fields[axisVel].chip(fieldEndIndex, axis);          
                break;

            case BC::uniform:
                faceVelocities[axisVel].chip(faceEndIndex, axis) = faceVelocities[axisVel].chip(faceEndIndex, axis).constant( boundaryConditions[axisVel][negativePatch].value );
                break;

            case BC::extrapolated:
                extrapFactor_p = mesh.extrapFactors[negativePatch].p;
                extrapFactor_a = mesh.extrapFactors[negativePatch].a;
                faceVelocities[axisVel].chip(faceEndIndex, axis) = fields[axisVel].chip(fieldEndIndex, axis) * fields[axisVel].chip(fieldEndIndex, axis).constant( extrapFactor_p )
                                                                 + fields[axisVel].chip(fieldEndIndex+1, axis) * fields[axisVel].chip(fieldEndIndex-1, axis).constant( extrapFactor_a );
                break;

            default:
                break;
        }
    }

}

} // end namespace CFD