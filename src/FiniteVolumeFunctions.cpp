#include "FiniteVolumeFunctions.h"
#include "Tensor"


namespace
{

    // Return linear interpolation factor for cell faces
    CFD::floatType InterpFactor( const CFD::array1D &cellFaces, const CFD::array1D &cellCenters, const CFD::iterType index) 
    {
        return ( cellFaces(index) - cellCenters(index-1) ) / ( cellCenters(index) - cellCenters(index-1) );
    }


    void UpdateFaceVelocites_x(CFD::array3D &faceVel, const CFD::array3D &cellVel, const CFD::Mesh &mesh)
    {
        using namespace CFD;

        CFD::floatType interpFactor;
        for (iterType k = 0; k != faceVel.dimension(2); k++ ) {
            for (iterType j = 0; j != faceVel.dimension(1); j++) {
                for (iterType i = 1; i != faceVel.dimension(0)-1; i++) {

                    // Linear interpolation
                    interpFactor = InterpFactor(mesh.cellFaces_x, mesh.cellCenters_x, i);
                    faceVel(i, j, k) = (1 - interpFactor)*cellVel(i-1, j, k) + interpFactor*cellVel(i, j, k);

                }
            }
        }
    }


    void UpdateFaceVelocites_y(CFD::array3D &faceVel, const CFD::array3D &cellVel, const CFD::Mesh &mesh)
    {
        using namespace CFD;

        CFD::floatType interpFactor;
        for (iterType k = 0; k != faceVel.dimension(2); k++ ) {
            for (iterType j = 1; j != faceVel.dimension(1)-1; j++) {
                for (iterType i = 0; i != faceVel.dimension(0); i++) {

                    // Linear interpolation
                    interpFactor = InterpFactor(mesh.cellFaces_x, mesh.cellCenters_x, i);
                    faceVel(i, j, k) = (1 - interpFactor)*cellVel(i, j-1, k) + interpFactor*cellVel(i, j, k);

                }
            }
        }
    }


    void UpdateFaceVelocites_z(CFD::array3D &faceVel, const CFD::array3D &cellVel, const CFD::Mesh &mesh)
    {
        using namespace CFD;

        CFD::floatType interpFactor;
        for (iterType k = 1; k != faceVel.dimension(2)-1; k++ ) {
            for (iterType j = 0; j != faceVel.dimension(1); j++) {
                for (iterType i = 0; i != faceVel.dimension(0); i++) {

                    // Linear interpolation
                    interpFactor = InterpFactor(mesh.cellFaces_x, mesh.cellCenters_x, i);
                    faceVel(i, j, k) = (1 - interpFactor)*cellVel(i, j, k-1) + interpFactor*cellVel(i, j, k);

                }
            }
        }
    }


}   // end anonymous namespace


void CFD::UpdateFaceVelocities( ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const Mesh &mesh, const ArrayAllocator<CFD::Fields::ENUMDATA> &fields, 
    const std::vector< std::vector<InputData::BoundaryConditionStruct> > &boundaryConditions)
{

    using F = CFD::Fields::ENUMDATA;
    using BC = CFD::BoundaryConditions::ENUMDATA;
    using BP = CFD::BoundaryPatches::ENUMDATA;

    // Non-boundary cells
    UpdateFaceVelocites_x( faceVelocities[F::U], fields[F::U], mesh);
    UpdateFaceVelocites_y( faceVelocities[F::V], fields[F::V], mesh);
    UpdateFaceVelocites_z( faceVelocities[F::W], fields[F::W], mesh);


    // +x boundary
    switch ( boundaryConditions[F::U][BP::xPositive].type ) {
        case BC::zeroGradient:
            break;

        case BC::uniform:
            break;

        case BC::extrapolated:
            break;

        default:
            break;
    }


    // -x boundary
    switch ( boundaryConditions[F::U][BP::xNegative].type ) {
        case BC::zeroGradient:
            break;

        case BC::uniform:
            break;

        case BC::extrapolated:
            break;

        default:
            break;
    }


    // +y boundary
    switch ( boundaryConditions[F::V][BP::yPositive].type ) {
        case BC::zeroGradient:
            break;

        case BC::uniform:
            break;

        case BC::extrapolated:
            break;

        default:
            break;
    }


    // -y boundary
    switch ( boundaryConditions[F::V][BP::yNegative].type ) {
        case BC::zeroGradient:
            break;

        case BC::uniform:
            break;

        case BC::extrapolated:
            break;

        default:
            break;
    }


    // +z boundary
    switch ( boundaryConditions[F::W][BP::zPositive].type ) {
        case BC::zeroGradient:
            break;

        case BC::uniform:
            break;

        case BC::extrapolated:
            break;

        default:
            break;
    }


    // -z boundary
    switch ( boundaryConditions[F::W][BP::zNegative].type ) {
        case BC::zeroGradient:
            break;

        case BC::uniform:
            break;

        case BC::extrapolated:
            break;

        default:
            break;
    }

}