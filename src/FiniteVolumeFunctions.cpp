#include "FiniteVolumeFunctions.h"
#include "Tensor"

#include <iostream>


namespace
{

    /* ------------------------------------------ Face velocity interpolation ------------------------------------------*/

    void UpdateFaceVelocites_x(CFD::array3D &faceVel, const CFD::array3D &cellVel, const CFD::Mesh &mesh)
    {
        using namespace CFD;

        floatType interpFactor;
        for (iterType k = 0; k != faceVel.dimension(2); k++ ) {
            for (iterType j = 0; j != faceVel.dimension(1); j++) {
                for (iterType i = 1; i != faceVel.dimension(0)-1; i++) {

                    // Linear interpolation
                    interpFactor = mesh.interpFactors_x(i);
                    faceVel(i, j, k) = (1 - interpFactor)*cellVel(i-1, j, k) + interpFactor*cellVel(i, j, k);

                }
            }
        }
    }


    void UpdateFaceVelocites_y(CFD::array3D &faceVel, const CFD::array3D &cellVel, const CFD::Mesh &mesh)
    {
        using namespace CFD;

        floatType interpFactor;
        for (iterType k = 0; k != faceVel.dimension(2); k++ ) {
            for (iterType j = 1; j != faceVel.dimension(1)-1; j++) {
                for (iterType i = 0; i != faceVel.dimension(0); i++) {

                    // Linear interpolation
                    interpFactor = mesh.interpFactors_y(j);
                    faceVel(i, j, k) = (1 - interpFactor)*cellVel(i, j-1, k) + interpFactor*cellVel(i, j, k);

                }
            }
        }
    }


    void UpdateFaceVelocites_z(CFD::array3D &faceVel, const CFD::array3D &cellVel, const CFD::Mesh &mesh)
    {
        using namespace CFD;

        floatType interpFactor;
        for (iterType k = 1; k != faceVel.dimension(2)-1; k++ ) {
            for (iterType j = 0; j != faceVel.dimension(1); j++) {
                for (iterType i = 0; i != faceVel.dimension(0); i++) {

                    // Linear interpolation
                    interpFactor = mesh.interpFactors_z(k);
                    faceVel(i, j, k) = (1 - interpFactor)*cellVel(i, j, k-1) + interpFactor*cellVel(i, j, k);

                }
            }
        }
    }



    /* ------------------------------------------ Boundary conditions on x normal faces ------------------------------------------*/

    void setZeroGradient_x( CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &fields, 
                            const CFD::iterType faceIndex, const CFD::iterType fieldIndex) 
    {
        using namespace CFD;
        using F = CFD::Fields::ENUMDATA;

        for (iterType k = 0; k != faceVelocities[F::U].dimension(2); k++ ) {
            for (iterType j = 0; j != faceVelocities[F::U].dimension(1); j++) {
                faceVelocities[F::U](faceIndex, j, k) = fields[F::U](fieldIndex, j, k);
            }
        }   
    }


    void setUniform_x( CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const CFD::iterType faceIndex, const CFD::floatType value) 
    {
        using namespace CFD;
        using F = CFD::Fields::ENUMDATA;

        for (iterType k = 0; k != faceVelocities[F::U].dimension(2); k++ ) {
            for (iterType j = 0; j != faceVelocities[F::U].dimension(1); j++) {
                faceVelocities[F::U](faceIndex, j, k) = value;
            }
        }   
    }


    void setExtrapolated_x( CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &fields, const CFD::Mesh &mesh, 
                            const CFD::iterType faceIndex, const CFD::iterType fieldIndex_p, const CFD::iterType fieldIndex_f) 
    {
        using namespace CFD;
        using F = CFD::Fields::ENUMDATA;

        floatType extrapFactor_p = ( 2.0*mesh.cellLengths_x(fieldIndex_p) + mesh.cellLengths_x(fieldIndex_f) )
                           / ( mesh.cellLengths_x(fieldIndex_p) + mesh.cellLengths_x(fieldIndex_f) );
        floatType extrapFactor_f = extrapFactor_p - ( mesh.cellLengths_x(fieldIndex_p) + mesh.cellLengths_x(fieldIndex_f) );

        for (iterType k = 0; k != faceVelocities[F::U].dimension(2); k++ ) {
            for (iterType j = 0; j != faceVelocities[F::U].dimension(1); j++) {
                faceVelocities[F::U](faceIndex, j, k) = extrapFactor_p*fields[F::U](fieldIndex_p, j, k)
                                                            - extrapFactor_f*fields[F::U](fieldIndex_f, j, k);
            }
        } 
    }



    /* ------------------------------------------ Boundary conditions on y normal faces ------------------------------------------*/

    void setZeroGradient_y( CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &fields, 
                            const CFD::iterType faceIndex, const CFD::iterType fieldIndex) 
    {
        using namespace CFD;
        using F = CFD::Fields::ENUMDATA;

        for (iterType k = 0; k != faceVelocities[F::V].dimension(2); k++ ) {
            for (iterType i = 0; i != faceVelocities[F::V].dimension(0); i++) {
                faceVelocities[F::V](i, faceIndex, k) = fields[F::V](i, fieldIndex, k);
            }
        }   
    }

    void setUniform_y( CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const CFD::iterType faceIndex, const CFD::floatType value) 
    {
        using namespace CFD;
        using F = CFD::Fields::ENUMDATA;

        for (iterType k = 0; k != faceVelocities[F::V].dimension(2); k++ ) {
            for (iterType i = 0; i != faceVelocities[F::V].dimension(0); i++) {
                faceVelocities[F::V](i, faceIndex, k) = value;
            }
        }   
    }

    void setExtrapolated_y( CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &fields, const CFD::Mesh &mesh, 
                            const CFD::iterType faceIndex, const CFD::iterType fieldIndex_p, const CFD::iterType fieldIndex_f) 
    {
        using namespace CFD;
        using F = CFD::Fields::ENUMDATA;

        floatType extrapFactor_p = ( 2.0*mesh.cellLengths_y(fieldIndex_p) + mesh.cellLengths_y(fieldIndex_f) )
                           / ( mesh.cellLengths_y(fieldIndex_p) + mesh.cellLengths_y(fieldIndex_f) );
        floatType extrapFactor_f = extrapFactor_p - ( mesh.cellLengths_y(fieldIndex_p) + mesh.cellLengths_y(fieldIndex_f) );

        for (iterType k = 0; k != faceVelocities[F::V].dimension(2); k++ ) {
            for (iterType i = 0; i != faceVelocities[F::V].dimension(0); i++) {
                faceVelocities[F::V](i, faceIndex, k) = extrapFactor_p*fields[F::V](i, fieldIndex_p, k)
                                                            - extrapFactor_f*fields[F::V](i, fieldIndex_f, k);
            }
        } 
    }



    /* ------------------------------------------ Boundary conditions on z normal faces ------------------------------------------*/

    void setZeroGradient_z( CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &fields, 
                            const CFD::iterType faceIndex, const CFD::iterType fieldIndex) 
    {
        using namespace CFD;
        using F = CFD::Fields::ENUMDATA;

        for (iterType j = 0; j != faceVelocities[F::W].dimension(1); j++ ) {
            for (iterType i = 0; i != faceVelocities[F::W].dimension(0); i++) {
                faceVelocities[F::W](i, j, faceIndex) = fields[F::W](i, j, fieldIndex);
            }
        }   
    }

    void setUniform_z( CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const CFD::iterType faceIndex, const CFD::floatType value) 
    {
        using namespace CFD;
        using F = CFD::Fields::ENUMDATA;

        for (iterType j = 0; j != faceVelocities[F::W].dimension(1); j++ ) {
            for (iterType i = 0; i != faceVelocities[F::W].dimension(0); i++) {
                faceVelocities[F::W](i, j, faceIndex) = value;
            }
        }   
    }

    void setExtrapolated_z( CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const CFD::ArrayAllocator<CFD::Fields::ENUMDATA> &fields, const CFD::Mesh &mesh, 
                            const CFD::iterType faceIndex, const CFD::iterType fieldIndex_p, const CFD::iterType fieldIndex_f) 
    {
        using namespace CFD;
        using F = CFD::Fields::ENUMDATA;

        floatType extrapFactor_p = ( 2.0*mesh.cellLengths_z(fieldIndex_p) + mesh.cellLengths_z(fieldIndex_f) )
                           / ( mesh.cellLengths_z(fieldIndex_p) + mesh.cellLengths_z(fieldIndex_f) );
        floatType extrapFactor_f = extrapFactor_p - ( mesh.cellLengths_z(fieldIndex_p) + mesh.cellLengths_z(fieldIndex_f) );

        for (iterType j = 0; j != faceVelocities[F::W].dimension(1); j++ ) {
            for (iterType i = 0; i != faceVelocities[F::W].dimension(0); i++) {
                faceVelocities[F::W](i, j, faceIndex) = extrapFactor_p*fields[F::W](i, j, fieldIndex_p)
                                                            - extrapFactor_f*fields[F::W](i, j, fieldIndex_f);
            }
        } 
    }


}   // end anonymous namespace


// Calculate cell normal face velocities using linear interpolation from cell centers
void CFD::UpdateFaceVelocities( ArrayAllocator<CFD::Fields::ENUMDATA> &faceVelocities, const CFD::Mesh &mesh, const ArrayAllocator<CFD::Fields::ENUMDATA> &fields, 
    const std::vector< std::vector<InputData::BoundaryConditionStruct> > &boundaryConditions)
{

    using F = CFD::Fields::ENUMDATA;
    using BC = CFD::BoundaryConditions::ENUMDATA;
    using BP = CFD::BoundaryPatches::ENUMDATA;

    // Non-boundary cells
    UpdateFaceVelocites_x( faceVelocities[F::U], fields[F::U], mesh);
    UpdateFaceVelocites_y( faceVelocities[F::V], fields[F::V], mesh);
    UpdateFaceVelocites_z( faceVelocities[F::W], fields[F::W], mesh);

    // Just to assign the boundary index of the face flux and field arrays
    iterType faceEndIndex, fieldEndIndex;

    // +x boundary
    faceEndIndex = faceVelocities[F::U].dimension(0)-1;
    fieldEndIndex = fields[F::U].dimension(0)-1;
    switch ( boundaryConditions[F::U][BP::xPositive].type ) {
        
        case BC::zeroGradient:
            setZeroGradient_x(faceVelocities, fields, faceEndIndex, fieldEndIndex);          
            break;

        case BC::uniform:
            setUniform_x(faceVelocities, faceEndIndex, boundaryConditions[F::U][BP::xPositive].value);
            break;

        case BC::extrapolated:
            setExtrapolated_x(faceVelocities, fields, mesh, faceEndIndex, fieldEndIndex, fieldEndIndex-1);
            break;

        default:
            break;
    }


    // -x boundary
    faceEndIndex = 0;
    fieldEndIndex = 0;
    switch ( boundaryConditions[F::U][BP::xNegative].type ) {
        case BC::zeroGradient:
            setZeroGradient_x(faceVelocities, fields, faceEndIndex, fieldEndIndex); 
            break;

        case BC::uniform:
            setUniform_x(faceVelocities, faceEndIndex, boundaryConditions[F::U][BP::xNegative].value);
            break;

        case BC::extrapolated:
            setExtrapolated_x(faceVelocities, fields, mesh, faceEndIndex, fieldEndIndex, fieldEndIndex+1);
            break;

        default:
            break;
    }


    // +y boundary
    faceEndIndex = faceVelocities[F::V].dimension(1)-1;
    fieldEndIndex = fields[F::V].dimension(1)-1;
    switch ( boundaryConditions[F::V][BP::yPositive].type ) {
        case BC::zeroGradient:
            setZeroGradient_y(faceVelocities, fields, faceEndIndex, fieldEndIndex); 
            break;

        case BC::uniform:
            setUniform_y(faceVelocities, faceEndIndex, boundaryConditions[F::V][BP::yPositive].value);
            break;

        case BC::extrapolated:
            setExtrapolated_y(faceVelocities, fields, mesh, faceEndIndex, fieldEndIndex, fieldEndIndex-1);
            break;

        default:
            break;
    }


    // -y boundary
    faceEndIndex = 0;
    fieldEndIndex = 0;
    switch ( boundaryConditions[F::V][BP::yNegative].type ) {
        case BC::zeroGradient:
            setZeroGradient_y(faceVelocities, fields, faceEndIndex, fieldEndIndex); 
            break;

        case BC::uniform:
            setUniform_y(faceVelocities, faceEndIndex, boundaryConditions[F::V][BP::yNegative].value);
            break;

        case BC::extrapolated:
            setExtrapolated_y(faceVelocities, fields, mesh, faceEndIndex, fieldEndIndex, fieldEndIndex+1);
            break;

        default:
            break;
    }


    // +z boundary
    faceEndIndex = faceVelocities[F::W].dimension(2)-1;
    fieldEndIndex = fields[F::W].dimension(2)-1;
    switch ( boundaryConditions[F::W][BP::zPositive].type ) {
        case BC::zeroGradient:
            setZeroGradient_z(faceVelocities, fields, faceEndIndex, fieldEndIndex); 
            break;

        case BC::uniform:
            setUniform_z(faceVelocities, faceEndIndex, boundaryConditions[F::W][BP::zPositive].value);
            break;

        case BC::extrapolated:
            setExtrapolated_z(faceVelocities, fields, mesh, faceEndIndex, fieldEndIndex, fieldEndIndex-1);
            break;

        default:
            break;
    }


    // -z boundary
    faceEndIndex = 0;
    fieldEndIndex = 0;
    switch ( boundaryConditions[F::W][BP::zNegative].type ) {
        case BC::zeroGradient:
            setZeroGradient_z(faceVelocities, fields, faceEndIndex, fieldEndIndex); 
            break;

        case BC::uniform:
            setUniform_z(faceVelocities, faceEndIndex, boundaryConditions[F::W][BP::zNegative].value);
            break;

        case BC::extrapolated:
            setExtrapolated_z(faceVelocities, fields, mesh, faceEndIndex, fieldEndIndex, fieldEndIndex+1);
            break;

        default:
            break;
    }

}