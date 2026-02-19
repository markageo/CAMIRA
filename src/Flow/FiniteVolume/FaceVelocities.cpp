#include "FiniteVolume.h"

#include "Core/FVTools.h"
#include "Core/FVLookups.h"
#include "FaceInterpolatedVelocity.h"


namespace CAMIRA
{

using namespace CORE;
using namespace FVT;
namespace FLOW
{

namespace
{

void LinearInterpInteriorFaceVelocitiesWithMWI( EnumVector<Axis, Tensor3D> &faceVelocities, 
                                                const FieldData<Tensor3D> &cellFields, 
                                                const FVCoefficients &fvCoeffs,
                                                const Mesh &mesh, 
                                                const Axis::ENUMDATA axis,
                                                const Axis::ENUMDATA velocityComponent )
{
    using enum Axis::ENUMDATA;
    using FVT::G;

    Tensor3D &faceVel = faceVelocities[ axis ];
    const Tensor3D &cellVel                         = cellFields.U[ velocityComponent ];
    const Tensor3D &cellPressure                    = cellFields.P;
    const Tensor3D &momentumDiagCoeff               = fvCoeffs.Mom[axis].AU[TransportCoefficients::p];
    const std::array<Tensor1D, 4> &mwiSparseCoeffs  = fvCoeffs.mwiSparseCoeffs[axis];
    const std::array<Tensor1D, 2> &mwiCompactCoeffs = fvCoeffs.mwiCompactCoeffs[axis];

    iArray3 startIndex = {0, 0, 0}; startIndex(axis)++;
    iArray3 nFaces = mesh.nFacesNormal[axis]; nFaces(axis)--;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++ ) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D idx = {i, j, k},
                              hiIndex = idx,
                              loIndex = idx,
                              loloIndex = idx,
                              hihiIndex = idx;
                loIndex[axis]   -= 1;
                loloIndex[axis] -= 2;
                hihiIndex[axis] += 1;

                faceVel(idx) = FaceInterpolatedVelocity<AdvectionSchemes::Central, +1>( cellVel, fvCoeffs, mesh, axis, hiIndex, loIndex );

                // Add MWI correction
                floatType d = ( 1.0f - mesh.interpFactors[axis](idx[axis]) ) * ( 1.0f / momentumDiagCoeff( G(loIndex) ) )  
                            +  mesh.interpFactors[axis](idx[axis])           *  (1.0f / momentumDiagCoeff( G(hiIndex) ) );
                floatType coeff0 = d *   mwiSparseCoeffs[0]( idx[axis] ),
                          coeff1 = d * ( mwiSparseCoeffs[1]( idx[axis] ) + mwiCompactCoeffs[0]( idx[axis] ) ),
                          coeff2 = d * ( mwiSparseCoeffs[2]( idx[axis] ) + mwiCompactCoeffs[1]( idx[axis] ) ),
                          coeff3 = d *   mwiSparseCoeffs[3]( idx[axis] );

                faceVel( idx ) += coeff0 * cellPressure( G(loloIndex) )
                                + coeff1 * cellPressure( G(loIndex) )
                                + coeff2 * cellPressure( G(hiIndex) )
                                + coeff3 * cellPressure( G(hihiIndex) );

            }
        }
    }
}


[[maybe_unused]]
void LinearInterpInteriorFaceVelocities( EnumVector<Axis, Tensor3D> &faceVelocities, 
                                         const Tensor3D &cellVelocities, 
                                         const Mesh &mesh, 
                                         const Axis::ENUMDATA axis)
{
    using enum Axis::ENUMDATA;

    Tensor3D &faceVel = faceVelocities[ axis ];

    iArray3 startIndex = {0, 0, 0}; startIndex(axis)++;
    iArray3 nFaces = mesh.nFacesNormal[axis]; nFaces(axis)--;

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++ ) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D idx = {i, j, k},
                              hiIndexG = G(idx),
                              loIndexG = G(idx);
                loIndexG[axis] -= 1;
                
                const floatType interpFactor = mesh.interpFactors[ axis ]( idx[axis] );
                faceVel( idx ) = (1 - interpFactor) * cellVelocities( loIndexG ) + interpFactor*cellVelocities( hiIndexG );

            }
        }
    }
}


// Performance optimised version
template<Axis::ENUMDATA axis>
void LinearInterpInteriorFaceVelocities2( EnumVector<Axis, Tensor3D> &faceVelocities, 
                                          const Tensor3D &cellVelocities, 
                                          const Mesh &mesh)
{
    using enum Axis::ENUMDATA;

    Tensor3D &faceVel = faceVelocities[ axis ];

    iArray3 startIndex = {0, 0, 0}; startIndex(axis)++;
    iArray3 nFaces = mesh.nFacesNormal[axis]; nFaces(axis)--;

    floatType interpFactor;
    TensorIndex3D hiIndexG, loIndexG;


    for (intType k = startIndex[Z]; k != nFaces[Z]; k++ ) {

        hiIndexG[Z] = G(k);
        loIndexG[Z] = G(k);

        if constexpr ( axis == Z ) {
            interpFactor = mesh.interpFactors[ axis ]( k );
            loIndexG[axis] -= 1;
        }

        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {

            hiIndexG[Y] = G(j);
            loIndexG[Y] = G(j);

            if constexpr ( axis == Y ) {
                interpFactor = mesh.interpFactors[ axis ]( j );
                loIndexG[axis] -= 1;
            }

            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                hiIndexG[X] = G(i);
                loIndexG[X] = G(i);

                if constexpr ( axis == X ) {
                    interpFactor = mesh.interpFactors[ axis ]( i );
                    loIndexG[axis] -= 1;
                }
                
                faceVel( i, j, k ) = (1 - interpFactor) * cellVelocities( loIndexG ) + interpFactor*cellVelocities( hiIndexG );

            }
        }
    }
}



void BoundaryFaceVelocities( EnumVector<Axis, Tensor3D> &faceVelocities, 
                            const EnumVector<Axis, Tensor3D> &cellVelocities, 
                            const Mesh &mesh, 
                            const EnumVector< Axis, BoundaryConditionData::Patches > &momentumBoundaryPatchData,
                            const BoundaryPatches::ENUMDATA boundaryPatch,
                            const Axis::ENUMDATA velocityComponent )
{
    using BC = BoundaryConditions::ENUMDATA;
    
    Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[ boundaryPatch ];
    
    static constexpr TensorIndex3D offsets = {nGhost, nGhost, nGhost};
    TensorIndex3D extents = {mesh.nCells(Axis::X), mesh.nCells(Axis::Y), mesh.nCells(Axis::Z)};
    
    intType faceEndIndex, fieldEndIndex;
    if ( boundaryPatch == LUT::PositivePatch[ axis ] ) {
        faceEndIndex = mesh.nCells(axis);
        fieldEndIndex = mesh.nCells(axis)-1;
    } else {
        faceEndIndex = 0;
        fieldEndIndex = 0;
    }

    switch ( momentumBoundaryPatchData[velocityComponent][boundaryPatch].type ) 
    {    
        case BC::zeroGradient:
        {
            faceVelocities[axis].chip(faceEndIndex, axis) = cellVelocities[velocityComponent].slice(offsets, extents).chip(fieldEndIndex, axis);          
            break;
        }
            

        case BC::fixed:
        {
            faceVelocities[axis].chip(faceEndIndex, axis) = momentumBoundaryPatchData[velocityComponent][boundaryPatch].value;
            break;
        }
            

        case BC::extrapolated:
        {
            floatType extrapFactor_p = mesh.extrapFactors[boundaryPatch].p;
            floatType extrapFactor_a = mesh.extrapFactors[boundaryPatch].a;
            intType fieldEndIndex_a  = ( boundaryPatch == LUT::PositivePatch[ axis ] ) ? fieldEndIndex - 1 : fieldEndIndex + 1 ;
            auto cellField_p = cellVelocities[velocityComponent].slice(offsets, extents).chip(fieldEndIndex  , axis);
            auto cellField_a = cellVelocities[velocityComponent].slice(offsets, extents).chip(fieldEndIndex_a, axis);
            faceVelocities[axis].chip(faceEndIndex, axis) = cellField_p * cellField_p.constant( extrapFactor_p )
                                                          + cellField_a * cellField_a.constant( extrapFactor_a );
            break;
        }


        case BC::periodic:
        {
            intType loCellIndex = 0,
                    hiCellIndex = mesh.nCells(axis)-1;
            floatType interpFactor = mesh.cellLengths[axis](loCellIndex) / ( mesh.cellLengths[axis](loCellIndex) + mesh.cellLengths[axis](hiCellIndex) );
            auto cellFieldHi = cellVelocities[velocityComponent].slice(offsets, extents).chip(hiCellIndex, axis);
            auto cellFieldLo = cellVelocities[velocityComponent].slice(offsets, extents).chip(loCellIndex, axis);
            faceVelocities[axis].chip(faceEndIndex, axis) = cellFieldLo * cellFieldLo.constant( 1.0f - interpFactor )
                                                          + cellFieldHi * cellFieldHi.constant( interpFactor );
        }
            

        default:
            break;
    }


}

}   // end anonymous namespace



// ---------------------------------------- Face Fluxes ----------------------------------------

// Faces are staggered in the negative direction:
//   cellFaceFlux[X](i, j, k) -> u(i-1/2, j    , k    )
//   cellFaceFlux[Y](i, j, k) -> u(i    , j-1/2, k    )
//   cellFaceFlux[Z](i, j, k) -> u(i    , j    , k-1/2)
// Subscript indicates the normal direction of the face.


// Calculates face velocity fluxes. i.e. normal component of velocity on faces
void UpdateFaceFluxes( EnumVector< Axis, Tensor3D > &faceFluxes, 
                       const Mesh &mesh, 
                       const EnumVector< Axis, Tensor3D > &cellVelocities, 
                       const BoundaryConditionData &bcData )
{
    // Internal faces
    // EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
    //     LinearInterpInteriorFaceVelocities( faceFluxes, cellVelocities[axis], mesh, axis);
    // } );
    using enum Axis::ENUMDATA;
    LinearInterpInteriorFaceVelocities2<X>( faceFluxes, cellVelocities[X], mesh);
    LinearInterpInteriorFaceVelocities2<Y>( faceFluxes, cellVelocities[Y], mesh);
    LinearInterpInteriorFaceVelocities2<Z>( faceFluxes, cellVelocities[Z], mesh);

    // Boundary faces
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA boundaryPatch) {

        Axis::ENUMDATA velocityComponent = LUT::BoundaryPatchAxis[ boundaryPatch ];
        BoundaryFaceVelocities( faceFluxes, cellVelocities, mesh, bcData.fields.U, boundaryPatch, velocityComponent );

    } );
}


// Calculates face velocity fluxes. i.e. normal component of velocity on faces with MWI correction
void UpdateFaceFluxesWithMWI( EnumVector< Axis, Tensor3D > &faceFluxes, 
                              const Mesh &mesh, 
                              const FieldData<Tensor3D> &cellFields,
                              const FVCoefficients &fvCoeffs, 
                              const BoundaryConditionData &bcData )
{
    // Internal faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        LinearInterpInteriorFaceVelocitiesWithMWI( faceFluxes, cellFields, fvCoeffs, mesh, axis, axis);
    } );
    
    // Boundary faces
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA boundaryPatch) {

        Axis::ENUMDATA velocityComponent = LUT::BoundaryPatchAxis[ boundaryPatch ];
        BoundaryFaceVelocities( faceFluxes, cellFields.U, mesh, bcData.fields.U, boundaryPatch, velocityComponent );

    } );
}


// ------------------------------------ Immersed Boundary Face Fluxes ------------------------------------

void SetIBFaceFluxes( EnumVector<Axis, Tensor3D> &faceFluxes,
                      const IBData &ibData ) 
{
    for ( auto &ibCellComponent : ibData.ibCells ) {

        for ( auto &ibCell : ibCellComponent ) { 
            for ( auto &sourceTermData : ibCell.sourceTermsData ) {

                Axis::ENUMDATA axis = sourceTermData.direction;
                TensorIndex3D faceIndex = ibCell.cellIndex;    
                faceIndex[axis] += sourceTermData.faceDirectionIndex;

                faceFluxes[axis](faceIndex) = sourceTermData.faceValues.U[axis];
            }
        }
    }
}





}   // end namespace FLOW

}   // end namespace CAMIRA