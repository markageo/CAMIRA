#include "FiniteVolume.h"
#include "../Tools/FVTools.h"
#include "../Tools/FVLookups.h"

namespace CFD
{

using namespace FVT;

namespace
{

void LinearInterpInteriorFaceVelocityWithMWI( EnumVector<Axis, Tensor3D> &faceVelocities, 
                                              const FieldData<Tensor3D> &cellFields, 
                                              const FVCoefficients &fvCoeffs,
                                              const Mesh &mesh, 
                                              const Axis::ENUMDATA axis,
                                              const Axis::ENUMDATA velocityComponent )
{
    using enum Axis::ENUMDATA;

    Tensor3D &faceVel = faceVelocities[ axis ];
    const Tensor3D &cellVel = cellFields.U[ velocityComponent ];
    const Tensor3D &cellPressure = cellFields.P;
    const Tensor3D &momentumDiagCoeffInv = fvCoeffs.Mom[axis].diagCoeffInv;
    const std::array<Tensor1D, 4> &mwiSparseCoeffs  = fvCoeffs.Cont.mwiSparseCoeffs[axis];
    const std::array<Tensor1D, 2> &mwiCompactCoeffs = fvCoeffs.Cont.mwiCompactCoeffs[axis];

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++ ) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D idx = {i, j, k},
                             HiIndex = idx,
                             LoIndex = idx,
                             LoLoIndex = idx,
                             HiHiIndex = idx;
                LoIndex[axis]   -= 1;
                LoLoIndex[axis] -= 2;
                HiHiIndex[axis] += 1;

                floatType interpFactor = mesh.interpFactors[ axis ]( idx[axis] );
                faceVel( idx ) = (1 - interpFactor)*cellVel( G(LoIndex) ) + interpFactor*cellVel( G(HiIndex) );

                // Add MWI correction
                floatType d = 0.5f * ( momentumDiagCoeffInv( HiIndex )  +  momentumDiagCoeffInv( LoIndex ) );
                floatType coeff0 = d *   mwiSparseCoeffs[0]( idx[axis] ),
                          coeff1 = d * ( mwiSparseCoeffs[1]( idx[axis] ) + mwiCompactCoeffs[0]( idx[axis] ) ),
                          coeff2 = d * ( mwiSparseCoeffs[2]( idx[axis] ) + mwiCompactCoeffs[1]( idx[axis] ) ),
                          coeff3 = d *   mwiSparseCoeffs[3]( idx[axis] );

                faceVel( idx ) += coeff0 * cellPressure( G(LoLoIndex) )
                                + coeff1 * cellPressure( G(LoIndex) )
                                + coeff2 * cellPressure( G(HiIndex) )
                                + coeff3 * cellPressure( G(HiHiIndex) );

            }
        }
    }
}



void LinearInterpInteriorFaceVelocity( EnumVector<Axis, Tensor3D> &faceVelocities, 
                                       const Tensor3D &cellVelocities, 
                                       const Mesh &mesh, 
                                       const Axis::ENUMDATA axis)
{
    using enum Axis::ENUMDATA;

    Tensor3D &faceVel = faceVelocities[ axis ];

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++ ) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D idx = {i, j, k},
                             HiIndex = idx,
                             LoIndex = idx;
                LoIndex[axis] -= 1;

                floatType interpFactor = mesh.interpFactors[ axis ]( idx[axis] );
                faceVel( idx ) = (1 - interpFactor)*cellVelocities( G(LoIndex) ) + interpFactor*cellVelocities( G(HiIndex) );

            }
        }
    }
}




void UpwindInteriorFaceVelocity( EnumVector<Axis, Tensor3D> &faceVelocities, 
                                 const Tensor3D &cellVelocities, 
                                 const EnumVector<Axis, Tensor3D> &faceFluxes,
                                 const Mesh &mesh, 
                                 const Axis::ENUMDATA axis)
{
    using enum Axis::ENUMDATA;

    Tensor3D &faceVel = faceVelocities[ axis ];

    auto [startIndex, nFaces] = FaceInternalIndices(mesh, axis);

    for (intType k = startIndex[Z]; k != nFaces[Z]; k++ ) {
        for (intType j = startIndex[Y]; j != nFaces[Y]; j++) {
            for (intType i = startIndex[X]; i != nFaces[X]; i++) {

                TensorIndex3D idx = {i, j, k};

                if ( faceFluxes[axis](idx) >= 0.0f ) {
                    idx[axis] -= 1;
                    faceVel( idx ) = cellVelocities( G(idx) );
                } else {
                    faceVel( idx ) = cellVelocities( G(idx) );
                }

            }
        }
    }
}





void BoundaryFaceVelocitiy( EnumVector<Axis, Tensor3D> &faceVelocities, 
                            const EnumVector<Axis, Tensor3D> &cellVelocities, 
                            const Mesh &mesh, 
                            const EnumVector< Axis, BoundaryConditionData >& boundaryConditions,
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

    switch ( boundaryConditions[velocityComponent][boundaryPatch].type ) 
    {    
        case BC::zeroGradient:
        {
            faceVelocities[axis].chip(faceEndIndex, axis) = cellVelocities[velocityComponent].slice(offsets, extents).chip(fieldEndIndex, axis);          
            break;
        }
            

        case BC::fixed:
        {
            faceVelocities[axis].chip(faceEndIndex, axis) = boundaryConditions[velocityComponent][boundaryPatch].value;
            break;
        }
            

        case BC::extrapolated:
        {
            floatType extrapFactor_p = mesh.extrapFactors[boundaryPatch].p;
            floatType extrapFactor_a = mesh.extrapFactors[boundaryPatch].a;
            faceVelocities[axis].chip(faceEndIndex, axis) = cellVelocities[velocityComponent].slice(offsets, extents).chip(fieldEndIndex  , axis) 
                                                                * cellVelocities[velocityComponent].slice(offsets, extents).chip(fieldEndIndex  , axis).constant( extrapFactor_p )
                                                          + cellVelocities[velocityComponent].slice(offsets, extents).chip(fieldEndIndex+1, axis) 
                                                                * cellVelocities[velocityComponent].slice(offsets, extents).chip(fieldEndIndex+1, axis).constant( extrapFactor_a );
            break;
        }
            

        default:
            break;
    }


}

}   // end anonymous namespace


// ---------------------------------------- Face Advected Velocities ----------------------------------------

// Calculates advected face velocities
void UpdateFaceAdvectedVelocities( EnumVector< Axis, EnumVector< Axis, Tensor3D > > &faceAdvectedVelocities, 
                                   const Mesh &mesh, 
                                   const EnumVector< Axis, Tensor3D > &cellVelocities, 
                                   const EnumVector< Axis, Tensor3D > &faceFluxes,
                                   const FieldData< BoundaryConditionData > &bcData )
{
    // Internal faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        EnumFor<Axis>( [&] (Axis::ENUMDATA velocityComponent) {
            UpwindInteriorFaceVelocity( faceAdvectedVelocities[velocityComponent], cellVelocities[velocityComponent], faceFluxes, mesh, axis);
        } );
    } );
    
    // Boundary faces
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA boundaryPatch) {
        EnumFor<Axis>( [&] (Axis::ENUMDATA velocityComponent) {
            BoundaryFaceVelocitiy( faceAdvectedVelocities[velocityComponent], cellVelocities, mesh, bcData.U, boundaryPatch, velocityComponent );
        } );
    } );
}



EnumVector< Axis, EnumVector<Axis, Tensor3D> > InitialiseAdvectedFaceVelocities( const Mesh &mesh, 
                                                                                const EnumVector<Axis, Tensor3D> &cellVelocities, 
                                                                                const EnumVector<Axis, Tensor3D> &faceFluxes,
                                                                                const FieldData< BoundaryConditionData > &bcData)
{
    // First index is the velocity component. Second index is the face normal.
    // Faces are staggered in the negative direction:
    //   cellFaceVelocity[X][X](i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceVelocity[X][Y](i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceVelocity[X][Z](i, j, k) -> u(i    , j    , k-1/2)
    EnumVector< Axis, EnumVector<Axis, Tensor3D> > faceAdvectedVelocities( EnumVector<Axis, Tensor3D> ( {{Axis::X, {mesh.nCells(0) + 1, mesh.nCells(1)    , mesh.nCells(2)    }},
                                                                                                       {Axis::Y, {mesh.nCells(0)    , mesh.nCells(1) + 1, mesh.nCells(2)    }},
                                                                                                       {Axis::Z, {mesh.nCells(0)    , mesh.nCells(1)    , mesh.nCells(2) + 1}}} ) );

                                                     
    UpdateFaceAdvectedVelocities(faceAdvectedVelocities, mesh, cellVelocities, faceFluxes, bcData);

    return faceAdvectedVelocities;
}





// ---------------------------------------- Face Fluxes ----------------------------------------

// Calculates face velocity fluxes. i.e. normal component of velocity on faces
void UpdateFaceFluxes( EnumVector< Axis, Tensor3D > &faceFluxes, 
                       const Mesh &mesh, 
                       const EnumVector< Axis, Tensor3D > &cellVelocities, 
                       const FieldData< BoundaryConditionData > &bcData )
{
    // Internal faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        LinearInterpInteriorFaceVelocity( faceFluxes, cellVelocities[axis], mesh, axis);
    } );
    
    // Boundary faces
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA boundaryPatch) {

        Axis::ENUMDATA velocityComponent = LUT::BoundaryPatchAxis[ boundaryPatch ];
        BoundaryFaceVelocitiy( faceFluxes, cellVelocities, mesh, bcData.U, boundaryPatch, velocityComponent );

    } );
}


// Calculates face velocity fluxes. i.e. normal component of velocity on faces with MWI correction
void UpdateFaceFluxesWithMWI( EnumVector< Axis, Tensor3D > &faceFluxes, 
                              const Mesh &mesh, 
                              const FieldData<Tensor3D> &cellFields,
                              const FVCoefficients &fvCoeffs, 
                              const FieldData< BoundaryConditionData > &bcData )
{
    // Internal faces
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        LinearInterpInteriorFaceVelocityWithMWI( faceFluxes, cellFields, fvCoeffs, mesh, axis, axis);
    } );
    
    // Boundary faces
    EnumFor<BoundaryPatches>( [&] (BoundaryPatches::ENUMDATA boundaryPatch) {

        Axis::ENUMDATA velocityComponent = LUT::BoundaryPatchAxis[ boundaryPatch ];
        BoundaryFaceVelocitiy( faceFluxes, cellFields.U, mesh, bcData.U, boundaryPatch, velocityComponent );

    } );
}


EnumVector<Axis, Tensor3D> InitialiseFaceFluxes( const Mesh &mesh, 
                                                const EnumVector<Axis, Tensor3D> &cellVelocities, 
                                                const FieldData< BoundaryConditionData > &bcData)
{
    // Faces are staggered in the negative direction:
    //   cellFaceFlux[X](i, j, k) -> u(i-1/2, j    , k    )
    //   cellFaceFlux[Y](i, j, k) -> u(i    , j-1/2, k    )
    //   cellFaceFlux[Z](i, j, k) -> u(i    , j    , k-1/2)
    // Subscript indicates the normal direction of the face.
    EnumVector<Axis, Tensor3D> faceFluxes( {{Axis::X, {mesh.nCells(0) + 1, mesh.nCells(1)    , mesh.nCells(2)    }},
                                           {Axis::Y, {mesh.nCells(0)    , mesh.nCells(1) + 1, mesh.nCells(2)    }},
                                           {Axis::Z, {mesh.nCells(0)    , mesh.nCells(1)    , mesh.nCells(2) + 1}}} );
                                                     
    UpdateFaceFluxes(faceFluxes, mesh, cellVelocities, bcData);

    return faceFluxes;
}





// ------------------------------------ Immersed Boundary Face Fluxes ------------------------------------

void SetIBFaceFluxes( EnumVector<Axis, Tensor3D> &faceFluxes,
                      const IBData &ibData ) 
{
    using CFD::FVT::G;

    for ( auto &ibCell : ibData.ibCells ) { 
        for ( auto &sourceTermData : ibCell.sourceTermsData ) {

            Axis::ENUMDATA axis = sourceTermData.direction;
            TensorIndex3D faceIndex = ibCell.cellIndex;    
            faceIndex[axis] += sourceTermData.faceDirectionIndex;

            faceFluxes[axis](faceIndex) = sourceTermData.faceValues.U[axis];
        }
    }
}



} // end namespace CFD