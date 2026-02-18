#ifndef CAMIRA_GHOST_CELLS 
#define CAMIRA_GHOST_CELLS


#include "FiniteVolume.h"
#include "Core/Macros.h"
#include "Core/FVTools.h"
#include "Core/FVLookups.h"

#include <utility>

namespace CAMIRA
{
 
using namespace FVT;


// Set the ghost cells for all fields to give correct boundary condition
// Assumes ghost cells are the same dimension as the last interior cell
__attribute__((always_inline))
inline void SetGhostCells( FieldData<Tensor3D> &fields,
                    const Mesh &mesh,
                    const BoundaryConditionData &bcData )
{
    using BC = BoundaryConditions;

    static constexpr Eigen::array< std::pair< Eigen::Index, Eigen::Index >, 2 > ghostPaddings = { std::make_pair( nGhost, nGhost ), 
                                                                                                  std::make_pair( nGhost, nGhost ) };

    auto SetBoundaryPatchGhostCells = [&] (BoundaryPatches::ENUMDATA boundaryPatch) {
        Axis::ENUMDATA axis = LUT::BoundaryPatchAxis[ boundaryPatch ];

        intType iCell_p = ( boundaryPatch == LUT::NegativePatch[axis] ) ? 0  : mesh.nCells[axis] - 1,
                iCell_g = ( boundaryPatch == LUT::NegativePatch[axis] ) ? -1 : mesh.nCells[axis];

        ForAllFieldData( [&] ( intType f ) {

            switch ( bcData.fields[f][boundaryPatch].type ) {
                case BC::zeroGradient: 
                {
                    fields[f].chip( G(iCell_g), axis ) = fields[f].chip( G(iCell_p) , axis );
                    break;
                }

                case BC::fixed: 
                {
                    fields[f].chip( G(iCell_g), axis ) = - fields[f].chip( G(iCell_p), axis )
                                                       +   fields[f].chip( G(iCell_p), axis ).constant( 2.0f ) * bcData.fields[f][boundaryPatch].value.pad( ghostPaddings );
                    break;
                }
                    
                case BC::extrapolated: 
                {
                    intType iCell_a = ( boundaryPatch = LUT::NegativePatch[axis] ) ? 1 : mesh.nCells[axis] - 2;
                    
                    auto extrapolatedFaceValues = fields[f].chip( G(iCell_p), axis ) * fields[f].constant( mesh.extrapFactors[boundaryPatch].p )
                                                + fields[f].chip( G(iCell_a), axis ) * fields[f].constant( mesh.extrapFactors[boundaryPatch].a );
                    
                    fields[f].chip( G(iCell_g), axis ) = - fields[f].chip( G(iCell_p), axis )
                                                       +   fields[f].chip( G(iCell_p), axis ).constant( 2.0f ) * extrapolatedFaceValues;
                    break;
                }

                case BC::periodic:
                {
                    intType loCellIndex = 0,
                            hiCellIndex = mesh.nCells(axis) - 1;
                    floatType denominator = mesh.cellLengths[axis](loCellIndex) + mesh.cellLengths[axis](hiCellIndex);
                    floatType numerator   = ( boundaryPatch == LUT::PositivePatch[axis] ) ? 
                                            mesh.cellLengths[axis](loCellIndex) - mesh.cellLengths[axis](hiCellIndex) :
                                            mesh.cellLengths[axis](loCellIndex) + mesh.cellLengths[axis](loCellIndex) ; 
                    floatType interpFactor =  numerator / denominator;
                    auto cellFieldHi = fields[f].chip( G(hiCellIndex), axis );
                    auto cellFieldLo = fields[f].chip( G(loCellIndex), axis );
                    fields[f].chip( G(iCell_g), axis ) = cellFieldLo * cellFieldLo.constant( 1.0f - interpFactor )
                                                       + cellFieldHi * cellFieldHi.constant( interpFactor );
                }
                    
            }

        } );
    };


    // Each boundary patch will be done by a seperate thread
    #pragma omp parallel for
    for ( intType i = 0; i != BoundaryPatches::count; i++ ) {
        auto boundaryPatch = static_cast<BoundaryPatches::ENUMDATA>(i);
        SetBoundaryPatchGhostCells(boundaryPatch);
    }

}



}   // end namespace CAMIRA


#endif // CAMIRA_GHOST_CELLS
