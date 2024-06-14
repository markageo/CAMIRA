#include "Mesh.h"
#include "../Tools/FVTools.h"
#include "../Tools/FVLookups.h"

#include <cmath>
#include <utility>
#include <limits>

namespace CFD
{


namespace
{

std::vector<floatType> CalculateGrowthRates(const std::vector<InputData::MeshSegment> &meshSegments)
{
    size_t nSegments = meshSegments.size();
    std::vector<floatType> growthRates(nSegments);

    // Negative growth rate means shrinking grid
    for (size_t i = 0; i != nSegments; i++) {
        growthRates[i] = std::pow( std::abs(meshSegments[i].biasFactor) , 1.0f / static_cast<floatType>( meshSegments[i].nCells - 1 ) );
        if (meshSegments[i].biasFactor < 0 )   
            growthRates[i] = 1.0f/growthRates[i];
    }

    return growthRates;
}



void CalculateCellLengths(Tensor1D &cellLengths, 
                          const std::vector<InputData::MeshSegment> &meshSegments, 
                          const std::vector<floatType> &growthRates)
{
    size_t nSegments = meshSegments.size();

    floatType segmentLength, firstCellLength, geometricFactor;
    int cellIndex = 0;
    for (size_t s = 0; s != nSegments; s++) {    // Segments

        if (growthRates[s] != 1.0f) { 
            geometricFactor = (1.0f - std::pow( growthRates[s], static_cast<floatType>(meshSegments[s].nCells) )) / (1.0f - growthRates[s]);   // geometric series formula
        } else {
            geometricFactor = static_cast<floatType>( meshSegments[s].nCells );
        }
        segmentLength = meshSegments[s].endCoordinate - meshSegments[s].startCoordinate;
        firstCellLength = segmentLength / geometricFactor; 

        for (int i = 0; i != meshSegments[s].nCells; i++) {        // Cells within segment
            cellLengths( cellIndex ) = firstCellLength*std::pow( growthRates[s], static_cast<floatType>(i) );
            cellIndex++;
        }
    }

}



void CalculateCellCenters(Tensor1D &cellCenters, 
                          const Tensor1D &cellLengths,
                          const floatType startPosition)
{
    intType nCellsTotal = cellLengths.size();

    floatType previousCellPosition = startPosition, // just for the first iteration, this is actually cell face position
              previousCellLength = 0.0f;            // Needs to be zero so the above is true in the below formula
    for (intType i = 0; i != nCellsTotal; i++) {
        cellCenters(i) = previousCellPosition + previousCellLength/2.0f + cellLengths(i)/2.0f;
        previousCellPosition = cellCenters(i);
        previousCellLength = cellLengths(i);
    }

}



void CalculateCellCenterDiffInv(Tensor1D &cellCenterDiffInv, 
                                const Tensor1D &cellCenters)
{
    // First and last element dont correspond to valid values
    intType nFaces = cellCenters.size() + 1;

    for (intType i = 1; i != nFaces-1; i++) {
        cellCenterDiffInv(i) = 1.0f/( cellCenters(i) - cellCenters(i-1) );
    }
}



void CalculateCellFaces(Tensor1D &cellFaces, 
                        const Tensor1D &cellLengths,
                        const floatType startPosition)
{
    intType nFaces = cellFaces.size();

    cellFaces(0) = startPosition;
    for (intType i = 1; i != nFaces; i++) {
        cellFaces(i) = cellFaces(i-1) + cellLengths(i-1);
    }
}



void CalculateCellFaceAreas(Tensor2D &cellFaceAreas, 
                            const Tensor1D &cellLengths_x, 
                            const Tensor1D &cellLengths_y)
{
    for (int j = 0; j != cellLengths_y.dimension(0); j++) {
        for (int i = 0; i != cellLengths_x.dimension(0); i++) {
            cellFaceAreas(i, j) = cellLengths_x(i) * cellLengths_y(j);
        }
    }
}



void CalculateInterpolationFactors_WeightedLinear( Tensor1D &interpFactors, 
                                                   const Tensor1D &cellCenters, 
                                                   const Tensor1D &cellFaces) 
{
    for (int i = 1; i != interpFactors.size()-1; i++) {
        interpFactors(i) = ( cellFaces(i) - cellCenters(i-1) ) 
                         / ( cellCenters(i) - cellCenters(i-1) );
    }
}



void CalculateInterpolationFactors_Average( Tensor1D &interpFactors ) 
{
    for (int i = 1; i != interpFactors.size()-1; i++) {
        interpFactors(i) = 0.5f;
    }
}



Mesh::ExtrapFactorsStruct GetExtrapolationFactors(const Tensor1D &cellLengths, 
                                                  const intType fieldIndex_p, 
                                                  const intType fieldIndex_a)
{
    floatType extrapFactor_p = ( 2.0f*cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) )
                                / ( cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) );

    floatType extrapFactor_a = - ( cellLengths(fieldIndex_p) )
                                / ( cellLengths(fieldIndex_p) + cellLengths(fieldIndex_a) );
    
    return Mesh::ExtrapFactorsStruct{ extrapFactor_p, extrapFactor_a };
}



void CalculateExtrapolationFactors(EnumVector<BoundaryPatches, Mesh::ExtrapFactorsStruct > &extrapFactors, 
                                   const EnumVector<Axis, Tensor1D> &cellLengths, 
                                   const Axis::ENUMDATA axis)
{  

    using enum BoundaryPatches::ENUMDATA;
    using enum Axis::ENUMDATA;

    BoundaryPatches::ENUMDATA patchPositive = LUT::PositivePatch[ axis ], 
                              patchNegative = LUT::NegativePatch[ axis ];

    // If mesh is only 1 cell think (such as in a 2D simulation), make the extrapolatino equal to the single cell
    if ( cellLengths[axis].size() == 1 ) {
        extrapFactors[patchPositive].a = 0.0f;
        extrapFactors[patchPositive].p = 1.0f;

        extrapFactors[patchNegative].a = 0.0f;
        extrapFactors[patchNegative].p = 1.0f;
        return;
    }


    
    intType fieldIndex_p, fieldIndex_a; // Boundary cell node and the adjacent one

    // Positive patch boundary
    fieldIndex_p = cellLengths[axis].size() - 1;
    fieldIndex_a = fieldIndex_p - 1;
    extrapFactors[patchPositive] = GetExtrapolationFactors(cellLengths[axis], fieldIndex_p, fieldIndex_a);

    // Negative patch boundary
    fieldIndex_p = 0;
    fieldIndex_a = fieldIndex_p + 1;
    extrapFactors[patchNegative] = GetExtrapolationFactors(cellLengths[axis], fieldIndex_p, fieldIndex_a);
}



intType TotalCells(const std::vector<InputData::MeshSegment> &meshSegments)
{
    intType totalCells = 0;
    for (auto segment : meshSegments) {
        totalCells += segment.nCells;
    }
    return totalCells;
}



iArray3 NumberOfFaces( const iArray3 &nCells,
                       Axis::ENUMDATA axis)
{
    iArray3 nFaces;
    EnumFor<Axis>( [&] (Axis::ENUMDATA a) {
        nFaces(a) = nCells(a);
    } );
    nFaces(axis) += 1;  // There is one more faces than cells in the normal direction
    return nFaces;
}



std::pair<floatType, floatType> MinMaxCellGrowthRatios( const Mesh &mesh,
                                                        const Axis::ENUMDATA axis )
{
    floatType minCellGrowthRatio = std::numeric_limits<floatType>::infinity(),
              maxCellGrowthRatio = 0.0f;

    for ( intType i = 1; i != mesh.nCells[axis]; i++ ) {

        floatType growthRatio = mesh.cellLengths[axis](i) / mesh.cellLengths[axis](i-1);
        if ( growthRatio < 1.0f ) {
            growthRatio = 1.0f / growthRatio;
        }

        if ( growthRatio <= minCellGrowthRatio ) {
            minCellGrowthRatio = growthRatio;
        }

        if ( growthRatio >= maxCellGrowthRatio ) {
            maxCellGrowthRatio = growthRatio;
        }
    }

    return {minCellGrowthRatio, maxCellGrowthRatio};
}



std::pair<floatType, floatType> MinMaxCellAspectRatio( const Mesh &mesh )
{
    using enum Axis::ENUMDATA;

    floatType minAspectRatio = std::numeric_limits<floatType>::infinity(),
              maxAspectRatio = 0.0f;

    for ( intType k = 0; k != mesh.nCells[Z]; k++ ) {
        for ( intType j = 0; j != mesh.nCells[Y]; j++ ) {
            for ( intType i = 0; i != mesh.nCells[X]; i++ ) {

                floatType maxCellLength = std::max( {mesh.cellLengths[X](i), mesh.cellLengths[Y](j), mesh.cellLengths[Z](k)} );
                floatType minCellLength = std::min( {mesh.cellLengths[X](i), mesh.cellLengths[Y](j), mesh.cellLengths[Z](k)} );

                floatType aspectRatio = maxCellLength / minCellLength;

                if ( aspectRatio <= minAspectRatio ) {
                    minAspectRatio = aspectRatio;
                }

                if ( aspectRatio >= maxAspectRatio ) {
                    maxAspectRatio = aspectRatio;
                }

            }
        }
    }

    return {minAspectRatio, maxAspectRatio};
}



void SetCoarsenedCellLengths( Tensor1D &coarseCellLengths, 
                              const Tensor1D &fineCellLengths )
{

    intType startIndex, endIndexFine;
    if ( fineCellLengths.size() % 2 == 0 ) {    // Even number of cells

        startIndex = 0;
        endIndexFine = fineCellLengths.size();

    } else {                                    // Odd number of cells

        // Dont agglomorate the cell on the end that is the largest
        if ( fineCellLengths(0) >= fineCellLengths( fineCellLengths.size()-1 ) ) {
            startIndex = 1;
            endIndexFine = fineCellLengths.size();
            coarseCellLengths(0) = fineCellLengths(0);
        } else {
            startIndex = 0;
            endIndexFine = fineCellLengths.size() - 1;
            coarseCellLengths( coarseCellLengths.size()-1 ) = fineCellLengths( fineCellLengths.size()-1 );
        }

    }

    for ( intType iCoarse = startIndex, iFine = startIndex;  iFine != endIndexFine;  iCoarse++, iFine += 2 ) {
        coarseCellLengths(iCoarse) = fineCellLengths(iFine) + fineCellLengths(iFine + 1);
    }

}



}   // end anonymous namespace




// Constructor, allocates mesh given dimensions
Mesh::Mesh( const iArray3 &nCellsArg ) :
    nCells( nCellsArg ),

    nFacesNormal( { NumberOfFaces( nCells, Axis::X ), NumberOfFaces( nCells, Axis::Y ), NumberOfFaces( nCells, Axis::Z ) } ),

    cellCenters( {{Axis::ENUMDATA::X, nCells(0)},
                  {Axis::ENUMDATA::Y, nCells(1)},
                  {Axis::ENUMDATA::Z, nCells(2)}} ),

    cellFaces( {{Axis::ENUMDATA::X, nCells(0) + 1},
                {Axis::ENUMDATA::Y, nCells(1) + 1},
                {Axis::ENUMDATA::Z, nCells(2) + 1}} ),

    cellLengths( {{Axis::ENUMDATA::X, nCells(0)},
                  {Axis::ENUMDATA::Y, nCells(1)},
                  {Axis::ENUMDATA::Z, nCells(2)}} ),

    cellLengthsInv( {{Axis::ENUMDATA::X, nCells(0)},
                     {Axis::ENUMDATA::Y, nCells(1)},
                     {Axis::ENUMDATA::Z, nCells(2)}} ),

    cellCenterDiffInv( {{Axis::ENUMDATA::X, nCells(0) + 1},
                        {Axis::ENUMDATA::Y, nCells(1) + 1},
                        {Axis::ENUMDATA::Z, nCells(2) + 1}} ),

    interpFactors( {{Axis::ENUMDATA::X, nCells(0) + 1},
                    {Axis::ENUMDATA::Y, nCells(1) + 1},
                    {Axis::ENUMDATA::Z, nCells(2) + 1}} ),

    cellFaceAreas( {{Axis::ENUMDATA::X, {nCells(1), nCells(2)} },
                    {Axis::ENUMDATA::Y, {nCells(0), nCells(2)} },
                    {Axis::ENUMDATA::Z, {nCells(0), nCells(1)} }} ),

    extrapFactors()

    {};




// Constructor, creates the mesh from user inputdata
Mesh::Mesh(const InputData &inputData) :
    Mesh( { TotalCells(inputData.meshSegments[Axis::X]),  TotalCells(inputData.meshSegments[Axis::Y]), TotalCells(inputData.meshSegments[Axis::Z])} )

    { 
        std::vector< std::vector<floatType> > growthRates(Axis::count);

        EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

            growthRates[axis] = CalculateGrowthRates(inputData.meshSegments[axis]);

            CalculateCellLengths(cellLengths[axis], inputData.meshSegments[axis], growthRates[axis]);
            cellLengthsInv[axis] = cellLengths[axis].inverse();

            CalculateCellCenters(cellCenters[axis], cellLengths[axis], inputData.meshSegments[axis].front().startCoordinate);
            CalculateCellCenterDiffInv(cellCenterDiffInv[axis], cellCenters[axis]);

            CalculateCellFaces(cellFaces[axis], cellLengths[axis], inputData.meshSegments[axis].front().startCoordinate);

            switch ( inputData.schemes.faceInterpolationScheme ) {
                case FaceInterpolationSchemes::Average:
                    CalculateInterpolationFactors_Average(interpFactors[axis]);
                    break;

                case FaceInterpolationSchemes::WeightedLinear:
                    CalculateInterpolationFactors_WeightedLinear(interpFactors[axis], cellCenters[axis], cellFaces[axis]);
                    break;
            }
                        
            CalculateExtrapolationFactors(extrapFactors, cellLengths, axis);

        } );


        // Cell face areas should be calculated on their own since they depend on other axis
        EnumFor<Axis> ( [&] (Axis::ENUMDATA axis) {

            // Axis are ordered by numbering
            Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[ axis ];
            Axis::ENUMDATA axis2 = LUT::HiOrthogonalAxis[ axis ];
            CalculateCellFaceAreas(cellFaceAreas[axis], cellLengths[axis1], cellLengths[axis2]);

        } );

    };




Mesh CreateMesh( const InputData &inputData )
{
    using enum Axis::ENUMDATA;

    std::cout << "Generating mesh ... ";
    Mesh mesh(inputData);
    std::cout << "Success."
              << "\n";

    // Calculate the total number of cells
    intType nCells = mesh.nCells[X] * mesh.nCells[Y] * mesh.nCells [Z];

    // Maximum and minimum cell growth ratio in each coordinate direction
    EnumVector<Axis, std::pair<floatType, floatType>> minMaxCellGrowthRatios;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        minMaxCellGrowthRatios = MinMaxCellGrowthRatios( mesh, axis );
    } );

    // Maximum and minimum cell aspect ratio 
    auto [minAspectRatio, maxAspectRatio] = MinMaxCellAspectRatio( mesh );

    // Output this information to console
    std::cout << "Number of cells        : " << nCells << "\n"

              << "Min. cell growth ratios: " << "(" << minMaxCellGrowthRatios[X].first << ", " 
                                                       << minMaxCellGrowthRatios[Y].first << ", "
                                                       << minMaxCellGrowthRatios[Z].first << ")" << "\n" 

              << "Max. cell growth ratios: " << "(" << minMaxCellGrowthRatios[X].second << ", " 
                                                       << minMaxCellGrowthRatios[Y].second << ", "
                                                       << minMaxCellGrowthRatios[Z].second << ")" << "\n"
                
              << "Min. cell aspect ratio : " << minAspectRatio << "\n"
              
              << "Max. cell aspect ratio : " << maxAspectRatio << "\n\n";


    return mesh;
}



bool MeshCanBeCoarsened( const Mesh& mesh ) 
{
    for ( intType a = 0; a != Axis::count; a++ ) {
        Axis::ENUMDATA axis = static_cast<Axis::ENUMDATA>( a );
        if ( mesh.nCells(axis) > 2 )        // Need at least just one dimension to be coarsenable 
            return true;
    }
    return false;
}



Mesh CoarsenMesh( const Mesh &fineMesh,
                  const FaceInterpolationSchemes faceInterpolationScheme ) 
{

    // Determine dimensions of coarse mesh
    iArray3 coarseMeshDims;
    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {
        if ( fineMesh.nCells(axis) % 2 == 0 ) {
            coarseMeshDims(axis) = fineMesh.nCells(axis) / 2;
        } else {
            coarseMeshDims(axis) = 1  +  ( fineMesh.nCells(axis) - 1 ) / 2;
        }
    } );

    Mesh coarseMesh( coarseMeshDims );

    EnumFor<Axis>( [&] (Axis::ENUMDATA axis) {

        SetCoarsenedCellLengths( coarseMesh.cellLengths[axis], fineMesh.cellLengths[axis] );
        coarseMesh.cellLengthsInv[axis] = coarseMesh.cellLengths[axis].inverse();

        CalculateCellCenters(coarseMesh.cellCenters[axis], coarseMesh.cellLengths[axis], fineMesh.cellFaces[axis](0));
        CalculateCellCenterDiffInv(coarseMesh.cellCenterDiffInv[axis], coarseMesh.cellCenters[axis]);

        CalculateCellFaces(coarseMesh.cellFaces[axis], coarseMesh.cellLengths[axis], fineMesh.cellFaces[axis](0));

        switch ( faceInterpolationScheme ) {
            case FaceInterpolationSchemes::Average:
                CalculateInterpolationFactors_Average(coarseMesh.interpFactors[axis]);
                break;

            case FaceInterpolationSchemes::WeightedLinear:
                CalculateInterpolationFactors_WeightedLinear(coarseMesh.interpFactors[axis], coarseMesh.cellCenters[axis], coarseMesh.cellFaces[axis]);
                break;
        }
                    
        CalculateExtrapolationFactors(coarseMesh.extrapFactors, coarseMesh.cellLengths, axis);

    } );

    // Cell face areas should be calculated on their own since they depend on other axis
    EnumFor<Axis> ( [&] (Axis::ENUMDATA axis) {

        // Axis are ordered by numbering
        Axis::ENUMDATA axis1 = LUT::LoOrthogonalAxis[ axis ];
        Axis::ENUMDATA axis2 = LUT::HiOrthogonalAxis[ axis ];
        CalculateCellFaceAreas(coarseMesh.cellFaceAreas[axis], coarseMesh.cellLengths[axis1], coarseMesh.cellLengths[axis2]);

    } );
    

    return coarseMesh;
}


}   // end namespace CFD