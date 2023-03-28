#include "Types.h"
#include "FiniteVolumeStructures.h"
#include "Tensor"

#include <iostream>
#include <cmath>
#include <functional>


/*-------------------------------------------------------------------------------------*\
                                    Array Allocator
\*-------------------------------------------------------------------------------------*/

/* NULL */


/*-------------------------------------------------------------------------------------*\
                                         Mesh
\*-------------------------------------------------------------------------------------*/

namespace
{

    using namespace CFD;

    std::vector<CFD::floatType> CalculateGrowthRates(const std::vector<InputData::MeshSegment> &meshSegments)
    {
        int nSegments = meshSegments.size();
        std::vector<CFD::floatType> growthRates(nSegments);

        // Negative growth rate means shrinking grid
        for (int i = 0; i != nSegments; i++) {
            growthRates[i] = std::pow( std::abs(meshSegments[i].biasFactor) , 1.0/( meshSegments[i].nCells - 1 ) );
            if (meshSegments[i].biasFactor < 0 )   
                growthRates[i] = 1.0/growthRates[i];
        }

        return growthRates;
    }


    void CalculateCellLengths(array1D &cellLengths, const std::vector<InputData::MeshSegment> &meshSegments, const std::vector<CFD::floatType> &growthRates)
    {
        int nSegments = meshSegments.size();

        CFD::floatType segmentLength, firstCellLength, geometricFactor;
        int cellIndex = 0;
        for (int s = 0; s != nSegments; s++) {    // Segments

            if (growthRates[s] != 1.0) {
                geometricFactor = (1.0 - std::pow( growthRates[s], meshSegments[s].nCells )) / (1.0 - growthRates[s]);   // geometric series formula
            } else {
                geometricFactor = meshSegments[s].nCells;
            }
            segmentLength = meshSegments[s].upperBound - meshSegments[s].lowerBound;
            firstCellLength = segmentLength/geometricFactor; 

            for (int i = 0; i != meshSegments[s].nCells; i++) {        // Cells within segment
                cellLengths( cellIndex ) = firstCellLength*std::pow( growthRates[s], static_cast<CFD::floatType>(i) );
                cellIndex++;
            }
        }

    }


    void CalculateCellCenters(CFD::array1D &cellCenters, const array1D &cellLengths)
    {
        int nCellsTotal = cellLengths.size();

        CFD::floatType previousCellPosition = 0.0, previousCellLength = 0.0;
        for (int i = 0; i != nCellsTotal; i++) {
            cellCenters(i) = previousCellPosition + previousCellLength/2.0 + cellLengths(i)/2.0;
            previousCellPosition = cellCenters(i);
            previousCellLength = cellLengths(i);
        }

    }


    void CalculateCellFaces(CFD::array1D &cellFaces, const array1D &cellLengths)
    {
        int nFaces = cellLengths.size() + 1;

        CFD::floatType previousFacePosition = 0.0, previousCellLength = 0.0;
        for (int i = 0; i != nFaces; i++) {
            cellFaces(i) = previousFacePosition + previousCellLength;
            previousFacePosition = cellFaces(i);
            previousCellLength = cellLengths(i);
        }

    }


    void CalculateInterpolationFactors(CFD::array1D &interpFactors, const CFD::array1D &cellCenters, const CFD::array1D &cellFaces) 
    {
        for (int i = 1; i != interpFactors.dimension(0)-1; i++) {
            interpFactors(i) = ( cellFaces(i) - cellCenters(i-1) ) / ( cellCenters(i) - cellCenters(i-1) );
        }
    }


    void CalculateCellFaceAreas(CFD::array2D &cellFaceAreas, const array1D &cellLengths_x, const array1D &cellLengths_y)
    {
        for (int j = 0; j != cellLengths_y.dimension(0); j++) {
            for (int i = 0; i != cellLengths_x.dimension(0); i++) {
                cellFaceAreas(i, j) = cellLengths_x(i) * cellLengths_y(j);
            }
        }
    }


    CFD::intType TotalCells(const std::vector<InputData::MeshSegment> &meshSegments)
    {
        CFD::intType totalCells = 0;
        for (auto segment : meshSegments) {
            totalCells += segment.nCells;
        }
        return totalCells;
    }

}

// Constructor, creates the mesh
CFD::Mesh::Mesh(const InputData &inputData) :
    nCells( { TotalCells(inputData.meshSegments_x),  TotalCells(inputData.meshSegments_y), TotalCells(inputData.meshSegments_z)} ),

    cellCenters( {{Axis::ENUMDATA::X, nCells(0)},
                  {Axis::ENUMDATA::Y, nCells(1)},
                  {Axis::ENUMDATA::Z, nCells(2)}} ),

    cellFaces( {{Axis::ENUMDATA::X, nCells(0) + 1},
                {Axis::ENUMDATA::Y, nCells(1) + 1},
                {Axis::ENUMDATA::Z, nCells(2) + 1}} ),

    cellLengths( {{Axis::ENUMDATA::X, nCells(0)},
                  {Axis::ENUMDATA::Y, nCells(1)},
                  {Axis::ENUMDATA::Z, nCells(2)}} ),

    interpFactors( {{Axis::ENUMDATA::X, nCells(0) + 1},
                    {Axis::ENUMDATA::Y, nCells(1) + 1},
                    {Axis::ENUMDATA::Z, nCells(2) + 1}} ),

    cellFaceAreas( {{Axis::ENUMDATA::X, {nCells(1), nCells(2)} },
                    {Axis::ENUMDATA::Y, {nCells(2), nCells(0)} },
                    {Axis::ENUMDATA::Z, {nCells(0), nCells(1)} }} )

    { 
        std::vector<CFD::floatType> growthRates_x = CalculateGrowthRates(inputData.meshSegments_x);
        std::vector<CFD::floatType> growthRates_y = CalculateGrowthRates(inputData.meshSegments_y);
        std::vector<CFD::floatType> growthRates_z = CalculateGrowthRates(inputData.meshSegments_z);

        using enum Axis::ENUMDATA;

        CalculateCellLengths(cellLengths[X], inputData.meshSegments_x, growthRates_x);
        CalculateCellLengths(cellLengths[Y], inputData.meshSegments_y, growthRates_y);
        CalculateCellLengths(cellLengths[Z], inputData.meshSegments_z, growthRates_z);

        CalculateCellCenters(cellCenters[X], cellLengths[X]);
        CalculateCellCenters(cellCenters[Y], cellLengths[Y]);
        CalculateCellCenters(cellCenters[Z], cellLengths[Z]);

        CalculateCellFaces(cellFaces[X], cellLengths[X]);
        CalculateCellFaces(cellFaces[Y], cellLengths[Y]);
        CalculateCellFaces(cellFaces[Z], cellLengths[Z]);

        CalculateInterpolationFactors(interpFactors[X], cellCenters[X], cellFaces[X]);
        CalculateInterpolationFactors(interpFactors[Y], cellCenters[Y], cellFaces[Y]);
        CalculateInterpolationFactors(interpFactors[Z], cellCenters[Z], cellFaces[Z]);

        CalculateCellFaceAreas(cellFaceAreas[X], cellLengths[Y], cellLengths[Z]);
        CalculateCellFaceAreas(cellFaceAreas[Y], cellLengths[Z], cellLengths[X]);
        CalculateCellFaceAreas(cellFaceAreas[Z], cellLengths[X], cellLengths[Y]);
    };
