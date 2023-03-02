#include "MeshStructure.h"
#include <iostream>

#include <cmath>
#include <functional>

// Helper functions
namespace
{
    std::vector<SIM::floatType> CalculateGrowthRates(const std::vector<InputData::MeshSegment> &meshSegments)
    {
        int nSegments = meshSegments.size();
        std::vector<SIM::floatType> growthRates(nSegments);

        // Negative growth rate means shrinking grid
        for (int i = 0; i != nSegments; i++) {
            growthRates[i] = std::pow( std::abs(meshSegments[i].biasFactor) , 1.0/( meshSegments[i].nCells - 1 ) );
            if (meshSegments[i].biasFactor < 0 )   
                growthRates[i] = 1.0/growthRates[i];
        }

        return growthRates;
    }


    void CalculateCellCenters(Eigen::Tensor<SIM::floatType, 1> &cellCenters, const std::vector<SIM::floatType> &growthRates, 
        const std::vector<InputData::MeshSegment> &meshSegments)
    {
        int nSegments = meshSegments.size();

        SIM::floatType segmentLength, firstCellLength, geometricFactor;
        SIM::floatType currentCellWidth, previousCellWidth = 0.0, previousCellPosition = 0.0;
        SIM::intType gridIndex = 0;
        for (int s = 0; s != nSegments; s++) {    // Segments

            if (growthRates[s] != 1.0) {
                geometricFactor = (1.0 - std::pow( growthRates[s], meshSegments[s].nCells )) / (1.0 - growthRates[s]);   // geometric series formaula
            } else {
                geometricFactor = meshSegments[s].nCells;
            }
            segmentLength = meshSegments[s].upperBound - meshSegments[s].lowerBound;
            firstCellLength = segmentLength/geometricFactor; 

            for (int i = 0; i != meshSegments[s].nCells; i++) {        // Cells within segment
                currentCellWidth = firstCellLength*std::pow( growthRates[s], static_cast<SIM::floatType>(i) );

                cellCenters(gridIndex) = previousCellPosition + previousCellWidth/2.0 + currentCellWidth/2.0;

                previousCellPosition = cellCenters(gridIndex);
                previousCellWidth = currentCellWidth;
                gridIndex++;
            }
        }
    }


    SIM::intType TotalCells(const std::vector<InputData::MeshSegment> &meshSegments)
    {
        SIM::intType totalCells = 0;
        for (auto segment : meshSegments) {
            totalCells += segment.nCells;
        }
        return totalCells;
    }

}

// Constructor
MeshStructure::MeshStructure(const InputData &inputData) :
    cellCenters_x( TotalCells(inputData.meshSegments_x) ),
    cellCenters_y( TotalCells(inputData.meshSegments_y) ),
    cellCenters_z( TotalCells(inputData.meshSegments_z) )
    { 
        std::vector<SIM::floatType> growthRates_x = CalculateGrowthRates(inputData.meshSegments_x);
        std::vector<SIM::floatType> growthRates_y = CalculateGrowthRates(inputData.meshSegments_y);
        std::vector<SIM::floatType> growthRates_z = CalculateGrowthRates(inputData.meshSegments_z);

        CalculateCellCenters(cellCenters_x, growthRates_x, inputData.meshSegments_x);
        CalculateCellCenters(cellCenters_y, growthRates_y, inputData.meshSegments_y);
        CalculateCellCenters(cellCenters_z, growthRates_z, inputData.meshSegments_z);
    };


