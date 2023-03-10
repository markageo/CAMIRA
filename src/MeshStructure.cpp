#include "MeshStructure.h"
#include <iostream>

#include <cmath>
#include <functional>

// Helper functions
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


    std::vector<CFD::floatType> CalculateCellLengths(const std::vector<InputData::MeshSegment> &meshSegments, const std::vector<CFD::floatType> &growthRates)
    {
        int nSegments = meshSegments.size();
        std::vector<CFD::floatType> cellLengths;

        CFD::floatType segmentLength, firstCellLength, geometricFactor;
        for (int s = 0; s != nSegments; s++) {    // Segments

            if (growthRates[s] != 1.0) {
                geometricFactor = (1.0 - std::pow( growthRates[s], meshSegments[s].nCells )) / (1.0 - growthRates[s]);   // geometric series formula
            } else {
                geometricFactor = meshSegments[s].nCells;
            }
            segmentLength = meshSegments[s].upperBound - meshSegments[s].lowerBound;
            firstCellLength = segmentLength/geometricFactor; 

            for (int i = 0; i != meshSegments[s].nCells; i++) {        // Cells within segment
                cellLengths.push_back( firstCellLength*std::pow( growthRates[s], static_cast<CFD::floatType>(i) ) );
            }
        }

        return cellLengths;
    }


    void CalculateCellCenters(Eigen::Tensor<CFD::floatType, 1> &cellCenters, const std::vector<CFD::floatType> &cellLengths)
    {
        int nCellsTotal = cellLengths.size();

        CFD::floatType previousCellPosition = 0.0, previousCellLength = 0.0;
        for (int i = 0; i != nCellsTotal; i++) {
            cellCenters(i) = previousCellPosition + previousCellLength/2.0 + cellLengths[i]/2.0;
            previousCellPosition = cellCenters(i);
            previousCellLength = cellLengths[i];
        }

    }

    void CalculateCellFaceAreas(Eigen::Tensor<CFD::floatType, 2> &cellFaceAreas, const std::vector<CFD::floatType> &cellLengths_x, 
        const std::vector<CFD::floatType> &cellLengths_y)
    {
        using v_size_type = std::vector<CFD::floatType>::size_type;

        for (v_size_type i = 0; i != cellLengths_x.size(); i++) {
            for (v_size_type j = 0; j != cellLengths_y.size(); j++) {
                cellFaceAreas(i, j) = cellLengths_x[i] * cellLengths_y[j];
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
CFD::MeshStructure::MeshStructure(const InputData &inputData) :
    cellCenters_x( TotalCells(inputData.meshSegments_x) ),
    cellCenters_y( TotalCells(inputData.meshSegments_y) ),
    cellCenters_z( TotalCells(inputData.meshSegments_z) ),
    cellFaceAreas_x( cellCenters_y.dimension(0), cellCenters_z.dimension(0) ),
    cellFaceAreas_y( cellCenters_z.dimension(0), cellCenters_x.dimension(0) ),
    cellFaceAreas_z( cellCenters_x.dimension(0), cellCenters_y.dimension(0) )
    { 
        std::vector<CFD::floatType> growthRates_x = CalculateGrowthRates(inputData.meshSegments_x);
        std::vector<CFD::floatType> growthRates_y = CalculateGrowthRates(inputData.meshSegments_y);
        std::vector<CFD::floatType> growthRates_z = CalculateGrowthRates(inputData.meshSegments_z);

        std::vector<CFD::floatType> cellLengths_x = CalculateCellLengths(inputData.meshSegments_x, growthRates_x);
        std::vector<CFD::floatType> cellLengths_y = CalculateCellLengths(inputData.meshSegments_y, growthRates_y);
        std::vector<CFD::floatType> cellLengths_z = CalculateCellLengths(inputData.meshSegments_z, growthRates_z);

        CalculateCellCenters(cellCenters_x, cellLengths_x);
        CalculateCellCenters(cellCenters_y, cellLengths_y);
        CalculateCellCenters(cellCenters_z, cellLengths_z);

        CalculateCellFaceAreas(cellFaceAreas_x, cellLengths_y, cellLengths_z);
        CalculateCellFaceAreas(cellFaceAreas_y, cellLengths_z, cellLengths_x);
        CalculateCellFaceAreas(cellFaceAreas_z, cellLengths_x, cellLengths_y);
    };


