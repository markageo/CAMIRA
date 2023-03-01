#include "MeshStructure.h"
#include <iostream>

#include <cmath>


namespace
{
    std::vector<SIM::floatType> CalculateGrowthRates(const std::vector<SIM::floatType> &biasFactors, 
    const std::vector<SIM::intType> &nCells)
    {
        int nSegments = nCells.size();
        std::vector<SIM::floatType> growthRates(nSegments);

        // Negative growth rate means shrinking grid
        for (int i = 0; i != nSegments; i++) {
            growthRates[i] = std::pow( std::abs(biasFactors[i]) , 1.0/( nCells[i] - 1 ) );
            if (biasFactors[i] < 0 )   
                growthRates[i] = 1.0/growthRates[i];
        }

        return growthRates;
    }


    void CalculateCellCenters(Eigen::Tensor<SIM::floatType, 1> &cellCenters, const std::vector<SIM::floatType> &growthRates, const std::vector<SIM::intType> &nCells, 
            const std::vector<std::pair<SIM::floatType, SIM::floatType>> &segmentBounds)
    {
        int nSegments = nCells.size();

        SIM::floatType segmentLength, firstCellLength, geometricFactor;
        SIM::floatType currentCellWidth, previousCellWidth = 0.0, previousCellPosition = 0.0;
        SIM::intType gridIndex = 0;
        for (int s = 0; s != nSegments; s++) {    // Segments

            if (growthRates[s] != 1.0) {
                geometricFactor = (1.0 - std::pow( growthRates[s], nCells[s] )) / (1.0 - growthRates[s]);   // geometric series formaula
            } else {
                geometricFactor = nCells[s];
            }
            segmentLength = segmentBounds[s].second - segmentBounds[s].first;
            firstCellLength = segmentLength/geometricFactor; 

            for (int i = 0; i != nCells[s]; i++) {        // Cells within segment
                currentCellWidth = firstCellLength*std::pow( growthRates[s], static_cast<SIM::floatType>(i) );

                cellCenters(gridIndex) = previousCellPosition + previousCellWidth/2.0 + currentCellWidth/2.0;

                previousCellPosition = cellCenters(gridIndex);
                previousCellWidth = currentCellWidth;
                gridIndex++;
            }
        }
    }



}


MeshStructure::MeshStructure(const InputData &inputData) :
    cellCenters_x(  std::reduce(inputData.mesh.nCells_x.begin(), inputData.mesh.nCells_x.end())  ),
    cellCenters_y(  std::reduce(inputData.mesh.nCells_y.begin(), inputData.mesh.nCells_y.end())  ),
    cellCenters_z(  std::reduce(inputData.mesh.nCells_z.begin(), inputData.mesh.nCells_z.end())  )
    { CreateMesh(inputData);};



void MeshStructure::CreateMesh(const InputData &inputData)
{
    
    std::vector<SIM::floatType> growthRates_x = CalculateGrowthRates(inputData.mesh.biasFactors_x, inputData.mesh.nCells_x);
    std::vector<SIM::floatType> growthRates_y = CalculateGrowthRates(inputData.mesh.biasFactors_y, inputData.mesh.nCells_y);
    std::vector<SIM::floatType> growthRates_z = CalculateGrowthRates(inputData.mesh.biasFactors_z, inputData.mesh.nCells_z);

    CalculateCellCenters(cellCenters_x, growthRates_x, inputData.mesh.nCells_x, inputData.mesh.segmentBounds_x);
    CalculateCellCenters(cellCenters_y, growthRates_y, inputData.mesh.nCells_y, inputData.mesh.segmentBounds_y);
    CalculateCellCenters(cellCenters_z, growthRates_z, inputData.mesh.nCells_z, inputData.mesh.segmentBounds_z);

}

