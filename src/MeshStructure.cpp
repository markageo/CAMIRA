#include "MeshStructure.h"
#include <iostream>

#include <cmath>

MeshStructure::MeshStructure(const InputData &inputData) :
    cellCenters_x(  std::reduce(inputData.mesh.nCells_x.begin(), inputData.mesh.nCells_x.end())  ),
    cellCenters_y(  std::reduce(inputData.mesh.nCells_y.begin(), inputData.mesh.nCells_y.end())  ),
    cellCenters_z(  std::reduce(inputData.mesh.nCells_z.begin(), inputData.mesh.nCells_z.end())  )
    { CreateMesh(inputData);};



void MeshStructure::CreateMesh(const InputData &inputData)
{
    int nSegments_x = inputData.mesh.nCells_x.size();
    int nSegments_y = inputData.mesh.nCells_y.size();
    int nSegments_z = inputData.mesh.nCells_z.size();

    // Calculate growth rates in each segment
    std::vector<SIM::floatType> growthRates_x( inputData.mesh.biasFactors_x.size() ), 
                                growthRates_y( inputData.mesh.biasFactors_y.size() ), 
                                growthRates_z( inputData.mesh.biasFactors_z.size() );


    // Negative growth rate means shrinking grid
    for (int i = 0; i != nSegments_x; i++) {
        growthRates_x[i] = std::pow( std::abs(inputData.mesh.biasFactors_x[i]) , 1.0/( inputData.mesh.nCells_x[i] - 1 ) );
        if (inputData.mesh.biasFactors_x[i] < 0 )   
            growthRates_x[i] = 1.0/growthRates_x[i];
    }

    for (int i = 0; i != nSegments_y; i++) {
        growthRates_y[i] = std::pow( std::abs(inputData.mesh.biasFactors_y[i]) , 1.0/( inputData.mesh.nCells_y[i] - 1 ) );
        if (inputData.mesh.biasFactors_y[i] < 0)
            growthRates_y[i] = 1.0/growthRates_y[i];
    }
    
    for (int i = 0; i != nSegments_z; i++) {
        growthRates_z[i] = std::pow( std::abs(inputData.mesh.biasFactors_z[i]) , 1.0/( inputData.mesh.nCells_z[i] - 1 ) );
        if (inputData.mesh.biasFactors_z[i] < 0)
            growthRates_z[i] = 1.0/growthRates_z[i];
    }



    // x
    SIM::floatType segmentLength, firstCellLength, geometricFactor;
    SIM::floatType currentCellWidth, previousCellWidth = 0.0, previousCellPosition = 0.0;
    SIM::intType gridIndex = 0;
    for (int s = 0; s != nSegments_x; s++) {    // Segments


        if (growthRates_x[s] != 1) {
            geometricFactor = (1 - std::pow( growthRates_x[s], inputData.mesh.nCells_x[s] )) / (1 - growthRates_x[s]);   // geometric series
        } else {
            geometricFactor = inputData.mesh.nCells_x[s];
        }
        segmentLength = inputData.mesh.segmentBounds_x[s].second - inputData.mesh.segmentBounds_x[s].first;
        firstCellLength = segmentLength/geometricFactor; 

        for (int i = 0; i != inputData.mesh.nCells_x[s]; i++) {        // Cells within segment
            currentCellWidth = firstCellLength*std::pow( growthRates_x[s], static_cast<SIM::floatType>(i) );
            cellCenters_x(gridIndex)= previousCellPosition + previousCellWidth/2.0 + currentCellWidth/2.0;
            previousCellPosition = cellCenters_x(gridIndex);
            previousCellWidth = currentCellWidth;
            gridIndex++;
        }
    }


    // y
    previousCellWidth = 0.0;
    previousCellPosition = 0.0;
    gridIndex = 0;
    for (int s = 0; s != nSegments_y; s++) {    // Segments


        if (growthRates_y[s] != 1) {
            geometricFactor = (1 - std::pow( growthRates_y[s], inputData.mesh.nCells_y[s] )) / (1 - growthRates_y[s]);   // geometric series
        } else {
            geometricFactor = inputData.mesh.nCells_y[s];
        }
        segmentLength = inputData.mesh.segmentBounds_y[s].second - inputData.mesh.segmentBounds_y[s].first;
        firstCellLength = segmentLength/geometricFactor; 

        for (int i = 0; i != inputData.mesh.nCells_y[s]; i++) {        // Cells within segment
            currentCellWidth = firstCellLength*std::pow( growthRates_y[s], static_cast<SIM::floatType>(i) );
            cellCenters_y(gridIndex)= previousCellPosition + previousCellWidth/2.0 + currentCellWidth/2.0;
            previousCellPosition = cellCenters_y(gridIndex);
            previousCellWidth = currentCellWidth;
            gridIndex++;
        }
    }


    // z
    previousCellWidth = 0.0;
    previousCellPosition = 0.0;
    gridIndex = 0;
    for (int s = 0; s != nSegments_z; s++) {    // Segments

        if (growthRates_z[s] != 1) {
            geometricFactor = (1 - std::pow( growthRates_z[s], inputData.mesh.nCells_z[s] )) / (1 - growthRates_z[s]);   // geometric series
        } else {
            geometricFactor = inputData.mesh.nCells_z[s];
        }
        segmentLength = inputData.mesh.segmentBounds_z[s].second - inputData.mesh.segmentBounds_z[s].first;
        firstCellLength = segmentLength/geometricFactor; 

        for (int i = 0; i != inputData.mesh.nCells_z[s]; i++) {        // Cells within segment
            currentCellWidth = firstCellLength*std::pow( growthRates_z[s], static_cast<SIM::floatType>(i) );
            cellCenters_z(gridIndex) = previousCellPosition + previousCellWidth/2.0 + currentCellWidth/2.0;
            previousCellPosition = cellCenters_z(gridIndex);
            previousCellWidth = currentCellWidth;
            gridIndex++;
        }
    }

}