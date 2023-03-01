#include "MeshStructure.h"


MeshStructure::MeshStructure(const InputData &inputData) :
    cellCenters_x(  std::reduce(inputData.mesh.nCells_x.begin(), inputData.mesh.nCells_x.end())  ),
    cellCenters_y(  std::reduce(inputData.mesh.nCells_y.begin(), inputData.mesh.nCells_y.end())  ),
    cellCenters_z(  std::reduce(inputData.mesh.nCells_z.begin(), inputData.mesh.nCells_z.end())  )
    { CreateMesh(inputData); };



void MeshStructure::CreateMesh(const InputData &)
{
    
}