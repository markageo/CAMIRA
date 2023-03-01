#ifndef MESH_STRUCTURE
#define MESH_STRUCTURE

#include "SimulationParameters.h"
#include "InputProcessing.h"
#include "Tensor"

// Rectilinear mesh structure
class MeshStructure
{

    public:

        MeshStructure(const InputData &);

        Eigen::Tensor<SIM::floatType, 1> cellCenters_x, cellCenters_y, cellCenters_z;

    private:

        
        void CreateMesh(const InputData &);

};

#endif  // MESH_STRCTURE