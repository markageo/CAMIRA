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

        // TODO: These should maybe be encapsulated in getters so that they cannot be changed from the outside
        Eigen::Tensor<SIM::floatType, 1> cellCenters_x, cellCenters_y, cellCenters_z;
        Eigen::Tensor<SIM::floatType, 2> cellFaceAreas_x, cellFaceAreas_y, cellFaceAreas_z; // Index by right hand rule

    private:

};

#endif  // MESH_STRCTURE