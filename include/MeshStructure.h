#ifndef MESH_STRUCTURE
#define MESH_STRUCTURE

#include "SimulationParameters.h"
#include "InputProcessing.h"
#include "Tensor"


namespace CFD
{

// Rectilinear mesh structure
class MeshStructure
{

    public:

        MeshStructure(const CFD::InputData &);

        // TODO: These should maybe be encapsulated in getters so that they cannot be changed from the outside
        Eigen::Tensor<CFD::floatType, 1> cellCenters_x, cellCenters_y, cellCenters_z;
        Eigen::Tensor<CFD::floatType, 2> cellFaceAreas_x, cellFaceAreas_y, cellFaceAreas_z; // Index by right hand rule

    private:

};

}   // end namespace CFD

#endif  // MESH_STRCTURE