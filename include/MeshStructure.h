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

        array1D cellCenters_x, cellCenters_y, cellCenters_z;
        array2D cellFaceAreas_x, cellFaceAreas_y, cellFaceAreas_z; // Index by right hand rule

    private:

};

}   // end namespace CFD

#endif  // MESH_STRCTURE