#ifndef CFD_MULTIGRID
#define CFD_MULTIGRID

#include "../Types.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../IO/InputProcessing.h"


namespace CFD
{

struct GridLevelData
{
    intType level;
    Mesh mesh;
    FieldData<Tensor3D> fineGridApproximationFields,
                        coarseGridFields,
                        residual;
    EnumVector<Axis, Tensor3D> faceFluxes;
    EnumVector< Axis, EnumVector< Axis, Tensor3D> > faceAdvectedVelocities;
    IBData ibData;
    FVCoefficients fvCoeffs;
};



// Create initial heirachy of grids
std::vector<GridLevelData> CreateMGLevels( const InputData & );


}

#endif  // CFD_MULTIGRID