#ifndef CAMIRA_TURBULENCE_MODEL_DATA
#define CAMIRA_TURBULENCE_MODEL_DATA

#include "TurbulenceModels.h"
#include "Core/Types.h"
#include "Core/Mesh/Mesh.h"
#include "Flow/InputProcessing/InputProcessing.h"
#include "Flow/FiniteVolume/BoundaryConditionData.h"
#include "Flow/ImmersedBoundary/ImmersedBoundary.h"

#include <memory>

namespace CAMIRA
{


struct TurbulenceModelData
{
    TurbulenceModels model;
    std::unique_ptr<TurbulenceModelInterface> turbModel;
};

void SetTurbulenceModelData( TurbulenceModelData &,
                             const InputData &, 
                             const AxisTransformationMap &,
                             const Mesh &, 
                             const IBData &,
                             const BoundaryConditionData & );


}   // end namespace CAMIRA


#endif // CAMIRA_TURBULENCE_MODEL_DATA