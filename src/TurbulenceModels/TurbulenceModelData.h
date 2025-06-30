#ifndef CAMIRA_TURBULENCE_MODEL_DATA
#define CAMIRA_TURBULENCE_MODEL_DATA

#include "TurbulenceModels.h"
#include "../Core/Types.h"
#include "../IO/InputProcessing.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/BoundaryConditionData.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"

#include <memory>

namespace CAMIRA
{


struct TurbulenceModelData
{
    TurbulenceModels model;
    std::unique_ptr<TurbulenceModelInterface> turbModel;
};

TurbulenceModelData CreateTurbulenceModelData( const InputData &, 
                                               const Mesh &, 
                                               const IBData &,
                                               const BoundaryConditionData & );


}   // end namespace CAMIRA


#endif // CAMIRA_TURBULENCE_MODEL_DATA