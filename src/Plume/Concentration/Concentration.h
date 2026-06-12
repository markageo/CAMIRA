#ifndef CAMIRA_PLUME_CONCENTRATION
#define CAMIRA_PLUME_CONCENTRATION

#include "Core/Types.h"
#include "Core/Mesh/Mesh.h"
#include "Plume/Particles/Particles.h"
#include "Plume/ConfigEnums.h"


#include <vector>
#include <utility>
#include <map>
#include <tuple>
#include <memory>

namespace CAMIRA
{

using namespace CORE;

namespace PLUME
{

void UpdateConcentrationField( Tensor3D &,
                               const Particles &,
                               const Mesh & );


}   // end namespace PLUME

}   // end namespace CAMIRA

#endif  // CAMIRA_PLUME_CONCENTRATION