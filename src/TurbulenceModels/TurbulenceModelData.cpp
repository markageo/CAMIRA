#include "TurbulenceModelData.h"
#include "TurbulenceModels.h"
#include "../Geometry/Geometry.h"

#include <memory>

namespace CAMIRA
{

TurbulenceModelData CreateTurbulenceModelData( const InputData &inputData, 
                                               const Mesh &mesh,  
                                               const BoundaryConditionData &bcData )
{
    TurbulenceModelData turbModelData;

    turbModelData.model = inputData.turbulenceModel;

    Polyhedron geometry = MakeGeometry( inputData );

    switch ( turbModelData.model ) {
        case TurbulenceModels::Laminar:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::Laminar> >(); 
            break;
        case TurbulenceModels::ChenAndXuZeroEquation:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::ChenAndXuZeroEquation> >(); 
            break;
        case TurbulenceModels::PrandtlZeroEquation:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::PrandtlZeroEquation> >();
            break;
        case TurbulenceModels::Null:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::Laminar> >();     // Laminar just sets tubulence viscosity to zero
            break;
    }
    turbModelData.turbModel->SetTurbulenceModelData( mesh, geometry, bcData );

    return turbModelData;
}


}   // end namespace CAMIRA