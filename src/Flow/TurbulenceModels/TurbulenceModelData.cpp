#include "TurbulenceModelData.h"
#include "TurbulenceModels.h"
#include "Core/Geometry/Geometry.h"

#include <memory>

namespace CAMIRA
{

void SetTurbulenceModelData( TurbulenceModelData &turbModelData, 
                             const InputData &inputData, 
                             const AxisTransformationMap &axisTransformation,
                             const Mesh &mesh,  
                             const IBData &ibData,
                             const BoundaryConditionData &bcData )
{
    turbModelData.model = inputData.turbulenceModel;

    switch ( turbModelData.model ) {
        case TurbulenceModels::Laminar:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::Laminar> >(); 
            break;
        case TurbulenceModels::PrandtlZeroEquation:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::PrandtlZeroEquation> >();
            break;
        case TurbulenceModels::ZEQ0:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::ZEQ0> >(); 
            break;
        case TurbulenceModels::ZEQ1:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::ZEQ1> >(); 
            break;
        case TurbulenceModels::ZEQ2:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::ZEQ2> >(); 
            break;
        case TurbulenceModels::ZEQ3:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::ZEQ3> >(); 
            break;
        case TurbulenceModels::ZEQ4:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::ZEQ4> >(); 
            break;
        case TurbulenceModels::Null:
            turbModelData.turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::Laminar> >();     // Laminar just sets tubulence viscosity to zero
            break;
    }
    turbModelData.turbModel->SetTurbulenceModelData( inputData, axisTransformation, mesh, ibData, bcData );
}


}   // end namespace CAMIRA