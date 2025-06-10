#include "TurbulenceModels.h"

#include <memory>

namespace CAMIRA
{

TurbulenceModelData::TurbulenceModelData(const TurbulenceModels turbulenceModel, 
                                         const Mesh &mesh, 
                                         const Polyhedron &geometry, 
                                         const BoundaryConditionData &bcData ) :
    model( turbulenceModel )
{
    switch ( model ) {
        case TurbulenceModels::Laminar:
            turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::Laminar> >(); 
            break;
        case TurbulenceModels::ChenAndXuZeroEquation:
            turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::ChenAndXuZeroEquation> >(); 
            break;
        case TurbulenceModels::PrandtlZeroEquation:
            turbModel = std::make_unique< TurbulenceModel<TurbulenceModels::PrandtlZeroEquation> >();
            break;
    }
    turbModel->SetTurbulenceModelData( mesh, geometry, bcData );
}


}   // end namespace CAMIRA