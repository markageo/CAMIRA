#include "TurbulenceModels.h"

namespace CAMIRA
{

void TurbulenceModel<TurbulenceModels::Laminar>::SetTurbulenceModelData( [[ maybe_unused ]] const Mesh &mesh,
                                                                         [[ maybe_unused ]] const Polyhedron &geometry,
                                                                         [[ maybe_unused ]] const BoundaryConditionData &bcData )
{ /* NULL */ }


void TurbulenceModel<TurbulenceModels::Laminar>::SetTurbulenceViscosityField( EnumVector<Axis, Tensor3D> &nuTurbulent,
                                                                              [[ maybe_unused ]] const FieldData<Tensor3D> &fields,
                                                                              [[ maybe_unused ]] const Mesh &mesh )
{
    EnumFor<Axis>( [&] (Axis::ENUMDATA faceNormal) {
        nuTurbulent[faceNormal].setZero();
    } );
}

}   // end namespace CAMIRA