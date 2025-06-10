#ifndef CAMIRA_TURBULENCE_MODELS
#define CAMIRA_TURBULENCE_MODELS

#include "../Core/Types.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/BoundaryConditionData.h"
#include "../Geometry/Geometry.h"

#include <memory>

namespace CAMIRA
{


struct TurbulenceModelInterface
{
    virtual void SetTurbulenceModelData(const Mesh &, const Polyhedron &, const BoundaryConditionData &) = 0; 
    virtual void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const Mesh &) = 0;
};


template< TurbulenceModels T >
struct TurbulenceModel;


template<>
struct TurbulenceModel< TurbulenceModels::Laminar > : public TurbulenceModelInterface
{
    void SetTurbulenceModelData(const Mesh &, const Polyhedron &, const BoundaryConditionData &);
    void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const Mesh &);
};


template<>
struct TurbulenceModel< TurbulenceModels::ChenAndXuZeroEquation > : public TurbulenceModelInterface
{
    void SetTurbulenceModelData(const Mesh &, const Polyhedron &, const BoundaryConditionData &);
    void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const Mesh &);

    private:
        floatType m_proportionalityConstant;
        EnumVector<Axis, Tensor3D> m_wallDistance; // Stored at cell faces
};


template<>
struct TurbulenceModel< TurbulenceModels::PrandtlZeroEquation > : public TurbulenceModelInterface
{
    void SetTurbulenceModelData(const Mesh &, const Polyhedron &, const BoundaryConditionData &);
    void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const Mesh &);

    private:
        floatType m_vonKarmanConstant;
        EnumVector<Axis, Tensor3D> m_wallDistance; // Stored at cell faces
        Tensor3D m_velocityDeformationRate; // Stored at cell centers, to be interpolated to faces
};


}   // end namespace CAMIRA


#endif // CAMIRA_TURBULENCE_MODELS