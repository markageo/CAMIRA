#ifndef CAMIRA_TURBULENCE_MODELS
#define CAMIRA_TURBULENCE_MODELS

#include "../Core/Types.h"
#include "../FiniteVolume/Mesh.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../FiniteVolume/BoundaryConditionData.h"
#include "../Geometry/Geometry.h"

#include <memory>

namespace CAMIRA
{


struct TurbulenceModelInterface
{
    virtual void SetTurbulenceModelData(const InputData &, const Mesh &, const IBData &, const BoundaryConditionData &) = 0; 
    virtual void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const IBData &, const Mesh &) = 0;
};


template< TurbulenceModels T >
struct TurbulenceModel;

// Laminar
template<>
struct TurbulenceModel< TurbulenceModels::Laminar > : public TurbulenceModelInterface
{
    void SetTurbulenceModelData(const InputData &, const Mesh &, const IBData &, const BoundaryConditionData &);
    void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const IBData &, const Mesh &);
};


// Prandtl's mixing length model
template<>
struct TurbulenceModel< TurbulenceModels::PrandtlZeroEquation > : public TurbulenceModelInterface
{
    void SetTurbulenceModelData(const InputData &, const Mesh &, const IBData &, const BoundaryConditionData &);
    void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const IBData &, const Mesh &);

    private:
        floatType m_vonKarmanConstant;
        EnumVector<Axis, Tensor3D> m_wallDistance; // Stored at cell faces
        Tensor3D m_velocityDeformationRate; // Stored at cell centers, to be interpolated to faces
};


// ZEQ0 (Chen and Xu)
template<>
struct TurbulenceModel< TurbulenceModels::ZEQ0 > : public TurbulenceModelInterface
{
    void SetTurbulenceModelData(const InputData &, const Mesh &, const IBData &, const BoundaryConditionData &);
    void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const IBData &, const Mesh &);

    private:
        floatType m_proportionalityConstant;
        EnumVector<Axis, Tensor3D> m_wallDistance; // Stored at cell faces
};


// ZEQ1
template<>
struct TurbulenceModel< TurbulenceModels::ZEQ1 > : public TurbulenceModelInterface
{
    void SetTurbulenceModelData(const InputData &, const Mesh &, const IBData &, const BoundaryConditionData &);
    void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const IBData &, const Mesh &);

    private:
        floatType m_reynoldsNumberBuildingHeight,                   // Inflow reynolds number at average building height
                  m_inflowTurbulenceIntensityBuildingHeight;        // Inflow turbulence intensity at average building height
        EnumVector<Axis, Tensor3D> m_wallDistance; // Stored at cell faces
};


// ZEQ2
template<>
struct TurbulenceModel< TurbulenceModels::ZEQ2 > : public TurbulenceModelInterface
{
    void SetTurbulenceModelData(const InputData &, const Mesh &, const IBData &, const BoundaryConditionData &);
    void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const IBData &, const Mesh &);

    private:
        floatType m_averageBuildingHeight,                            // Average building height
                  m_inflowVelocityBuildingHeight,                     // Inflow velocity magnitude at average building height
                  m_inflowTKEBuildingHeight,                          // Inflow turbulence kinetic energy at average building height
                  m_inflowIntergralTimeScaleBuildingHeight,           // Inflow integral timescale at average building height
                  m_wallDistanceLengthScale,                          // Reference length for wall distance, usually just 1.
                  m_nu;                                               // Kinematic viscosity
        EnumVector<Axis, Tensor3D> m_wallDistance; // Stored at cell faces
};


// ZEQ3
template<>
struct TurbulenceModel< TurbulenceModels::ZEQ3 > : public TurbulenceModelInterface
{
    void SetTurbulenceModelData(const InputData &, const Mesh &, const IBData &, const BoundaryConditionData &);
    void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const IBData &, const Mesh &);

    private:
        Axis::ENUMDATA m_heightAxis;                                  // The axis which corresponds to height from ground level
        floatType m_averageBuildingHeight,                            // Average building height
                  m_inflowVelocityBuildingHeight,                     // Inflow velocity magnitude at average building height
                  m_roughnessLength,                                  // Roughness length
                  m_alpha,                                            // Part in roughness height parameter
                  m_zh,                                               // Part of roughness height parameter
                  m_b;                                                // Urban morphological parameter
        EnumVector<Axis, Tensor3D> m_wallDistance; // Stored at cell faces
};


// ZEQ4
template<>
struct TurbulenceModel< TurbulenceModels::ZEQ4 > : public TurbulenceModelInterface
{
    void SetTurbulenceModelData(const InputData &, const Mesh &, const IBData &, const BoundaryConditionData &);
    void SetTurbulenceViscosityField(EnumVector<Axis, Tensor3D> &, const FieldData<Tensor3D> &, const IBData &, const Mesh &);

    private:
        Axis::ENUMDATA m_heightAxis;                                  // The axis which corresponds to height from ground level
        floatType m_averageBuildingHeight,                            // Average building height
                  m_averageBuildingWidth,                             // Average building width
                  m_referenceHeight,                                  // Reference height for city
                  m_Cmu,                                              // Model constant
                  m_Ig,                                               // Model constant given by AIJ
                  m_alpha;                                            // Ground roughness parameter for city    
        EnumVector<Axis, Tensor3D> m_wallDistance; // Stored at cell faces
        Tensor3D m_velocityDeformationRate; // Stored at cell centers, to be interpolated to faces
};


}   // end namespace CAMIRA


#endif // CAMIRA_TURBULENCE_MODELS