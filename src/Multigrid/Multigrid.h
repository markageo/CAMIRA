#ifndef CAMIRA_MULTIGRID
#define CAMIRA_MULTIGRID

#include "../Core/Types.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../TurbulenceModels/TurbulenceModelData.h"
#include "../IO/InputProcessing.h"
#include "../Solver/LinearSolver.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"

#include<memory>


namespace CAMIRA
{

enum class MultigridEquation {
    NoTauCorrection,
    TauCorrection
};

template< MomentumInterpolation MI >
struct GridLevelData
{
    
    size_t level;
    bool isCoarsestLevel, isFinestLevel;
    Mesh mesh;
    FieldData<Tensor3D> fields,                     // Solution at latest iteration (steady) or timestep t (transient)
                        fieldsOld,                  // Old estimate of fields, used for under-relaxation
                        fieldsPrevTime,             // Solution at previous timestep t-1 (only used in transient)
                        fieldsPrevPrevTime,         // Solution at timestep before previous timestep t-2 (only used in transient)
                        fieldsRestricted,
                        residuals,
                        residualsRestricted,
                        fineGridCorrection,
                        coarseGridRightHandSide;    // Only used on coarse grids
    EnumVector<Axis, Tensor3D> faceFluxes;
    BoundaryConditionData bcData;
    IBData ibData;
    TurbulenceModelData turbModelData;
    std::unique_ptr< FVCoefficients > fvCoeffsPtr;  
    std::unique_ptr< LinearSolverInterface<MI> > linearSolverPtr;
};



// Create initial heirachy of grids
template< MomentumInterpolation MI >
void SetMGLevels( std::vector< GridLevelData<MI > > &, 
                  const InputData &,
                  const AxisTransformationMap & );


void RestrictField( Tensor3D &,
                    const Tensor3D &, 
                    const Mesh &,
                    const Mesh & );


void ProlongateField( Tensor3D &,
                      const Tensor3D &,
                      const Mesh &,
                      const Mesh & );

void RestrictFields( FieldData<Tensor3D> &,
                     const FieldData<Tensor3D> &,
                     const Mesh &, 
                     const Mesh &, 
                     const Tensor3D &mask );

void ComputeFineGridCorrection( FieldData<Tensor3D> &,
                                const FieldData<Tensor3D> &,
                                const FieldData<Tensor3D> &,
                                const Mesh &,
                                const Mesh &,
                                const Tensor3D & );


template< MomentumInterpolation MI >
void TransformToCoarseGridEquations( FVCoefficients &,
                                     const FieldData<Tensor3D> &,
                                     const FieldData<Tensor3D> &,
                                     const Tensor3D & );


template< MomentumInterpolation MI >
void CalculateCoarseGridRightHandSide( FieldData<Tensor3D> &,
                                       const FVCoefficients &,
                                       const FieldData<Tensor3D> &,
                                       const FieldData<Tensor3D> &,
                                       const Tensor3D & );

}

#endif  // CAMIRA_MULTIGRID