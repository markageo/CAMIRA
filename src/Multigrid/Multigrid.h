#ifndef CFD_MULTIGRID
#define CFD_MULTIGRID

#include "../Core/Types.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolumeStructures.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../IO/InputProcessing.h"
#include "../Solver/LinearSolver.h"
#include "../CoordinateTransformations/AxisTransformationFunctions.h"

#include<memory>


namespace CFD
{

enum class MultigridEquation {
    NoTauCorrection,
    TauCorrection
};

template< MomentumInterpolation MI, Linearisation LI >
struct GridLevelData
{
    
    size_t level;
    bool isCoarsestLevel, isFinestLevel;
    Mesh mesh;
    FieldData<Tensor3D> fields,                 // Solution at latest iteration (steady) or timestep t (transient)
                        fieldsOld,              // Old estimate of fields, used for under-relaxation
                        fieldsPrevTime,         // Solution at previous timestep t-1 (only used in transient)
                        fieldsPrevPrevTime,     // Solution at timestep before previous timestep t-2 (only used in transient)
                        fieldsRestricted,
                        residualsRestricted;
    BoundaryConditionData bcData;
    EnumVector<Axis, Tensor3D> faceFluxes;
    EnumVector< Axis, EnumVector< Axis, Tensor3D> > faceAdvectedVelocities;
    IBData ibData;
    FVCoefficients fvCoeffs;
    std::unique_ptr< LinearSolverInterface<MI, LI> > linearSolver;
};



// Create initial heirachy of grids
template< MomentumInterpolation MI, Linearisation LI >
void SetMGLevels( std::vector< GridLevelData<MI, LI> > &, 
                  const InputData &,
                  const AxisTransformationMap & );


Tensor3D RestrictField( const Tensor3D &, 
                        const Mesh &,
                        const Mesh & );


Tensor3D ProlongateField( const Tensor3D &,
                          const Mesh &,
                          const Mesh & );

FieldData<Tensor3D> RestrictFields( const FieldData<Tensor3D>,
                                   const Mesh &, 
                                   const Mesh &, 
                                   const Tensor3D &mask );

FieldData<Tensor3D> ComputeFineGridCorrection( const FieldData<Tensor3D> &,
                                               const FieldData<Tensor3D> &,
                                               const Mesh &,
                                               const Mesh &,
                                               const Tensor3D & );


template< MomentumInterpolation MI, Linearisation LI >
void TransformToCoarseGridEquations( FVCoefficients &,
                                     const FieldData<Tensor3D> &,
                                     const FieldData<Tensor3D> &,
                                     const Tensor3D & );


template< MomentumInterpolation MI, Linearisation LI >
FieldData<Tensor3D> CalculateCoarseGridRightHandSide( FVCoefficients &,
                                                      const FieldData<Tensor3D> &,
                                                      const FieldData<Tensor3D> &,
                                                      const Tensor3D & );

}

#endif  // CFD_MULTIGRID