#ifndef CFD_MULTIGRID
#define CFD_MULTIGRID

#include "../Types.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../IO/InputProcessing.h"
#include "../Solver/LinearSolver.h"

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
    
    intType level;
    bool isCoarsestLevel, isFinestLevel;
    Mesh mesh;
    FieldData<Tensor3D> fields,
                        fieldsOld,
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
                  const InputData & );


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