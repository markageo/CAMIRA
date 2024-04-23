#ifndef CFD_MULTIGRID
#define CFD_MULTIGRID

#include "../Types.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/FiniteVolume.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../IO/InputProcessing.h"
#include "../Solver/LinearSolver.h"


namespace CFD
{

template< MomentumInterpolation MI, Linearisation LI >
struct GridLevelData
{
    intType level;
    bool isCoarsestLevel, isFinestLevel;
    Mesh mesh;
    FieldData<Tensor3D> fields,
                        fieldsOld,
                        fieldsRestricted,
                        residuals;
    BoundaryConditionData bcData;
    EnumVector<Axis, Tensor3D> faceFluxes;
    EnumVector< Axis, EnumVector< Axis, Tensor3D> > faceAdvectedVelocities;
    IBData ibData;
    FVCoefficients fvCoeffs;
    LinearSolver<MI, LI> *linearSolver;
};



// Create initial heirachy of grids
template< MomentumInterpolation MI, Linearisation LI >
std::vector< GridLevelData<MI, LI> > CreateMGLevels( const InputData & );


Tensor3D RestrictField( const Tensor3D &, 
                        const Mesh &,
                        const Mesh & );

Tensor3D ProlongateField( const Tensor3D &,
                          const Mesh &,
                          const Mesh & );

}

#endif  // CFD_MULTIGRID