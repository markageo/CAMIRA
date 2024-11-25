#ifndef CFD_FINITE_VOLUME_STRUCTURES   
#define CFD_FINITE_VOLUME_STRUCTURES

#include "../Types.h"

namespace CFD
{

// -------------------------------------- Definition in FiniteVolumeStructures.cpp -------------------------------------- //


struct MomentumEquation {
    MomentumEquation() {};
    MomentumEquation(const Axis::ENUMDATA, const iArray3 &, Linearisation);
    EnumVector< Axis, EnumVector< TransportCoefficients, Tensor3D > > AU;     // Velocity coefficients (LHS), has dummy cells
    EnumVector<TransportCoefficients, Tensor1D> AP;                           // Pressure coefficients (LHS), has dummy cells
    Tensor3D B;                                                               // Constants that come from boundary conditions and linearisation (LHS), has dummy cells
    Tensor3D F;                                                               // Source terms (RHS), has dummy cells
    Tensor3D diagCoeffInv;                                                    // Inverse of diagonal coefficient, has dummy cells
    EnumVector< Axis, EnumVector<TransportCoefficients, Tensor1D> > diff;     // Diffusion coefficients (LHS)
    EnumVector< BoundaryPatches, floatType > diffBoundary;                    // Diffusion coefficients for constant boundary conditions (LHS)
    EnumVector< BoundaryPatches, Tensor2D   > BUBoundary, BPBoundary;         // Constant terms that come from fixed BC (LHS)
    struct HiOrderAdvectionCoeffs {                                           // Precomputed high
        // SOU     : phi_f = g1 * phi_U  +  g2 * phi_UU
        // QUICK   : phi_f = phi_U  +  g1 * ( phi_D - phi_U )  +  g2 * ( phi_U - phi_UU )
        EnumVector< Axis, Tensor1D > g1, g2;                        
    };
    HiOrderAdvectionCoeffs positiveFluxHiOrderAdvectionCoeffs, negativeFluxHiOrderAdvectionCoeffs;
    floatType relaxation;
    Axis::ENUMDATA component;                                                 // The momentum component
    Linearisation linearisation;
    AdvectionSchemes advectionScheme;
    TimeSchemes timeScheme;
    floatType timeStep;
    floatType advectionBlendingFactor;
};


struct ContinuityEquation {
    ContinuityEquation() {};
    ContinuityEquation(const iArray3 &, MomentumInterpolation);
    EnumVector< Axis, EnumVector< TransportCoefficients, Tensor1D > > AU;    // Velocity coefficients (LHS), has dummy cells
    EnumVector<TransportCoefficients, Tensor3D> AP;                          // Pressure coefficients (LHS), has dummy cells
    Tensor3D B;                                                              // Constants that come from boundary conditions and linearisation (LHS), has dummy cells
    Tensor3D F;                                                              // Source terms (RHS), has dummy cells
    EnumVector< Axis, std::array< Tensor1D, 4 > > mwiSparseCoeffs;           // Unweighted MWI coefficients from the sparse pressure gradient (LHS)
    EnumVector< Axis, std::array< Tensor1D, 2 > > mwiCompactCoeffs;          // Unweighted MWI coefficients from the compact pressure gradient (LHS)
    EnumVector< BoundaryPatches, Tensor2D   > BUBoundary, BPBoundary;        // Constant terms that come from fixed BC (LHS)
    floatType relaxation; 
    MomentumInterpolation momentumInterpolation;
};


// Structure to store all domain boundary condition information
struct BoundaryConditionData {
    struct Patch {
        BoundaryConditions::ENUMDATA type;
        Tensor2D value;
    };
    using Patches = EnumVector< BoundaryPatches, Patch >;
    FieldData< Patches > fields; 
    bool pressureFieldIsFloating;
};


// Structure to store finite volume discrete equation coefficients
struct FVCoefficients
{
    FVCoefficients() {};
    FVCoefficients(const iArray3 &, Linearisation, MomentumInterpolation);
    
    EnumVector<Axis, MomentumEquation > Mom;
    ContinuityEquation Cont;
    iArray3 nCells;
};


} // end namespace CFD

#endif // CFD_FINITE_VOLUME_STRUCTURES