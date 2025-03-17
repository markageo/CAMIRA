#ifndef CFD_FINITE_VOLUME_STRUCTURES   
#define CFD_FINITE_VOLUME_STRUCTURES

#include "../Core/Types.h"

namespace CFD
{

// -------------------------------------- Definition in FiniteVolumeStructures.cpp -------------------------------------- //

struct MomentumEquations {

    MomentumEquations();
    MomentumEquations(const iArray3 &);
    MomentumEquations(const MomentumEquations &);
    MomentumEquations &operator=(MomentumEquations);
    MomentumEquations(MomentumEquations&&) noexcept;

private:
    // Velocity coefficients are shared for each momentum equation since they are the same, and this has significant memory savings.
    // To allow a seperated interface between momentum equation, they all hold references to this variable.
    EnumVector< TransportCoefficients, Tensor3D > m_AU;
    Tensor3D m_diagCoeffInv;     

public:

    struct Coeffs {
        Axis::ENUMDATA component;
        EnumVector<TransportCoefficients, Tensor3D > &AU;                    // Velocity coefficients (LHS), has dummy cells
        EnumVector<TransportCoefficients, Tensor1D> AP;                      // Pressure coefficients (LHS), has dummy cells
        Tensor3D B;                                                          // Constants that come from boundary conditions, immersed boundary, etc., has dummy cells
        Tensor3D F;                                                          // Source terms (RHS), has dummy cells
        Tensor3D &diagCoeffInv;                                              // Inverse of diagonal coefficient, has dummy cells
    };
    EnumVector<Axis, Coeffs> coeffs;

    EnumVector< Axis, EnumVector<TransportCoefficients, Tensor1D> > diff;     // Diffusion coefficients (LHS)
    struct HiOrderAdvectionCoeffs {                                           // Precomputed high
        // SOU     : phi_f = g1 * phi_U  +  g2 * phi_UU
        // QUICK   : phi_f = phi_U  +  g1 * ( phi_D - phi_U )  +  g2 * ( phi_U - phi_UU )
        EnumVector< Axis, Tensor1D > g1, g2;                        
    };
    HiOrderAdvectionCoeffs positiveFluxHiOrderAdvectionCoeffs, negativeFluxHiOrderAdvectionCoeffs;
    AdvectionSchemes advectionScheme;
    TimeSchemes timeScheme;
    floatType timeStep;
    floatType advectionBlendingFactor;
};


struct ContinuityEquation {
    ContinuityEquation() {};
    ContinuityEquation(const iArray3 &, MomentumInterpolation);

    struct Coeffs {
        EnumVector< Axis, EnumVector< TransportCoefficients, Tensor1D > > AU;    // Velocity coefficients (LHS), has dummy cells
        EnumVector<TransportCoefficients, Tensor3D> AP;                          // Pressure coefficients (LHS), has dummy cells
        Tensor3D B;                                                              // Constants that come from boundary conditions and linearisation (LHS), has dummy cells
        Tensor3D F;                                                              // Source terms (RHS), has dummy cells
    };
    Coeffs coeffs;

    EnumVector< Axis, std::array< Tensor1D, 4 > > mwiSparseCoeffs;               // Unweighted MWI coefficients from the sparse pressure gradient (LHS)
    EnumVector< Axis, std::array< Tensor1D, 2 > > mwiCompactCoeffs;              // Unweighted MWI coefficients from the compact pressure gradient (LHS)
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
    FVCoefficients(const iArray3 &, MomentumInterpolation);
    
    MomentumEquations Mom;
    ContinuityEquation Cont;
    floatType nu, rho;
    iArray3 nCells;
};


} // end namespace CFD

#endif // CFD_FINITE_VOLUME_STRUCTURES