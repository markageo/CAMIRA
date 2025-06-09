#ifndef CAMIRA_FINITE_VOLUME_STRUCTURES   
#define CAMIRA_FINITE_VOLUME_STRUCTURES

#include "../Core/Types.h"

namespace CAMIRA
{

// -------------------------------------- Definition in FiniteVolumeStructures.cpp -------------------------------------- //


// Momentum and continuity equations may share some memory. This shared memory is kept here.
// Uses Views as the underlying array data structure.
class FVCoefficients
{

private:
    // Velocity coefficients are shared for each momentum equation since they are the same, and this has significant memory savings.
    // To allow a seperated interface between momentum equation, they all hold references to this variable.
    EnumVector< TransportCoefficients, View3D > m_momentumVelocityCoeffs;

    // Pressure gradient coefficients in momentum equations are same as divergence coefficients in continuity equation (apart from 
    // density, which is absorbed into the pressure field). They will hold references to this variable. 
    EnumVector< Axis, EnumVector< TransportCoefficients, View1D > > m_gradientCoeffs;


public:

    // Deleted since RAJA::Views are not copy assignable, they contant a const member.
    FVCoefficients(const std::string &, const iArray3 &, MomentumInterpolation);
    FVCoefficients(const FVCoefficients &);
    FVCoefficients(FVCoefficients&&) noexcept;
    ~FVCoefficients();


    // Continuity equation data
    struct ContinuityEquation
    {
        EnumVector< Axis, EnumVector< TransportCoefficients, View1D > > &AU;   // Velocity coefficients (LHS), has dummy cells
        EnumVector<TransportCoefficients, View3D> AP;                          // Pressure coefficients (LHS), has dummy cells
        View3D B;                                                              // Constants that come from boundary conditions and linearisation (LHS), has dummy cells
        View3D F;                                                              // Source terms (RHS), has dummy cells
    };
    ContinuityEquation Cont;

    // Continuity equation auxiliary data
    EnumVector< Axis, View1D[4] > mwiSparseCoeffs;                             // Unweighted MWI coefficients from the sparse pressure gradient (LHS)
    EnumVector< Axis, View1D[2] > mwiCompactCoeffs;                            // Unweighted MWI coefficients from the compact pressure gradient (LHS)
    MomentumInterpolation momentumInterpolation;


    // Momentum equation data
    struct MomentumEquation
    {
        Axis::ENUMDATA component;
        EnumVector<TransportCoefficients, View3D > &AU;                        // Velocity coefficients (LHS), has dummy cells
        EnumVector<TransportCoefficients, View1D> &AP;                         // Pressure coefficients (LHS), has dummy cells
        View3D B;                                                              // Constants that come from boundary conditions, immersed boundary, etc., has dummy cells
        View3D F;                                                              // Source terms (RHS), has dummy cells
    };
    EnumVector<Axis, MomentumEquation> Mom;

    // Momentum equation auxiliary data
    struct HiOrderAdvectionCoeffs {                                              // Precomputed high
        // SOU     : phi_f = g1 * phi_U  +  g2 * phi_UU
        // QUICK   : phi_f = phi_U  +  g1 * ( phi_D - phi_U )  +  g2 * ( phi_U - phi_UU )
        EnumVector< Axis, View1D > g1, g2;                        
    };
    HiOrderAdvectionCoeffs positiveFluxHiOrderAdvectionCoeffs, negativeFluxHiOrderAdvectionCoeffs;
    AdvectionSchemes advectionScheme;
    TimeSchemes timeScheme;
    floatType timeStep;
    floatType advectionBlendingFactor;


    // General data
    floatType nu, rho;
    iArray3 nCells;
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



} // end namespace CAMIRA

#endif // CAMIRA_FINITE_VOLUME_STRUCTURES