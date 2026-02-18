#ifndef CAMIRA_FVCOEFFICIENTS  
#define CAMIRA_FVCOEFFICIENTS

#include "../../Core/Types.h"
#include "../TurbulenceModels/TurbulenceModels.h"
#include <memory>

namespace CAMIRA
{

// -------------------------------------- Definition in FiniteVolumeStructures.cpp -------------------------------------- //


// Momentum and continuity equations may share some memory. This shared memory is kept here.
class FVCoefficients
{

private:
    // Velocity coefficients are shared for each momentum equation since they are the same, and this has significant memory savings.
    // To allow a seperated interface between momentum equation, they all hold references to this variable.
    EnumVector< TransportCoefficients, Tensor3D > m_momentumVelocityCoeffs;

    // Pressure gradient coefficients in momentum equations are same as divergence coefficients in continuity equation (apart from 
    // density, which is absorbed into the pressure field). They will hold references to this variable. 
    EnumVector< Axis, EnumVector< TransportCoefficients, Tensor1D > > m_gradientCoeffs;


public:

    FVCoefficients();
    FVCoefficients(const iArray3 &, MomentumInterpolation);

    // Delete copy and move constructors
    // Copy and move for this class are not neccessary
    // In the case of NUMA aware allocations, copy and move can break the NUMA layout
    FVCoefficients(const FVCoefficients &) = delete;
    FVCoefficients &operator=(const FVCoefficients &) = delete;
    FVCoefficients(FVCoefficients &&) = delete;
    FVCoefficients &operator=(FVCoefficients &&) = delete;



    // Continuity equation data
    struct ContinuityEquation
    {
        EnumVector< Axis, EnumVector< TransportCoefficients, Tensor1D > > &AU;   // Velocity coefficients (LHS), has dummy cells
        EnumVector<TransportCoefficients, Tensor3D> AP;                          // Pressure coefficients (LHS), has dummy cells
        Tensor3D B;                                                              // Constants that come from boundary conditions and linearisation (LHS), has dummy cells
        Tensor3D F;                                                              // Source terms (RHS), has dummy cells
    };
    ContinuityEquation Cont;

    // Continuity equation auxiliary data
    EnumVector< Axis, std::array< Tensor1D, 4 > > mwiSparseCoeffs;               // Unweighted MWI coefficients from the sparse pressure gradient (LHS)
    EnumVector< Axis, std::array< Tensor1D, 2 > > mwiCompactCoeffs;              // Unweighted MWI coefficients from the compact pressure gradient (LHS)
    MomentumInterpolation momentumInterpolation;


    // Momentum equation data
    struct MomentumEquation
    {
        Axis::ENUMDATA component;
        EnumVector<TransportCoefficients, Tensor3D > &AU;                        // Velocity coefficients (LHS), has dummy cells
        EnumVector<TransportCoefficients, Tensor1D> &AP;                         // Pressure coefficients (LHS), has dummy cells
        Tensor3D B;                                                              // Constants that come from boundary conditions, immersed boundary, etc., has dummy cells
        Tensor3D F;                                                              // Source terms (RHS), has dummy cells
    };
    EnumVector<Axis, MomentumEquation> Mom;

    // Momentum equation auxiliary data
    struct HiOrderAdvectionCoeffs {                                              // Precomputed high
        // SOU     : phi_f = g1 * phi_U  +  g2 * phi_UU
        // QUICK   : phi_f = phi_U  +  g1 * ( phi_D - phi_U )  +  g2 * ( phi_U - phi_UU )
        EnumVector< Axis, Tensor1D > g1, g2;                        
    };
    HiOrderAdvectionCoeffs positiveFluxHiOrderAdvectionCoeffs, negativeFluxHiOrderAdvectionCoeffs;
    AdvectionSchemes advectionScheme;
    TimeSchemes timeScheme;
    floatType timeStep;
    floatType advectionBlendingFactor;


    // Turbulence Viscosity
    Tensor3D nuTurb;                                                            // Stored at cell centers, has ghost cells for BCs
    

    // General data
    floatType nu, rho;
    iArray3 nCells;
};



} // end namespace CAMIRA

#endif // CAMIRA_FVCOEFFICIENTS