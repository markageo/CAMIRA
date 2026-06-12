#ifndef CAMIRA_FLOW_CONFIG_ENUMS
#define CAMIRA_FLOW_CONFIG_ENUMS

namespace CAMIRA 
{

namespace FLOW
{

enum class BoundaryConditions {
    zeroGradient, fixed, extrapolated, periodic
};

enum class Smoothers {
    nestedLineSymmetricSerial, domainSymmetricSerial, domainSymmetricParallel
};

enum class GeometryBoundaryTreatment {
    Staircase, DirectionalImmersedBoundary
};

enum class MomentumInterpolation {
    Implicit, SemiExplicit
};

enum class AdvectionSchemes {
    Upwind, Central, SOU, QUICK
};

enum class TimeSchemes {
    Steady, BackwardsEuler, BackwardsThreeLevel
};

enum class TurbulenceModels {
    Null, Laminar, PrandtlZeroEquation, ZEQ0, ZEQ1, ZEQ2, ZEQ3, ZEQ4
};

enum class MultigridCycleType {
    V, F, W
};


}   // end namespace FLOW

}   // end namespace CAMIRA

#endif // CAMIRA_FLOW_CONFIG_ENUMS