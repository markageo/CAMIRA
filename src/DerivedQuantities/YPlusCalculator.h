#ifndef CAMIRA_YPLUS_CALCULATOR
#define CAMIRA_YPLUS_CALCULATOR

#include "../Core/Types.h"
#include "../ImmersedBoundary/ImmersedBoundary.h"
#include "../FiniteVolume/Mesh.h"
#include "../FiniteVolume/BoundaryConditionData.h"

namespace CAMIRA
{

// Calculates some y+ statistics
class YPlusCalculator
{
    public: 
        YPlusCalculator( const InputData &,
                         const AxisTransformationMap &,
                         const BoundaryConditionData &, 
                         const IBData &, 
                         const Mesh &, 
                         const FieldData<Tensor3D> & );

        void Update();
        floatType minYPlus, maxYPlus, averageYPlus;

    private:
        const FieldData<Tensor3D> &m_fields;
        const floatType m_rho, m_nu;

        struct WallCellData {
            floatType yPlus;
            TensorIndex3D cellIndex;
            floatType wallDistance;
            fVector3 normalVector;
            fVector3 wallTangentialVelocity;
        };

        std::vector< WallCellData > m_wallCells;
};



}   // end namespace CAMIRA

#endif // CAMIRA_YPLUS_CALCULATOR