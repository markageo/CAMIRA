target_sources( libs 
    PRIVATE 
        FaceVelocities.cpp 
        FiniteVolumeCoefficients.cpp 
        InitialConditions.cpp
        BoundaryConditions.cpp 
        FVCoefficients.cpp
        VertexValues.cpp
    PRIVATE 
        FVCoefficients.h
        BoundaryConditionData.h
        FiniteVolumeFunctions.h
        FiniteVolume.h
        GhostCells.h
        FaceInterpolatedVelocity.h )