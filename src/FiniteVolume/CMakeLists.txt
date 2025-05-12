target_sources( libs 
    PRIVATE 
        FaceVelocities.cpp 
        FiniteVolumeCoefficients.cpp 
        FiniteVolumeStructures.cpp 
        VertexValues.cpp
        Mesh.cpp 
    PRIVATE 
        FiniteVolumeStructures.h
        FiniteVolumeFunctions.h
        FiniteVolume.h
        Mesh.h
        GhostCells.h
        FaceInterpolatedVelocity.h )