
# Source subdirectories
add_subdirectory(ImmersedBoundary)
add_subdirectory(DerivedQuantities)
add_subdirectory(FiniteVolume)
add_subdirectory(Multigrid)
add_subdirectory(InputProcessing)
add_subdirectory(Solver)
add_subdirectory(CoordinateTransformations)
add_subdirectory(Parallel)
add_subdirectory(TurbulenceModels)

# Main solver executable
add_executable(camira-flow main.cpp)
target_link_libraries(camira-flow PUBLIC libs)