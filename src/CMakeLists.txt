# Add libraries
add_library(libs InputProcessing.cpp
                 VTKWriter.cpp
                 FiniteVolumeStructures.cpp
                 FaceVelocities.cpp
                 FVCoefficients.cpp
                 Solver.cpp)

# Add the executable
add_executable(main main.cpp)

# Link libraries
target_link_libraries(main PUBLIC compiler_flags)
target_link_libraries(libs PUBLIC compiler_flags)

target_link_libraries(libs PUBLIC inputparser)
target_link_libraries(main PUBLIC libs)
