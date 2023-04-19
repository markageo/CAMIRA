# Add libraries
add_library(libs InputProcessing.cpp
                 VTKWriter.cpp
                 FiniteVolumeStructures.cpp
                 FaceVelocities.cpp
                 FVCoefficients.cpp
                 Solver.cpp)

# Model library
add_library(model INTERFACE)

# Model library
target_link_libraries(model INTERFACE compiler_flags)
target_link_libraries(libs PUBLIC compiler_flags)

target_link_libraries(libs PUBLIC inputparser)
target_link_libraries(model INTERFACE libs)
target_link_libraries(model INTERFACE fmt::fmt)

# Main executable
add_executable(main main.cpp)
target_link_libraries(main PUBLIC model)
