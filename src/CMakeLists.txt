# Add libraries
add_library(libs InputParser.cpp
                 InputProcessing.cpp
                 SweepTransformations.cpp
                 FiniteVolumeStructures.cpp
                #  FaceVelocities.cpp
                #  FiniteVolumeCoefficients.cpp
                #  Solver.cpp)
                )

# Model library
add_library(model INTERFACE)

target_link_libraries(libs PUBLIC fmt::fmt)
target_link_libraries(libs PUBLIC eigen)

# Model library
target_link_libraries(model INTERFACE compiler_flags)
target_link_libraries(libs PRIVATE compiler_flags)


target_link_libraries(model INTERFACE libs)
target_link_libraries(model INTERFACE fmt::fmt)
target_link_libraries(model INTERFACE eigen)

# Main executable
add_executable(main main.cpp)
target_link_libraries(main PUBLIC model)
