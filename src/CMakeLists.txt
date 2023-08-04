# Add libraries
add_library(libs "")

# Source subdirectories
add_subdirectory(FiniteVolume)
add_subdirectory(IO)
add_subdirectory(Solver)
add_subdirectory(Tools)

# Link libraries
target_link_libraries(libs PUBLIC fmt::fmt)
target_link_libraries(libs PUBLIC eigen)
target_link_libraries(libs PRIVATE warning_flags)
target_link_libraries(libs PUBLIC compiler_flags)

# Main solver executable
add_executable(main main.cpp)
target_link_libraries(main PUBLIC libs)

# Field transformation executable
add_executable(soltransform soltransform.cpp)
target_link_libraries(soltransform PUBLIC libs)