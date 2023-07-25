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
target_link_libraries(libs PRIVATE compiler_flags)

# Main executable
add_executable(main main.cpp)
target_link_libraries(main PUBLIC libs)