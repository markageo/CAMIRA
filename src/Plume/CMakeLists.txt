# Source subdirectories
add_subdirectory(InputProcessing)
add_subdirectory(Particle)
add_subdirectory(Sources)
add_subdirectory(Solver)
add_subdirectory(Concentration)

target_sources( libs 
    PRIVATE 
        ConfigEnums.h )

# Main solver executable
add_executable(camira-plume main.cpp)
target_link_libraries(camira-plume PUBLIC libs)