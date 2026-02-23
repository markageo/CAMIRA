# Source subdirectories
add_subdirectory(InputProcessing)

target_sources( libs 
    PRIVATE 
        ConfigEnums.h )

# Main solver executable
add_executable(camira-plume main.cpp)
target_link_libraries(camira-plume PUBLIC libs)