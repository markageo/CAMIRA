# Add libraries
add_library(libs InputParser.cpp)

# Add the executable
add_executable(main main.cpp)

# Add include directories
target_include_directories(libs PUBLIC ../include)

target_include_directories(libs PUBLIC ../Eigen)
target_include_directories(libs PUBLIC ../Eigen/unsupported)
target_include_directories(libs PUBLIC ../Eigen/unsupported/CXX11)

target_include_directories(libs PUBLIC ../Boost/boost/property_tree)

target_include_directories(libs PUBLIC ../include)

# Link libraries
target_link_libraries(main PUBLIC compiler_flags)
target_link_libraries(libs PUBLIC compiler_flags)
target_link_libraries(main PUBLIC libs)

