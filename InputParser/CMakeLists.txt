# Add libraries
add_library(inputparser InputParser.cpp)

# Add include directories
target_include_directories(inputparser PUBLIC ../Boost/boost/property_tree)

# Link libraries
target_link_libraries(inputparser PUBLIC compiler_flags)

