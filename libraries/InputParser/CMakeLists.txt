# Add libraries
add_library(inputparser InputParser.cpp)

# Link libraries
target_link_libraries(inputparser PUBLIC compiler_flags)

