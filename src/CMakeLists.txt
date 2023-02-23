# Add libraries
# add_library(libs )

# Add the executable
add_executable(main main.cpp)

# Add include directories
# target_include_directories(libs PUBLIC ../include)
# target_include_directories(libs PUBLIC ../include)

# Link libraries
target_link_libraries(main PUBLIC compiler_flags)
# target_link_libraries(libs PUBLIC compiler_flags)
# target_link_libraries(main PUBLIC libs)

