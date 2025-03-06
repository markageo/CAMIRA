# Add libraries
add_library(libs "")

# Source subdirectories
add_subdirectory(ImmersedBoundary)
add_subdirectory(DerivedQuantities)
add_subdirectory(FiniteVolume)
add_subdirectory(Multigrid)
add_subdirectory(Geometry)
add_subdirectory(IO)
add_subdirectory(Solver)
add_subdirectory(Core)
add_subdirectory(CoordinateTransformations)
add_subdirectory(Parallel)

# Link libraries
target_link_libraries(libs PUBLIC eigen)
target_link_libraries(libs PUBLIC RAJA)
target_link_libraries(libs PUBLIC CGAL::CGAL CGAL::CGAL_Basic_viewer)

if ( CFD_HAS_VTK_LIB )
    target_link_libraries(libs PUBLIC ${VTK_LIBRARIES})
    vtk_module_autoinit( TARGETS libs
                         MODULES ${VTK_LIBRARIES} )
endif()

target_link_libraries(libs PRIVATE warning_flags)
target_link_libraries(libs PUBLIC compiler_flags)

# Main solver executable
add_executable(camira main.cpp)
target_link_libraries(camira PUBLIC libs)


