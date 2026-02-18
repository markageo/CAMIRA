add_library(libs "")

# Link libraries
target_link_libraries(libs PUBLIC eigen)
target_link_libraries(libs PUBLIC CGAL::CGAL)
target_link_libraries(libs PUBLIC OpenMP::OpenMP_CXX)

if ( CAMIRA_HAS_VTK_LIB )
    target_link_libraries(libs PUBLIC ${VTK_LIBRARIES})
    vtk_module_autoinit( TARGETS libs
                         MODULES ${VTK_LIBRARIES} )
endif()

target_link_libraries(libs PRIVATE warning_flags)
target_link_libraries(libs PUBLIC compiler_flags)

add_subdirectory(Core)

add_subdirectory(Flow)

if ( CAMIRA_COMPILE_PLUME )
	add_subdirectory(Plume)
endif()



