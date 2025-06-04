# Add libraries
add_library(implementation_libs "")

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


# External dependencies
set(application_depends eigen
                        CGAL::CGAL CGAL::CGAL_Basic_viewer
                        camp
                        RAJA
                        umpire )

if ( CAMIRA_HAS_VTK_LIB )
    set(application_depends ${application_depends} ${VTK_LIBRARIES})
endif ()

if ( ENABLE_OPENMP )
    set(application_depends ${application_depends} openmp)
endif ()

if ( ENABLE_CUDA )
    set(application_depends ${application_depends} cuda)
    set(CUDA_PROPAGATE_HOST_FLAGS OFF)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --disable-warnings")
endif ()

target_link_libraries(implementation_libs PUBLIC ${application_depends})

blt_add_executable( NAME camira
                    SOURCES main.cpp
                    DEPENDS_ON ${application_depends})


target_link_libraries(camira PUBLIC implementation_libs)

blt_add_target_compile_flags( TO camira
                              SCOPE PUBLIC 
                              FLAGS "-Wall;-Wextra;-Wshadow;-Wunused;-Wconversion;-Wpedantic")

# TODO, figure out how to add -march=native and -mno-avc512f when using CUDA. Second one is needed for Eigen to avoid issues.
# blt_add_target_compile_flags( TO camira
#                               SCOPE PRIVATE 
#                               FLAGS "-std=c++20;-march=native;-mno-avx512f")

