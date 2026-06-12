cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(CAMIRA VERSION 1.0)

include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)

# Default executable location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# specify the C++ standard
set (CMAKE_CXX_STANDARD 20)

# For clangd
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Add compiler flags
add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_20)
target_compile_options(compiler_flags INTERFACE "-std=c++2a;-march=native;-mno-avx512f")

add_library(warning_flags INTERFACE)
target_compile_features(warning_flags INTERFACE cxx_std_20)
target_compile_options(warning_flags INTERFACE "-Wall;-Wextra;-Wshadow;-Wunused;-Wconversion;-Wpedantic")


# Add address sanitizer to the debug build
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")

if ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
	message ( WARNING "THE BUILD TYPE IS DEBUG. ADDRESS SANITIZER WILL BE USED" )
endif()

# Make release the default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" )
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build (Debug, Release etc." FORCE)
endif()

# Profiler
option ( CAMIRA_PROFILING "Enable CAMIRA Profiling Tickers" on )
if ( CAMIRA_PROFILING )
	message( STATUS "CAMIRA profiling tickers are on!" )
	add_compile_definitions( CAMIRA_PROFILING )
endif()


# Selection of solver precision
set(CAMIRA_PRECISION "double" CACHE STRING "Precision of floating point numbers to compile in.")
set_property(CACHE CAMIRA_PRECISION PROPERTY STRINGS single double)

if ( CAMIRA_PRECISION STREQUAL "double" )
	message( STATUS "Compiling in double precision" )
	add_compile_definitions( CAMIRA_DOUBLE_PRECISION )
elseif ( CAMIRA_PRECISION STREQUAL "single" )
	message( STATUS "Compiling in single precision" )
endif()



# Parallel smoother strategy
set(CAMIRA_PARALLEL_COLOURING_STRATEGY "planes" CACHE STRING "Method of colour ordering for parallel smoother")
set_property(CACHE CAMIRA_PARALLEL_COLOURING_STRATEGY PROPERTY STRINGS points lines planes)

if ( CAMIRA_PARALLEL_COLOURING_STRATEGY STREQUAL "points" )

	message( STATUS "Parallel colour strategy: points" )
	add_compile_definitions( CAMIRA_PARALLEL_COLOURING_POINTS )

elseif ( CAMIRA_PARALLEL_COLOURING_STRATEGY STREQUAL "lines" )

	message( STATUS "Parallel colour strategy: lines" )
	add_compile_definitions( CAMIRA_PARALLEL_COLOURING_LINES )

elseif ( CAMIRA_PARALLEL_COLOURING_STRATEGY STREQUAL "planes" )

	message( STATUS "Parallel colour strategy: planes" )
	add_compile_definitions( CAMIRA_PARALLEL_COLOURING_PLANES )

endif()

# Compile the Plume dispersion application
option ( CAMIRA_COMPILE_PLUME "Compile the plume dispersion application. Requires VTK library." on )

# Compile option for using VTK library, which is optional if only flow application is used
option ( CAMIRA_HAS_VTK_LIB "Some features in camira-flow require VTK library to be installed." on )

# Enforce that plume requires vtk to be used
if ( CAMIRA_COMPILE_PLUME AND NOT CAMIRA_HAS_VTK_LIB )
	message( FATAL_ERROR "Compiling the dispersion application (CAMIRA_COMPILE_PLUME) requires VTK library (CAMIRA_HAS_VTK_LIB)." )
endif()


# Allow the ability to read fields in (for example to read in an initial condition), this requires VTK library
if ( CAMIRA_HAS_VTK_LIB )

	find_package( VTK COMPONENTS 
				CommonCore
				CommonExecutionModel
				IOXMLParser
				CommonDataModel
				CommonMisc
				CommonSystem
				IOCore
				vtksys
				IOLegacy
				)

	if ( ${VTK_FOUND} )
		message( STATUS "VTK library found successfully." )
		add_compile_definitions( CAMIRA_HAS_VTK_LIB )
	else()
		message( FATAL_ERROR "Could not locate VTK library installation! This is required for some features, but not essential!" )
	endif()

else()
	message( STATUS "VTK library will not be used. Some features will be limited." )
endif()



# Add libraries
include_directories(SYSTEM libraries)	#  SYSTEM to suppress warnings

# Supress eigen warnings
add_library(eigen_internal INTERFACE)
target_include_directories(eigen_internal SYSTEM INTERFACE ${CMAKE_SOURCE_DIR}/libraries/eigen-3.4.0)
target_compile_options(eigen_internal INTERFACE -Wno-array-bounds)

# Public Eigen target without warnings suppressed
add_library(eigen INTERFACE)
target_link_libraries(eigen INTERFACE eigen_internal)

# Openmp, so that it works with clang
find_package(OpenMP REQUIRED)

# CGAL library
find_package(CGAL REQUIRED 
			 COMPONENTS Qt5
			 PATHS libraries/CGAL-5.6 NO_DEFAULT_PATH )
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

add_subdirectory(src)

# Link time optimisation
if ( supported )
    message(STATUS "IPO / LTO enabled")
    set_property(TARGET camira-flow PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
	
	if ( CAMIRA_COMPILE_PLUME )
		set_property(TARGET camira-plume PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
	endif()
	
else()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif()


