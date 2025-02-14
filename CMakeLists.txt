cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(Coupled-Sweep-NS VERSION 1.0)

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
target_compile_options(compiler_flags INTERFACE "-std=c++2a;-march=native")

add_library(warning_flags INTERFACE)
target_compile_features(warning_flags INTERFACE cxx_std_20)
target_compile_options(warning_flags INTERFACE "-Wall;-Wextra;-Wshadow;-Wunused;-Wconversion;-Wpedantic")


# Add address sanitizer to the debug build
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")

if ( CMAKE_BUILD_TYPE STREQUAL "DEBUG" )
	message ( WARNING "THE BUILD TYPE IS DEBUG. ADDRESS SANITIZER WILL BE USED" )
endif()

# Make release the default build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" )
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build (Debug, Release etc." FORCE)
endif()

# Ankith's profiler
option ( CFD_PROFILING "Enable CFD Profiling Tickers" on )
if ( CFD_PROFILING )
	message( STATUS "CFD profiling tickers are on!" )
	add_compile_definitions( CFD_PROFILING )
endif()


# Selection of solver precision
set(CFD_PRECISION "double" CACHE STRING "Precision of floating point numbers to compile in.")
set_property(CACHE CFD_PRECISION PROPERTY STRINGS single double)

if ( CFD_PRECISION STREQUAL "double" )
	message( STATUS "Compiling in double precision" )
	add_compile_definitions( CFD_DOUBLE_PRECISION )
elseif ( CFD_PRECISION STREQUAL "single" )
	message( STATUS "Compiling in single precision" )
endif()


# Use of functions which were written for the compiler to autovectorise
option ( CFD_USE_AUTOVEC_FUNCTIONS "Use functions written for the compiler to auto vectorize." off )
if ( CFD_USE_AUTOVEC_FUNCTIONS )
	message( STATUS "Using Auto Vectorization functions." )
	add_compile_definitions( CFD_USE_AUTOVEC_FUNCTIONS )
else()
	message( STATUS "Not using Auto Vectorization functions." )
endif()


# Allow the ability to read fields in (for example to read in an initial condition), this requires VTk library
option ( CFD_HAS_VTK_LIB "Allow extra features which require VTK library to be installed.")
if ( CFD_HAS_VTK_LIB )

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
		add_compile_definitions( CFD_HAS_VTK_LIB )
	else()
		message( FATAL_ERROR "Could not locate VTK library installation! This is required for some features, but not essential!" )
	endif()

else()
	message( STATUS "VTK library will not be used. Some features will be limited." )
endif()


# Add libraries
include_directories(SYSTEM libraries)	#  SYSTEM to suppress warnings

# Eigen library, need to do it this way to not get warnings from it when compiling with Clang
add_library(eigen INTERFACE)
target_sources(eigen INTERFACE libraries/eigen-3.4.0)
target_include_directories(eigen SYSTEM INTERFACE libraries/eigen-3.4.0)

# CGAL library
find_package(CGAL REQUIRED 
			 COMPONENTS Qt5
			 PATHS libraries/CGAL-5.6 NO_DEFAULT_PATH )
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
			 
add_subdirectory(src)

# Link time optimisation
if( supported )
    message(STATUS "IPO / LTO enabled")
    set_property(TARGET main PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif()


