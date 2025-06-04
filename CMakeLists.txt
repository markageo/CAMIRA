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

# BLT
include(blt/SetupBLT.cmake)

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


# Allow the ability to read fields in (for example to read in an initial condition), this requires VTk library
option ( CAMIRA_HAS_VTK_LIB "Allow extra features which require VTK library to be installed.")
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

# Dont build tests and examples
set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)

# Openmp, so that it works with clang
find_package(OpenMP REQUIRED)

# Eigen library
blt_import_library( NAME eigen 
					TREAT_INCLUDES_AS_SYSTEM ON
					INCLUDES libraries/eigen-3.4.0/
							 libraries/eigen-3.4.0/unsupported/
					COMPILE_FLAGS -Wno-array-bounds		# There are a few of these coming from Eigen				
					EXPORTABLE ON )

# CGAL library
find_package(CGAL REQUIRED 
			 COMPONENTS Qt5
			 PATHS libraries/CGAL-5.6 NO_DEFAULT_PATH )
set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)


# Although RAJA and Umpire have their own version of camp in their directories, Umpire seems to have a bug with the 
# camp_DIR variable definition. This gets around it.
add_subdirectory(libraries/camp)
blt_convert_to_system_includes(TARGET camp)

if ( ENABLE_CUDA )
    blt_add_target_definitions(
      TO camp
      SCOPE INTERFACE
      TARGET_DEFINITIONS CAMP_HAVE_CUDA)
endif ()

if(ENABLE_HIP)
    blt_add_target_definitions(
      TO camp
      SCOPE INTERFACE
      TARGET_DEFINITIONS CAMP_HAVE_HIP)
endif()

# Umpire
add_subdirectory(libraries/Umpire)
blt_convert_to_system_includes(TARGET umpire)

# RAJA
add_subdirectory(libraries/RAJA)
blt_convert_to_system_includes(TARGET RAJA)

set(RAJA_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(RAJA_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(RAJA_ENABLE_EXERCISES OFF CACHE BOOL "" FORCE)


set(CMAKE_CUDA_ARCHITECTURES native CACHE STRING "")

add_subdirectory(src)

# Link time optimisation
if( supported )
    message(STATUS "IPO / LTO enabled")
    set_property(TARGET camira PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif()


