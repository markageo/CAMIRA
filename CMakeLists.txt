cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(Coupled-Sweep-NS VERSION 1.0)

# Default executable location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# specify the C++ standard
add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_11)

# Helper variables for different compilers
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")

# Add compiler flags
target_compile_options(compiler_flags INTERFACE
    "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wunused;-std=c++2a;-O3>>"
    "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
)


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	message (WARNING "THE BUILD TYPE IS DEBUG. ADDRESS SANITIZER WILL BE USED")
endif()


# Ankiths profiler
option (PROFILING "Enable Profiling Tickers" off)
if (PROFILING)
	message(STATUS "Profiling tickers are on!")
	add_compile_definitions(PROFILING)
	# target_compile_definitions(main PRIVATE PROFILING)
endif()


# include directories
include_directories(${CMAKE_SOURCE_DIR})
include_directories(include)
include_directories(InputParser)

# Add directories
add_subdirectory(fmt)
add_subdirectory(InputParser)
add_subdirectory(src)

