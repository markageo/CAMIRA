cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(Coupled-Sweep-NS VERSION 1.0)

# Default executable location
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

# specify the C++ standard
add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_20)

# Add compiler flags
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
target_compile_options(compiler_flags INTERFACE
    "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wunused;-Wconversion;-Wpedantic;-std=c++2a;-fconcepts>>"
    "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
)


# Add address sanitizer to the debug build
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	message (WARNING "THE BUILD TYPE IS DEBUG. ADDRESS SANITIZER WILL BE USED")
endif()


# Ankith's profiler
option (PROFILING "Enable Profiling Tickers" off)
if (PROFILING)
	message(STATUS "Profiling tickers are on!")
	add_compile_definitions(PROFILING)
	# target_compile_definitions(main PRIVATE PROFILING)
endif()


# include directories
include_directories(include)
include_directories(libraries)

# Add directories
add_subdirectory(libraries/fmt)
add_subdirectory(libraries/InputParser)
add_subdirectory(src)
add_subdirectory(tests)

