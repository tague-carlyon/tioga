cmake_minimum_required(VERSION 3.5)

project(TIOGA CXX)

option(BUILD_SHARED_LIBS "Build shared libraries (default: off)" off)

option(TIOGA_HAS_NODEGID "Support node global IDs (default: on)" off)
option(TIOGA_ENABLE_TIMERS "Track timing information for TIOGA (default: off)" OFF)
option(TIOGA_OUTPUT_STATS "Output statistics for TIOGA holecutting (default: off)" OFF)

option(BUILD_CUDA_SUPPORT "Build with CUDA support (default: on)" on)
option(DEBUG "Build with debug flags (default: off)" on)

find_package(MPI REQUIRED)
# Always build libtioga
add_subdirectory(src)
if(DEBUG)
    message(STATUS "Debug mode enabled")
    set(CMAKE_BUILD_TYPE Debug)
    set(CMAKE_C_FLAGS "-fpic -O0 -g3 -Wall -Wextra -Wpedantic")
    set(CMAKE_CXX_FLAGS "-fpic -O0 -g3 -Wall -Wextra -Wpedantic")
    add_compile_definitions(DEBUG)
else()
    message(STATUS "Release mode enabled")
    SET(CMAKE_C_FLAGS "-O3 -flto -funroll-loops -ftree-vectorize -march=native")
    SET(CMAKE_CXX_FLAGS "-O3 -flto -funroll-loops -ftree-vectorize -march=native")
    set(CMAKE_BUILD_TYPE Release)
endif()

if (APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif()

if (BUILD_SHARED_LIBS)
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
endif()

if (TIOGA_HAS_NODEGID)
    add_definitions(-DTIOGA_HAS_NODEGID)
endif()

if (TIOGA_ENABLE_TIMERS)
    add_definitions(-DTIOGA_ENABLE_TIMERS)
endif()

if (TIOGA_OUTPUT_STATS)
    add_definitions(-DTIOGA_OUTPUT_STATS)
endif()

# CMake installation configuration

install(EXPORT TIOGALibraries
    DESTINATION lib/cmake/TIOGA
    FILE TIOGALibraries.cmake)

# Create TIOGA config so that other codes can find TIOGA
include(CMakePackageConfigHelpers)
set(INCLUDE_INSTALL_DIR include/)
set(LIB_INSTALL_DIR lib/)
configure_package_config_file(
                                cmake/TIOGAConfig.cmake.in
                                ${CMAKE_CURRENT_BINARY_DIR}/TIOGAConfig.cmake
                                INSTALL_DESTINATION lib/cmake/TIOGA
                                PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/TIOGAConfig.cmake
                DESTINATION lib/cmake/TIOGA)
