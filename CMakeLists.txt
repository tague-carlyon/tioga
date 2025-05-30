cmake_minimum_required(VERSION 3.5)

project(TIOGA CXX C Fortran)

# Automatically set CMAKE_CUDA_ARCHITECTURES based on the current PC configuration
include(FindCUDA/select_compute_arch)
# Manually specify CUDA architectures
set(CMAKE_CUDA_ARCHITECTURES 80 90)

# Explicitly enable CUDA language support
enable_language(CUDA)

option(BUILD_SHARED_LIBS "Build shared libraries (default: off)" off)
option(BUILD_TIOGA_EXE "Build tioga driver code (default: off)" off)
option(BUILD_GRIDGEN_EXE "Build grid generator code (default: off)" off)
option(TIOGA_HAS_NODEGID "Support node global IDs (default: on)" off)
option(TIOGA_ENABLE_TIMERS "Track timing information for TIOGA (default: off)" OFF)
option(TIOGA_OUTPUT_STATS "Output statistics for TIOGA holecutting (default: off)" OFF)
option(BUILD_PYTHON "Build python wrapper (default: on)" on)
option(BUILD_CUDA_SUPPORT "Build with CUDA support (default: on)" on)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Ensure C++11 standard is enabled
if (CMAKE_VERSION VERSION_LESS "3.1")
  set(CMAKE_CXX_FLAGS "-g -rdynamic -fPIC -std=c++0x ${CMAKE_CXX_FLAGS}")
else()
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_FLAGS "-O3 -fPIC ${CMAKE_CXX_FLAGS}")
endif()

# Set some default compilation settings for Fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -fbounds-check -fbacktrace -fdefault-real-8 -ffree-line-length-none -fPIC")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -double_size 128 -fPIC")
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

# Always build libtioga
add_subdirectory(src)

# Optionally build driver exe and gridGen if the user requests it
if (BUILD_TIOGA_EXE)
  add_subdirectory(driver)
endif()

if (BUILD_GRIDGEN_EXE)
  add_subdirectory(gridGen)
endif()

if (BUILD_PYTHON)
  add_subdirectory(python)
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
