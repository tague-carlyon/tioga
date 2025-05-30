cmake_minimum_required(VERSION 2.8)

project(TIOGA CXX C Fortran)

option(BUILD_SHARED_LIBS "Build shared libraries (default: off)" off)
option(BUILD_TIOGA_EXE "Build tioga driver code (default: off)" off)
option(BUILD_GRIDGEN_EXE "Build grid generator code (default: off)" off)
option(BUILD_PYTHON "Build Python wrapper" on)
option(BUILD_CUDA_SUPPORT "Build with support for GPU memory exchange" on)

#
# Find Libraries
#
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

set(MPI_CXX_COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS} -fpic -O3 -flto -funroll-loops -ftree-vectorize -march=native")
set(MPI_C_COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS} -fpic -O3 -flto -funroll-loops -ftree-vectorize -march=native")

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MPI_C_COMPILE_FLAGS} ${C_DEFAULT_FLAGS} -fpic -O3 -flto -funroll-loops -ftree-vectorize -march=native")
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS} ${Fortran_DEFAULT_FLAGS} -Ofast -O5")
SET(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILEFLAGS} ${CXX_DEFAULT_FLAGS} -fpic -O3 -flto -funroll-loops -ftree-vectorize -march=native")

# Set some default compilation settings for Fortran compiler
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -ffree-line-length-none -fPIC")
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -double_size 128 -fPIC")
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
