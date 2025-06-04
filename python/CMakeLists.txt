cmake_minimum_required(VERSION 3.12)


# modern cmake implimentation of python bindings
find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)

message(STATUS "Python3 include directories: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python3 libraries: ${Python3_LIBRARIES}")
message(STATUS "numpy include directories: ${Python3_NumPy_INCLUDE_DIRS}")


execute_process(COMMAND python3 -c "import mpi4py; print( mpi4py.get_include() )"
  OUTPUT_VARIABLE mpi4py_inc OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "mpi4py: ${mpi4py_inc}")
message(STATUS "numpy: ${Python3_NumPy_INCLUDE_DIRS}")

#
# Include Directories
#
set(INCLUDE_DIRS ${MPI_C_INCLUDE_PATH}
                 ${Python3_INCLUDE_DIRS}
		 ${mpi4py_inc}
		 ${Python3_NumPy_INCLUDE_DIRS}
		 ${TIOGA_SOURCE_DIR}/src)

file(GLOB PYTIOGA_SRC *.cpp)

if(BUILD_CUDA_SUPPORT)
  list(APPEND PYTIOGA_DEFINITIONS "-DUSE_CUDA")
endif()

message(STATUS "defs: ${PYTIOGA_DEFINITIONS}")

add_library(pytioga SHARED ${PYTIOGA_SRC})
if(BUILD_CUDA_SUPPORT)
  target_link_libraries(pytioga tioga ${MPI_LIBRARIES} ${Python3_LIBRARIES} ${CUDA_LIBRARIES})
else()
  target_link_libraries(pytioga tioga ${MPI_LIBRARIES} ${Python3_LIBRARIES})
endif()
target_include_directories(pytioga PUBLIC ${INCLUDE_DIRS})
target_compile_definitions(pytioga PRIVATE ${PYTIOGA_DEFINITIONS})
set_target_properties(pytioga PROPERTIES
  ARCHIVE_OUTPUT_DIRECOTY ${TIOGA_BINARY_DIR}/lib
  LIBRARY_OUTPUT_DIRECOTY ${TIOGA_BINARY_DIR}/lib
  COMPILE_FLAG "-fPIC"
  PREFIX "" # no lib prefix
  )

