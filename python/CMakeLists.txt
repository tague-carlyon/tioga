if(NOT ON_ZARATAN)
  message(STATUS "Building on idea system")
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
else()
  message(STATUS "Building on Zaratan system")
  find_package(Python3 3.10 EXACT REQUIRED COMPONENTS Interpreter Development NumPy)
endif()

execute_process(COMMAND python3 -c "import mpi4py; print( mpi4py.get_include() )"
   OUTPUT_VARIABLE mpi4py_inc OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "mpi4py: ${mpi4py_inc}")
message(STATUS "numpy: ${Python3_NumPy_INCLUDE_DIRS}")
message(STATUS "Python Include Dir: ${Python3_INCLUDE_DIRS}")
message(STATUS "Python Libraries: ${Python3_LIBRARIES}")

#
# Include Directories
#
set(INCLUDE_DIRS ${MPI_C_INCLUDE_PATH}
                 ${Python3_INCLUDE_DIRS}
                 ${mpi4py_inc}
                 ${Python3_NumPy_INCLUDE_DIRS}
                 ${TIOGA_SOURCE_DIR}/src)

file(GLOB PYTIOGA_SRC *.cpp)
message(STATUS "PYTIOGA_SRC: ${PYTIOGA_SRC}")

if(BUILD_CUDA_SUPPORT)
  list(APPEND PYTIOGA_DEFINITIONS "-DUSE_CUDA")
endif()

message(STATUS "defs: ${PYTIOGA_DEFINITIONS}")

# Create the pytioga library
add_library(pytioga SHARED ${PYTIOGA_SRC})

# Set properties for the pytioga target
target_include_directories(pytioga PUBLIC ${INCLUDE_DIRS})
target_compile_definitions(pytioga PRIVATE ${PYTIOGA_DEFINITIONS})

# Link libraries to the pytioga target
if(BUILD_CUDA_SUPPORT)
  target_link_libraries(pytioga tioga ${MPI_LIBRARIES} ${Python3_LIBRARIES} ${CUDA_LIBRARIES})
else()
  target_link_libraries(pytioga tioga ${MPI_LIBRARIES} ${Python3_LIBRARIES})
endif()

set_target_properties(pytioga PROPERTIES
  OUTPUT_NAME "pytioga"  # Ensures the shared library is named pytioga.so
  PREFIX ""              # Removes the default "lib" prefix
  SUFFIX ".so"           # Ensures the correct suffix for Python modules
  ARCHIVE_OUTPUT_DIRECTORY ${TIOGA_BINARY_DIR}/python
  LIBRARY_OUTPUT_DIRECTORY ${TIOGA_BINARY_DIR}/python
)

install(TARGETS pytioga
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.x/site-packages
)