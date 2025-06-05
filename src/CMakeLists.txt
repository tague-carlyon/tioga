

set(TIOGA_SOURCES
  # Fortran sources
  kaiser.f
  cellVolume.f90
  median.F90

  # C sources
  buildADTrecursion.c
  linklist.c
  math.c
  get_amr_index_xyz.c
  utils.c

  # CXX sources
  ADT.C
  CartBlock.C
  CartGrid.C
  MeshBlock.C
  bookKeeping.C
  cartOps.C
  checkContainment.C
  dataUpdate.C
  dataUpdateGPU.C
  exchangeAMRDonors.C
  exchangeBoxes.C
  exchangeDonors.C
  exchangeSearchData.C
  getCartReceptors.C
  highOrder.C
  holeMap.C
  parallelComm.C
  search.C
  searchADTrecursion.C
  tioga.C
  tiogaInterface.C
  )

file(GLOB TIOGA_HEADERS *.h)

if(BUILD_CUDA_SUPPORT)
  enable_language(CUDA)
  set(CMAKE_CUDA_STANDARD 14)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)
  set(CMAKE_CUDA_ARCHITECTURES 80)
  message(STATUS "CMAKE_CUDA_ARCHITECTURES: ${CMAKE_CUDA_ARCHITECTURES}")  # Load the CUDA cmake variables
  find_package(CUDAToolkit REQUIRED)
  message(STATUS "CUDAToolkit include directories: ${CUDAToolkit_INCLUDE_DIRS}")
  message(STATUS "CUDAToolkit libraries: ${CUDAToolkit_LIBRARIES}")

  # Basic CUDA settings
  set(CUDA_SEPARABLE_COMPILATION ON CACHE INTERNAL "Compile GPU code seperately")
  set(CUDA_PROPAGATE_HOST_FLAGS OFF CACHE INTERNAL "Do not propagate host flags")
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --compiler-options -O3 -fPIC")
  list(APPEND TIOGA_DEFINITIONS "-DUSE_CUDA")
  #list(APPEND CUDA_OPTS -ccbin ${CMAKE_C_COMPILER} --compiler-options -fPIC)

  # Get all the CUDA files
  file(GLOB cu_sources "*.cu")
  message(STATUS "cuda sources: ${cu_sources}")
  #cuda_include_directories( ${MPI_C_INCLUDE_PATH} )

  add_library(tioga SHARED ${TIOGA_SOURCES} ${cu_sources})
  target_link_libraries(tioga CUDA::cudart)
  #target_compile_options(tioga PRIVATE ${CUDA_OPTS})
  target_compile_definitions(tioga PRIVATE ${TIOGA_DEFINITIONS})
  set_target_properties(tioga PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON "Compile GPU code seperately"
    CUDA_PROPAGATE_HOST_FLAGS OFF "Do not propagate host flags"
    ARCHIVE_OUTPUT_DIRECTORY ${TIOGA_BINARY_DIR}/src
    LIBRARY_OUTPUT_DIRECTORY ${TIOGA_BINARY_DIR}/src
    )
  target_compile_definitions(tioga PRIVATE ${TIOGA_DEFINITIONS})
  target_include_directories(tioga PUBLIC ${MPI_C_INCLUDE_PATH})

  #
  # Static libraries dont actually link to external libraries... BUT
  # if you tell cmake it does, cmake will remember and add those
  # libraries when the ".a" file is linked somewhere else. If the ".a"
  # file is linked externally (like form another project), just make
  # sure the shared object or executable also links these libraries:
  #
  target_link_libraries(tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})

else()

  # add_library(tioga SHARED ${TIOGA_SOURCES})
  add_library(tioga ${TIOGA_SOURCES})
  target_compile_definitions(tioga PRIVATE ${TIOGA_DEFINITIONS})
  target_include_directories(tioga PUBLIC ${MPI_C_INCLUDE_PATH})
  target_link_libraries(tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})

endif()



# set_target_properties(tioga PROPERTIES COMPILE_FLAGS "-std=c++11 -fPIC")

if(MPI_COMPILE_FLAGS)
  set_target_properties(tioga PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(tioga PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

# if(BUILD_CUDA_SUPPORT)
#   install(TARGETS tioga gpufunc
#     EXPORT "${CMAKE_PROJECT_NAME}Libraries"
#     RUNTIME DESTINATION bin
#     ARCHIVE DESTINATION lib
#     LIBRARY DESTINATION lib)
# else()
  install(TARGETS tioga
    EXPORT "${CMAKE_PROJECT_NAME}Libraries"
    RUNTIME DESTINATION bin
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
# endif()

install(FILES ${TIOGA_HEADERS}
  DESTINATION include)
