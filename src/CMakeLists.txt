

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
  linCartInterp.C
  parallelComm.C
  search.C
  searchADTrecursion.C
  tioga.C
  tiogaInterface.C
  )

file(GLOB TIOGA_HEADERS *.h)

if(BUILD_CUDA_SUPPORT)
  # Basic CUDA settings
  set(CUDA_SEPARABLE_COMPILATION ON CACHE INTERNAL "Compile GPU code seperately")
  set(CUDA_PROPAGATE_HOST_FLAGS OFF CACHE INTERNAL "Do not propagate host flags")
  list(APPEND CUDA_OPTS "-O2")
  list(APPEND TIOGA_DEFINITIONS "-DUSE_CUDA")
  list(APPEND CUDA_OPTS ${CMAKE_C_COMPILER} -fPIC)
  # -ccbin 
  # Get all the CUDA files
  file(GLOB cu_sources "*.cu")
  message(STATUS "sources: ${cu_sources}")

  add_library(tioga STATIC ${TIOGA_SOURCES} ${cu_sources})
  target_include_directories(tioga PUBLIC ${MPI_C_INCLUDE_PATH})

  target_compile_options(tioga PRIVATE ${CUDA_OPTS})
  set_target_properties(tioga PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${TIOGA_BINARY_DIR}/src
    LIBRARY_OUTPUT_DIRECTORY ${TIOGA_BINARY_DIR}/src
    )
  target_compile_definitions(tioga PRIVATE ${TIOGA_DEFINITIONS})
  target_include_directories(tioga PUBLIC ${MPI_C_INCLUDE_PATH})
  set(CUDA_LIBRARIES ${CUDA_CUDART_LIBRARY} ${CUDA_CUBLAS_LIBRARIES} ${CUDA_CURAND_LIBRARIES})
  message(STATUS "cudalibs: ${CUDA_LIBRARIES}")
  #
  # Static libraries dont actually link to external libraries... BUT
  # if you tell cmake it does, cmake will remember and add those
  # libraries when the ".a" file is linked somewhere else. If the ".a"
  # file is linked externally (like form another project), just make
  # sure the shared object or executable also links these libraries:
  #
  target_link_libraries(tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS} ${CUDA_LIBRARIES})

else()
  message(STATUS "BUILDING WITHOUT CUDA SUPPORT")
  # add_library(tioga SHARED ${TIOGA_SOURCES})
  add_library(tioga ${TIOGA_SOURCES})
  target_compile_definitions(tioga PRIVATE ${TIOGA_DEFINITIONS})
  target_include_directories(tioga PUBLIC ${MPI_C_INCLUDE_PATH})
  target_link_libraries(tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS})

endif()


if(MPI_COMPILE_FLAGS)
  set_target_properties(tioga PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(tioga PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

install(TARGETS tioga
  EXPORT "${CMAKE_PROJECT_NAME}Libraries"
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(FILES ${TIOGA_HEADERS}
  DESTINATION include)
