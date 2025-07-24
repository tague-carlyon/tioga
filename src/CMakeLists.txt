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
  # Explicitly enable CUDA language support
  enable_language(CUDA)

  # Manually specify CUDA architectures
  set(CMAKE_CUDA_ARCHITECTURES 80)

  # Automatically set CMAKE_CUDA_ARCHITECTURES based on the current PC configuration
  include(FindCUDA/select_compute_arch)

  find_package(CUDAToolkit REQUIRED)

  # Basic CUDA settings
  set(CUDA_SEPARABLE_COMPILATION ON CACHE INTERNAL "Compile GPU code separately")
  set(CUDA_PROPAGATE_HOST_FLAGS OFF CACHE INTERNAL "Do not propagate host flags")
  list(APPEND CUDA_OPTS "-O3" "--compile -march=native")
  list(APPEND TIOGA_DEFINITIONS "-DUSE_CUDA")
  list(APPEND CUDA_OPTS --compiler-options -fPIC)

  # Get all the CUDA files
  file(GLOB cu_sources "*.cu")
  message(STATUS "CUDA sources: ${cu_sources}")

  # Create a separate library for CUDA files
  foreach(cu_file ${cu_sources})
    get_filename_component(cu_file_name ${cu_file} NAME_WE)
    add_library(${cu_file_name}_obj OBJECT ${cu_file})
    target_compile_options(${cu_file_name}_obj PRIVATE ${CUDA_OPTS})
    target_include_directories(${cu_file_name}_obj PUBLIC ${MPI_C_INCLUDE_PATH})
    list(APPEND cuda_objects $<TARGET_OBJECTS:${cu_file_name}_obj>)
  endforeach()

  # Combine CUDA objects into the main library
  add_library(tioga STATIC ${TIOGA_SOURCES} ${cuda_objects})
  target_compile_definitions(tioga PRIVATE ${TIOGA_DEFINITIONS})
  target_include_directories(tioga PUBLIC ${MPI_C_INCLUDE_PATH})
  target_link_libraries(tioga ${MPI_LIBRARIES} ${CMAKE_DL_LIBS} ${CUDA_LIBRARIES} CUDA::cudart)
else()
  message(STATUS "BUILDING WITHOUT CUDA SUPPORT")
  add_library(tioga STATIC ${TIOGA_SOURCES})
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
  RUNTIME DESTINATION bin/tioga
  ARCHIVE DESTINATION lib/tioga
  LIBRARY DESTINATION lib/tioga)

install(FILES ${TIOGA_HEADERS}
  DESTINATION include/tioga)
