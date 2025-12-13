project(arborescence)
cmake_minimum_required(VERSION 3.1.0)
cmake_policy(SET CMP0169 OLD)

set(CPM_USE_NAMED_CACHE_DIRECTORIES ON)

set(CPM_SOURCE_CACHE ${CPM_SOURCE_CACHE}/${CMAKE_BUILD_TYPE}-${CMAKE_C_COMPILER_ID})
include(cmake/CPM.cmake)
include(cmake/code-format.cmake)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/;${CMAKE_MODULE_PATH}")


unset(CMAKE_CUDA_COMPILER CACHE)


option(ENABLE_MEMORY_TEST "Enable memory test" OFF)
if(${ENABLE_MEMORY_TEST})
  message("ENABLE_MEMORY_TEST")
  add_compile_options("-fsanitize=address")
  add_link_options("-fsanitize=address")
endif(${ENABLE_MEMORY_TEST})


set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} -fdefault-real-8")

#STATIC or NOT
set(ZCY_COMPILE_AS_LIBRARY ON)
if(ZCY_COMPILE_AS_LIBRARY)
  add_compile_definitions(ZCY_STATIC_LIBRARY)
endif()


#O3 optimization
message("cmake build type: ${CMAKE_BUILD_TYPE}")
if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")
endif()
#enable all instructions
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wno-c++11-narrowing -Wno-unused-variable -Wno-unused-parameter -Wno-array-bounds -Wno-dangling-reference -Wno-tautological-compare")
set(DISABLED_WARNING "-Wno-unused-but-set-variable  -Wno-unused-parameter -Wno-c++11-narrowing -Wno-unused-variable  -Wno-array-bounds -Wno-tautological-compare -Wno-deprecated-literal-operator")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native ${DISABLED_WARNING}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -Wno-everything")



# PIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#maxima
set(MAXIMA_OUTPUT_DIR "${PROJECT_SOURCE_DIR}/maxima")


# #libigl
# include(libigl)
# igl_include(glfw)



#set(CMAKE_CUDA_ARCHITECTURES 75)
# List of all includes
set(INCLUDES
  # suitesparse
  # cereal
  # openmp
  # backward
  # boost
  eigen
  # spectra
  # blas
  # amgcl
  # progress_bar
  # fmt
  # spdlog
  # ipopt
  # # vtkio
  # autodiff
  # cgal
  # fast_matrix_market
  # nlopt
  # OpenVolumeMesh
  # libigl
  # tbb
  # mkl
  # mpi
)

foreach(INCLUDE IN LISTS INCLUDES)
  include(cmake/include_${INCLUDE}.cmake)
  message("Load ${INCLUDE} done.")
endforeach()













# add_subdirectory(3rd)
include_directories(${PROJECT_SOURCE_DIR}/src)
add_subdirectory(src)
# add_subdirectory(exec)
add_subdirectory(test)
# add_subdirectory(apps)


