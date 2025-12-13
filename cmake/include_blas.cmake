set(BLA_VENDOR OpenBLAS)
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  message("-- BLAS libs @ ${BLAS_LIBRARIES}")
endif(BLAS_FOUND)


find_package(LAPACK REQUIRED)

if(LAPACK_FOUND)
  message(STATUS "LAPACK found.")
  message(STATUS "LAPACK include directories: ${LAPACK_INCLUDE_DIRS}")
  message(STATUS "LAPACK libraries: ${LAPACK_LIBRARIES}")
else()
  message(FATAL_ERROR "LAPACKE not found.")
endif()
