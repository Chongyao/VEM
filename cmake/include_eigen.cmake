CPMAddPackage(
  NAME Eigen3
  GIT_REPOSITORY "https://gitlab.com/libeigen/eigen.git"
  GIT_TAG "3.4.0"
  PATCH_COMMAND
  git restore . && git apply "${PROJECT_SOURCE_DIR}/cmake/eigen3_skip_build_demo.patch"
  OPTIONS
  "BUILD_TESTING OFF"
  "EIGEN_BUILD_DOC OFF"
  "EIGEN_BUILD_PKGCONFIG OFF"

)

add_library(Eigen3_WITH_BLAS INTERFACE)
target_link_libraries(Eigen3_WITH_BLAS
  INTERFACE Eigen3::Eigen
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
)
target_compile_definitions(Eigen3_WITH_BLAS
  INTERFACE
  EIGEN_USE_BLAS
  EIGEN_USE_LAPACKE
  INTERFACE
  "UI"
)
