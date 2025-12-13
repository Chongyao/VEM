CPMAddPackage(
  NAME spectra
  GIT_REPOSITORY "https://github.com/yixuan/spectra.git"
  GIT_TAG "v1.0.1"
  PATCH_COMMAND
  git restore . && git apply "${PROJECT_SOURCE_DIR}/cmake/spectra_fix_eigen.patch"
  # git restore . && git apply "${PROJECT_SOURCE_DIR}/cmake/suitesparse_expose_metis.patch"

)
