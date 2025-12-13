CPMAddPackage(
  NAME libigl
  GIT_REPOSITORY "https://github.com/libigl/libigl.git"
  GIT_TAG "v2.5.0"

  PATCH_COMMAND
  git restore . && git apply "${PROJECT_SOURCE_DIR}/cmake/libigl_disable_macro_I.patch"

)
