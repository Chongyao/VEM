CPMAddPackage(
  NAME progressbar
  GIT_REPOSITORY https://github.com/gipert/progressbar.git
  GIT_TAG v2.1
)
add_library(progressbar INTERFACE)
target_include_directories(progressbar INTERFACE
  ${progressbar_SOURCE_DIR}/include
)
