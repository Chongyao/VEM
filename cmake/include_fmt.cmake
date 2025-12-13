CPMAddPackage(
  NAME fmt
  GIT_REPOSITORY "https://github.com/fmtlib/fmt"
  GIT_TAG "9.0.0"
  PATCH_COMMAND
  git restore . && git apply "${PROJECT_SOURCE_DIR}/cmake/fmt_disable_macro_I.patch"

)


# add_library(fmt_FixLapak INTERFACE)
# target_link_libraries(fmt_FixLapak
#   INTERFACE fmt::fmt
# )
# target_compile_definitions(fmt_FixLapak
#   PRIVATE
#   "-UI"
# )
