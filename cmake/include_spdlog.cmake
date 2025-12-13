CPMAddPackage(
  NAME spdlog
  GIT_REPOSITORY "https://github.com/gabime/spdlog.git"
  GIT_TAG "v1.13.0"
  OPTIONS
  "SPDLOG_FMT_EXTERNAL ON"
  # PATCH_COMMAND
  # git restore . && git apply "${PROJECT_SOURCE_DIR}/cmake/spdlog_disable_macro_I.patch"

)
