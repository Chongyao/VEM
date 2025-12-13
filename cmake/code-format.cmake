include_guard()
###############################################################################
# ClangFormat
###############################################################################
find_program(CLANG_FORMAT clang-format)

if(CLANG_FORMAT)
  message(STATUS "clang-format find: ${CLANG_FORMAT}")
else()
  message(FATAL_ERROR "clang-format not found! please set {CLANG_FORMAT_ROOT} correctly!")
endif()

function(code_format)
  file(GLOB_RECURSE FORMAT_SOURCES ${ARGN})
  add_custom_target(code-format-inplace
    COMMAND ${CLANG_FORMAT} -i ${FORMAT_SOURCES}
  )
  add_custom_target(code-format
    COMMAND ${CLANG_FORMAT} --dry-run ${FORMAT_SOURCES}
  )
endfunction()

set(source_path
  src
  examples
  apps
  tests
)

foreach(path ${source_path})
  set(source_pattern
    ${source_pattern}
    ${path}/*.cc ${path}/*.hh
    ${path}/*.[ch]pp
    ${path}/*.[ch]xx
    ${path}/*.[ch]
  )
endforeach()

message(STATUS "[chaos/tools]: Enable clang-format, add target: `code-format`")
code_format(${source_pattern})

unset(source_pattern)
unset(source_path)
