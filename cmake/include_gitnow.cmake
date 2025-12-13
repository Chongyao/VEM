find_program(MAKE_BIN make)
find_program(INSTALL_BIN install)
if (MAKE_BIN AND INSTALL_BIN)
  file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/git-now)
  execute_process(
    COMMAND ${MAKE_BIN} prefix=${PROJECT_BINARY_DIR}/git-now install
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/3rd/git-now
  )
endif()
