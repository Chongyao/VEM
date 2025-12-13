find_package(IPOPT REQUIRED MODULE)
if(IPOPT_FOUND)
  message("-- IPOPT @ ${IPOPT_INCLUDE_DIRS}")
  message("-- IPOPT lib @ ${IPOPT_LIBRARIES}")
  # include_directories(${IPOPT_INCLUDE_DIRS})
endif(IPOPT_FOUND)

add_library(Ipopt INTERFACE)
target_link_libraries(Ipopt INTERFACE
  ${IPOPT_LIBRARIES}
)
target_include_directories(Ipopt INTERFACE
  ${IPOPT_INCLUDE_DIRS}
)

# message("CMAKE INSTALL PREFIX = ${PROJECT_BINARY_DIR}")
# CPMAddPackage(
#   NAME Ipopt
#   GIT_REPOSITORY "https://ryon.ren/JinkengLIN/Ipopt-cmake-wrapper.git"
#   GIT_TAG "master"
#   OPTIONS "CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR}"
# )
