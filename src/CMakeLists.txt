add_subdirectory(io)
add_subdirectory(fem)

aux_source_directory(. source)
file(GLOB header *.h)

add_library(VEM SHARED ${source})
target_link_libraries(VEM
  Eigen3::Eigen
  spdlog::spdlog
  IO
  FEM
)


target_compile_features(VEM PUBLIC cxx_std_20)
