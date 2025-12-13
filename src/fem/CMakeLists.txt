aux_source_directory(. source)
file(GLOB header *.h)

add_library(FEM SHARED ${source})
target_link_libraries(FEM
  Eigen3::Eigen
)


target_compile_features(FEM PUBLIC cxx_std_20)
