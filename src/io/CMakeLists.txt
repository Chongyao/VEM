aux_source_directory(. source)
file(GLOB header *.h)

add_library(IO SHARED ${source})
target_link_libraries(IO
  Eigen3::Eigen
  spdlog::spdlog
  tinyxml2::tinyxml2
)


target_compile_features(IO PUBLIC cxx_std_20)
