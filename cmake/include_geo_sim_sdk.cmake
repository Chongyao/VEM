# include($ENV{HOME}/usr/share/cmake/Modules/geo_sim_sdk.cmake)
set(_GEOSIM_PREFIX_OLD geosim-2.0.0)
set(_GEOSIM_PREFIX_NEW geosim-2.0.1)
set(GEO_SIM_SDK_LIB_DIR $ENV{HOME}/usr/Linux/64/gcc12/${_GEOSIM_PREFIX_NEW}/lib)
set(GEO_SIM_SDK_INCLUDE_DIR $ENV{HOME}/usr/include/${_GEOSIM_PREFIX_NEW} $ENV{HOME}/usr/include)


# include_geo_sim_sdk()

# function(link_geo_sim_include_dir author_name package_name)
#   add_custom_commad(
#     OUTPUT geo_sim_sdk_dir_${author_name}
#     COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/geo_sim_sdk/${author_name}
#     PRE_BUILD
#   )
#   add_custom_command(
#     OUTPUT ${author_name}-${package_name}_include_dir
#     COMMAND ${CMAKE_COMMAND} -E create_symlink ${GEO_SIM_SDK_INCLUDE_DIR}/${author_name}/${package_name} ${PROJECT_BINARY_DIR}/geo_sim_sdk/${author_name}/${package_name}
#     DEPENDS geo_sim_sdk_dir_${author_name}
#     PRE_BUILD
# )
# endfunction(link_geo_sim_include_dir)

# function(add_geo_sim_sdk_package author_name package_name libname)
#   link_geo_sim_include_dir(${author_name} ${package_name})
#   add_library(${libname} INTERFACE)
#   target_link_libraries(${libname} INTERFACE
#     ${GEO_SIM_SDK_LIB_DIR}/lib${libname}.so
#   )
#   target_include_directories(${libname} INTERFACE
#     ${author_name}-${package_name}_include_dir
#   )
# endfunction(add_geo_sim_sdk_package)



# add_geo_sim_sdk_package(zjucad ptree zjucad-ptree)
# add_geo_sim_sdk_package(zcylib io zcylib)

#zcy-io
add_library(zcy-io INTERFACE)
target_link_libraries(zcy-io INTERFACE
  ${GEO_SIM_SDK_LIB_DIR}/libzcy-io.so 
)
target_include_directories(zcy-io INTERFACE
  ${GEO_SIM_SDK_INCLUDE_DIR}
)

add_library(hj-math INTERFACE)
target_link_libraries(hj-math INTERFACE
  ${GEO_SIM_SDK_LIB_DIR}/libhj-math.so
)
target_include_directories(hj-math INTERFACE
  ${GEO_SIM_SDK_INCLUDE_DIR}
)



add_library(zjucad-ptree INTERFACE)
target_link_libraries(zjucad-ptree INTERFACE
  ${GEO_SIM_SDK_LIB_DIR}/libzjucad-ptree.so 
)

target_include_directories(zjucad-ptree INTERFACE
  ${GEO_SIM_SDK_INCLUDE_DIR}
)

# add_library(zjucad-opt-suite INTERFACE)
# target_link_libraries(zjucad-opt-suite INTERFACE
#   ${GEO_SIM_SDK_LIB_DIR}/libzjucad-opt-suite.so 
# )
add_library(jtf-mesh INTERFACE)
target_link_libraries(jtf-mesh INTERFACE
  ${GEO_SIM_SDK_LIB_DIR}/libjtf-mesh.so
)

target_include_directories(jtf-mesh INTERFACE
  ${GEO_SIM_SDK_INCLUDE_DIR}
)


add_library(hj-sparse INTERFACE)
# target_link_libraries(hj-sparse INTERFACE
#   ${GEO_SIM_SDK_LIB_DIR}/libhj-sparse.so
# )
target_include_directories(hj-sparse INTERFACE
  ${GEO_SIM_SDK_INCLUDE_DIR}
)


add_library(zjucad-matrix INTERFACE)
target_include_directories(zjucad-matrix INTERFACE
  ${GEO_SIM_SDK_INCLUDE_DIR}
)


