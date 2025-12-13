add_executable(test_geo test_geo.cpp)
target_link_libraries(test_geo VEM)


# test_vem
add_executable(test_vem test_vem.cpp)
target_link_libraries(test_vem PRIVATE VEM)