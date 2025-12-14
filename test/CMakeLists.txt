set(TESTS
    geo
    vem
    solver
    patch_linear
    patch_single
    fem
    vem_vs_fem
    beam_gravity
    stability
    adjacent
    agglomeration_preview
)

foreach(TEST IN LISTS TESTS)
    message(${TEST})
    add_executable(test_${TEST} test_${TEST}.cpp)
    target_link_libraries(test_${TEST} PRIVATE VEM)
endforeach()
