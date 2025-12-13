find_package(OpenMP REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
message("CXX compile flags: ${CMAKE_CXX_FLAGS}")
