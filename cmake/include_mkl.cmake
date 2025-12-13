find_package(MKL CONFIG REQUIRED)

if(MKL_FOUND)
  message(STATUS "Found oneMKL")
  # MKL::mkl 目标包含了头文件路径，但如果需要，也可以手动添加
  # target_include_directories(your_executable PRIVATE ${MKL_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "oneMKL not found!")
endif()



