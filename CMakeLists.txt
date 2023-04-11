project(IWO_test)
cmake_minimum_required(VERSION 3.0)

include_directories(${PROJECT_SOURCE_DIR}/include)
find_package(Threads)

add_library(iwo 
  src/individual.cpp
  src/population.cpp
  src/optimizer.cpp
)

add_executable(test_normal
  test/test_normal.cpp
)

target_link_libraries(test_normal
  iwo
  ${CMAKE_THREAD_LIBS_INIT}
)