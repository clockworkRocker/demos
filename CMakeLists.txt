cmake_minimum_required(VERSION 3.0)
project(lab1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY lib)

include_directories(inc)

file(GLOB THERMAL2D_SRC src/thermal2d/*.cpp)
add_library(thermal2d STATIC ${THERMAL2D_SRC})

if (USING_OPENMP)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
target_compile_definitions(thermal2d PUBLIC 
  USING_OPENMP,
  OPENMP_THREADS=4
  )
endif()


add_executable(lab1 main.cpp)
target_link_libraries(lab1 thermal2d)