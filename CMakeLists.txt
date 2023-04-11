cmake_minimum_required(VERSION 3.21)
project(lab2_v2)

set(CMAKE_CXX_STANDARD 11)

# -------- Directories --------
set(${PROJECT_NAME}_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
set(${PROJECT_NAME}_SOURCES_DIR ${PROJECT_SOURCE_DIR}/src)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

# -------- Include Eigen --------
find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIRS})

# -------- Make executable --------
include_directories(${${PROJECT_NAME}_INCLUDE_DIR})
# file(GLOB_RECURSE ${PROJECT_NAME}_FILES ${${PROJECT_NAME}_SOURCES_DIR}/solver.cpp)
add_executable(lab2_v2 main.cpp ${${PROJECT_NAME}_FILES})
