cmake_minimum_required(VERSION 2.7)
project(apc_omp CXX)

#Include Eigen directory
include_directories(libEigen)

set(ENABLE_OMP On CACHE BOOL "Enable OpenMP")

if(ENABLE_OMP)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

# Note: -fopenmp automatically links libgomp, here we manually force its linking
#       so that our programs keep working even if OpenMP is disabled.

#Set flag optimization O3
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_definitions(-std=c++0x)

add_executable(apc_omp main_omp.cpp spmatfullomp.cpp spmatvwbbsomp.cpp)
target_link_libraries(apc_omp gomp)
