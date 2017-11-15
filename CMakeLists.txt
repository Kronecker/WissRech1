cmake_minimum_required(VERSION 3.6)
project(WissRech1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp jacobiIterationSolver.cpp jacobiIterationSolver.h)
add_executable(WissRech1 ${SOURCE_FILES})