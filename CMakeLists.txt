cmake_minimum_required(VERSION 3.6)
project(WissRech1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -msse4 -O3")

set(SOURCE_FILES main.cpp jacobiIterationSolver.cpp jacobiIterationSolver.h solveTask3.cpp solveTasks.h inverseIterationMethods.cpp solveTask0.cpp mUtil.cpp mUtil.h solveTask4.cpp solveTask5.cpp executionData.h Abgabe/Task6_7_8/solveTask6.cpp Abgabe/Task6_7_8/solveTask7.cpp Abgabe/Task6_7_8/solveTask8.cpp solveTask9.cpp memAlignOS.cpp memAlignOS.h)
add_executable(WissRech1 ${SOURCE_FILES})