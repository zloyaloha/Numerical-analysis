cmake_minimum_required(VERSION 3.10)
project(hyperbolic2d-equation)

include_directories(../utils/matrix)
include_directories(../lab1/tridiagonal-matrix)

add_executable(hyperbolic2d main.cpp)