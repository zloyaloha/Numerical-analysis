cmake_minimum_required(VERSION 3.10)
project(boundary-equation)

include_directories(../cauchy-equation)

include_directories(../../utils/matrix)
include_directories(../../lab1/tridiagonal-matrix)

add_executable(diff_boundary main.cpp)