cmake_minimum_required(VERSION 3.10)
project(cube-splines)

include_directories(../../utils/matrix)
include_directories(../../lab1/tridiagonal-matrix)

add_executable(splines main.cpp)