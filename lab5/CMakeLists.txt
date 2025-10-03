cmake_minimum_required(VERSION 3.10)
project(parabolic-equation)

include_directories(../utils/matrix)
include_directories(../lab1/tridiagonal-matrix)

add_executable(parabolic main.cpp)