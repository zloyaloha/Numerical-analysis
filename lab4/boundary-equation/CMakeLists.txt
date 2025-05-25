cmake_minimum_required(VERSION 3.10)
project(boundary-equation)

include_directories(../cauchy-equation)
add_executable(diff_boundary main.cpp)