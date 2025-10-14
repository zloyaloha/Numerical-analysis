cmake_minimum_required(VERSION 3.10)
project(tridiagonal)

find_package(Boost REQUIRED COMPONENTS system filesystem)

include_directories(../../utils/matrix)

add_executable(tridiagonal main.cpp)