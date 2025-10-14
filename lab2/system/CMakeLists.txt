cmake_minimum_required(VERSION 3.10)
project(system)

include_directories(../../utils/matrix)
include_directories(../../lab1/lu)

add_executable(system main.cpp)