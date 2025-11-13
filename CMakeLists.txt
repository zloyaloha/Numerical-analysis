project(NumericAnalysys)
cmake_minimum_required(VERSION 3.22)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(utils)

add_subdirectory(lab1)
add_subdirectory(lab2)
add_subdirectory(lab3)
add_subdirectory(lab4)
add_subdirectory(lab5)
add_subdirectory(lab6)