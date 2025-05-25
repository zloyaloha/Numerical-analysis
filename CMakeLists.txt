project(NumericAnalysys)
cmake_minimum_required(VERSION 3.22)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(
  -I/usr/include/c++/11
  -I/usr/include/x86_64-linux-gnu/c++/11
)


include_directories(utils)

add_subdirectory(lab1)
add_subdirectory(lab2)
add_subdirectory(lab3)
add_subdirectory(lab4)