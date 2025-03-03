cmake_minimum_required(VERSION 3.10)
project(lab1)

find_package(Boost REQUIRED COMPONENTS system filesystem)

include_directories(lu)
include_directories(../utils/matrix)

add_executable(main main.cpp)