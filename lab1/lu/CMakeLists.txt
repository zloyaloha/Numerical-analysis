cmake_minimum_required(VERSION 3.10)
project(lu)

find_package(Boost REQUIRED COMPONENTS system filesystem)

include_directories(../../utils/matrix)

add_executable(lu main.cpp)