cmake_minimum_required(VERSION 3.10)
project(rotating)

# Include directories
include_directories(../../utils/matrix)
include_directories(${Boost_INCLUDE_DIRS})

# Add the executable
add_executable(rotating main.cpp)

# Link Boost libraries to the target
target_link_libraries(rotating ${Boost_LIBRARIES})