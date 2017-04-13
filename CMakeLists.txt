cmake_minimum_required(VERSION 3.7)
project(miniRayTracing)

#set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wextra -Wpedantic")

add_executable(smallpt smallpt.cpp)