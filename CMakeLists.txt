cmake_minimum_required(VERSION 3.3)
project(color_detector)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp CIE_matching_functions.h)
add_executable(color_detector ${SOURCE_FILES})