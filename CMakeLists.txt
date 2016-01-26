cmake_minimum_required(VERSION 3.3)
project(color_detector)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(LIBRARY_SOURCE_FILES vecmat/Matrix.cpp vecmat/Vector.cpp vecmat/Utility.cpp)
add_library(core OBJECT ${LIBRARY_SOURCE_FILES})
add_library(static STATIC $<TARGET_OBJECTS:core>)
add_library(shared SHARED $<TARGET_OBJECTS:core>)

set(SOURCE_FILES main.cpp CIE_matching_functions.h)
add_executable(color_detector ${SOURCE_FILES} $<TARGET_OBJECTS:core>)