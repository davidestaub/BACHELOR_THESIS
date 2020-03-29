cmake_minimum_required(VERSION 3.5)


file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library(simulation ${sources})

target_include_directories(simulation PUBLIC ./)

