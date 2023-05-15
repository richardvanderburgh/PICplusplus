cmake_minimum_required(VERSION 3.5)
project(PIC++Main)

add_subdirectory(lib)

add_executable(${CMAKE_PROJECT_NAME} 
src/main.cpp
)

target_link_libraries( ${CMAKE_PROJECT_NAME} 
PUBLIC

PIC++Lib
)

