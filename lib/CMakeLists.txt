cmake_minimum_required(VERSION 3.5)
project(PIC++Lib)

# Define the library target
add_library(${PROJECT_NAME} STATIC

include/First_EE.hpp
include/init.hpp
include/fft.h
include/complex.h

fft.cpp
complex.cpp
)

# Set include directories for the library
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Set compile options for the library
#target_compile_options(my_lib PRIVATE
#    -Wall
#    -Wextra
#    -Werror
#)
