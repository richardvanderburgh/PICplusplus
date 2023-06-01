cmake_minimum_required(VERSION 3.5)
project(PIC++Lib)

find_package(Matplot++ REQUIRED HINTS $ENV{MATPLOTPP_DIR})
add_library(${PROJECT_NAME} STATIC

include/Accel.hpp
include/First_EE.hpp
include/Fields.hpp
include/init.hpp
include/fft.hpp
include/complex.hpp
include/SetRho.hpp

fft.cpp
complex.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries( ${PROJECT_NAME} 

Matplot++::matplot
)

