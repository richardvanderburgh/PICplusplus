cmake_minimum_required(VERSION 3.5)
project(PIC++Lib)

find_package(Matplot++ REQUIRED HINTS $ENV{MATPLOTPP_DIR})
#find_package(Matplot++ REQUIRED HINTS C:/Program\ Files\ (x86)/matplotplusplus\ 1.0.1/lib/cmake/Matplot++)
# Define the library target
add_library(${PROJECT_NAME} STATIC

include/First_EE.hpp
include/init.hpp
include/fft.hpp
include/complex.hpp

fft.cpp
complex.cpp
)

# Set include directories for the library
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
	#${CMAKE_CURRENT_SOURCE_DIR}/deps/include
	c:/Program\ Files/matplotplusplus\ 1.1.0/include/matplot
	c:/Program\ Files/matplotplusplus\ 1.1.0/include
)

# Set compile options for the library
#target_compile_options(my_lib PRIVATE
#    -Wall
#    -Wextra
#    -Werror
#)
target_link_libraries( ${PROJECT_NAME} 
PUBLIC

Matplot++::matplot
)

