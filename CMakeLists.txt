cmake_minimum_required(VERSION 3.5)
project(PIC++Main)

# Find Eigen package
#find_package(Eigen3 REQUIRED)

add_subdirectory(deps/matplotplusplus)
add_subdirectory(lib)

# Add executable target
add_executable(${CMAKE_PROJECT_NAME} 
src/main.cpp
)

# Link against Eigen library
target_link_libraries( ${CMAKE_PROJECT_NAME} 

PIC++Lib
#Eigen3::Eigen 

matplot)

