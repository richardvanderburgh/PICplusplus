cmake_minimum_required(VERSION 3.5)
project(PIC++Lib)

find_package(Matplot++ REQUIRED HINTS $ENV{MATPLOTPP_DIR})
# set(FFTW_DIR C:/Users/vande/.conan/data/fftw/3.3.10/_/_/package/9a115689d024a6e2cefc8c2de3eabb39c2da0603/lib)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../build)

find_package(FFTW3 REQUIRED)
find_package(nlohmann_json REQUIRED)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

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
	#c:/Program\ Files/matplotplusplus\ 1.1.0/include/matplot
	#c:/Program\ Files/matplotplusplus\ 1.1.0/include
	
	# ${EIGEN_INCLUDE_DIRS}
	${CONAN_INCLUDE_DIRS}
	
)

target_link_libraries( ${PROJECT_NAME} 

Matplot++::matplot

 ${CONAN_LIBS}
)

