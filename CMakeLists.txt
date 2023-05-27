cmake_minimum_required(VERSION 3.5)
project(PIC++Main)

add_subdirectory(lib)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/build)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()
find_package(boost REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(nlohmann_json REQUIRED)

add_executable(${CMAKE_PROJECT_NAME} 
src/main.cpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
	${CONAN_INCLUDE_DIRS}
)

target_link_libraries( ${CMAKE_PROJECT_NAME} 
PUBLIC
 ${CONAN_LIBS}

PIC++Lib
)