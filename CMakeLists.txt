# cmake_minimum_required(VERSION 3.5)
# project(PIC++Main)

# add_subdirectory(lib)
# add_subdirectory(test)

# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/build)

# include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
# conan_basic_setup()
# find_package(gtest REQUIRED)
# # find_package(boost REQUIRED)
# find_package(FFTW3 REQUIRED)
# find_package(nlohmann_json REQUIRED)

# add_executable(${CMAKE_PROJECT_NAME} 
# src/main.cpp
# )

# target_include_directories(${PROJECT_NAME} PUBLIC
	# ${CONAN_INCLUDE_DIRS}
# )

# target_link_libraries( ${CMAKE_PROJECT_NAME} 
# PUBLIC
 # ${CONAN_LIBS}

# PIC++Lib
# )

# project(PICplusplusTest)

# enable_testing()

# add_executable(${PROJECT_NAME} 
# test/PICplusplusTest.cpp
# test/include/PICTest.hpp
# )

# target_include_directories(${PROJECT_NAME} PUBLIC
	# ${CONAN_INCLUDE_DIRS}
# )

# target_link_libraries( ${PROJECT_NAME} 
# ${CONAN_LIBS}
# PIC++Lib
# )

# include(GoogleTest)
# gtest_discover_tests(${PROJECT_NAME} )
cmake_minimum_required(VERSION 3.5)
project(PIC++Main)

add_subdirectory(lib)
add_subdirectory(test)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/build)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

find_package(gtest REQUIRED)
find_package(FFTW3 REQUIRED)
find_package(nlohmann_json REQUIRED)

add_executable(${CMAKE_PROJECT_NAME}
    src/main.cpp
)

target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC
    ${CONAN_INCLUDE_DIRS}
)

target_link_libraries(${CMAKE_PROJECT_NAME}
    PUBLIC
    ${CONAN_LIBS}
    PIC++Lib
)

project(PICplusplusTest)

enable_testing()

add_executable(${PROJECT_NAME}
    test/PICplusplusTest.cpp
    test/include/PICTest.hpp
)

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CONAN_INCLUDE_DIRS}
)

link_directories(C:/Users/vande/.conan/data/gtest/1.13.0/_/_/package/10f3c0a3a24a1055f3670e6694da4a93d078e663/lib)

target_link_libraries(${PROJECT_NAME}
    PRIVATE
    ${CONAN_LIBS}
    PIC++Lib
    gtest
)

include(GoogleTest)
gtest_discover_tests(${PROJECT_NAME})

