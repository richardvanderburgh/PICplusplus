# CMakeLists.txt

cmake_minimum_required(VERSION 3.5)

project(example_wrapper)

set(BOOST_ROOT "C:/Program\ Files/boost_1_82_0")
set(BOOST_INCLUDE_DIR "C:/Program\ Files/boost_1_82_0/boost")
set(Boost_LIBRARY_DIR "C:/Program\ Files/boost_1_82_0/stage/lib")


# Find Boost and Python libraries
find_package(Boost 1.82.0)
find_package(PythonLibs REQUIRED)
find_package(Matplot++ REQUIRED HINTS $ENV{MATPLOTPP_DIR})

# Define wrapper library
add_library(example_wrapper SHARED pythonBind.cpp)

# Link the wrapper library with the C++ library and Boost libraries
target_link_libraries(example_wrapper PIC++Lib ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

# Set the include directories
target_include_directories(example_wrapper PUBLIC Matplot++::matplot
${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS}

C:/Users/vande/Programming/PICplusplus/lib/include
C:/Users/vande/Programming/PICplusplus/lib/

c:/Program\ Files/matplotplusplus\ 1.1.0/include/matplot
c:/Program\ Files/matplotplusplus\ 1.1.0/include
)
