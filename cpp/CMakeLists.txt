cmake_minimum_required(VERSION 3.22)
project(eSRIM)

set(CMAKE_CXX_STANDARD 17)

# Add the source files to the project
add_executable(main main.cpp energy_functions.cpp)

# Include the header files directory
include_directories(${CMAKE_SOURCE_DIR})

# Set the compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")