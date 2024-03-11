cmake_minimum_required(VERSION 3.22)
project(eSRIM)

set(CMAKE_CXX_STANDARD 20)

# Create a 'bin' directory for the build
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Set the path to the source directory
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

# Add all source files in the 'src' directory to the project
file(GLOB SOURCE_FILES ${SOURCE_DIR}/*.cpp)

# Add the source files to the project
add_executable(eSRIM ${SOURCE_FILES})

# Include the header files directory
include_directories(${SOURCE_DIR})

# Set the compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pthread")

# Set the output directory for the executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# Define a CMake option for enabling/disabling debug mode
option(DEBUG_MODE "Enable debug mode" OFF)

# Add compile definitions based on the DEBUG_MODE option
if(DEBUG_MODE)
    add_compile_definitions(DEBUG_MODE)

    # Add compiler flags for debugging with gdb
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
else()
    # Add optimization flags if DEBUG_MODE is OFF
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -flto")
endif()
