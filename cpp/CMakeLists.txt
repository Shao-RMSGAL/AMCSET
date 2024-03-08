cmake_minimum_required(VERSION 3.22)
project(eSRIM)

set(CMAKE_CXX_STANDARD 20)

# Create a 'bin' directory for the build
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Set the path to the source directory
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

# Add the source files to the project
add_executable(main ${SOURCE_DIR}/main.cpp
                    ${SOURCE_DIR}/constants.h
                    ${SOURCE_DIR}/eSRIM_classes.h
                    ${SOURCE_DIR}/eSRIM_classes.cpp)

# Include the header files directory
include_directories(${SOURCE_DIR})

# Set the compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# Set the output directory for the executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# Define a CMake option for enabling/disabling debug mode
option(DEBUG_MODE "Enable debug mode" OFF)

# Add compile definitions based on the DEBUG_MODE option
if(DEBUG_MODE)
    add_compile_definitions(DEBUG_MODE)
endif()