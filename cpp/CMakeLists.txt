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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -pthread")

# Set the output directory for the executable
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# Define a CMake option for enabling/disabling debug mode
option(DEBUG_MODE "Enable debug mode" OFF)

# Define a CMake option for enabling/disabling performance analysis mode
option(PERF "Enable performance analysis mode" OFF)

# Define a CMake option for enabling/disabling production mode
option(PROD "Enable production mode" ON)

# Add compile definitions based on the DEBUG_MODE option
if(DEBUG_MODE)
    message("Debug mode is enabled.")
    add_compile_definitions(DEBUG_MODE)
    # Add compiler flags for debugging with gdb
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# Add compile definitions based on the PERF option
if(PERF)
    message("Performance analysis mode is enabled. Run using Valgrind. Output disabled.")
    # Add compiler flags for performance analysis mode. This is NOT performance 
    # mode. For maximum performance, pass -DPROD=ON.
    add_compile_definitions(NO_OUTPUT)
endif()

# Add compile definitions based on the PROD option
if(PROD)
    message("Production mode is enabled.")
    # Add optimization flags 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -ffast-math -finline-functions -fomit-frame-pointer -DNDEBUG")
endif()
