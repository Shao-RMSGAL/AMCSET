cmake_policy(SET CMP0077 NEW)

cmake_minimum_required(VERSION 3.22.1)
project(eSRIM)

set(CMAKE_CXX_STANDARD 20)

# Define compiler options
set(COMPILER "gcc" CACHE STRING "Choose the compiler: 'gcc' or 'clang' (default is 'gcc')")

# Check for clang compiler
if(COMPILER STREQUAL "clang")
    find_program(CLANG_FOUND clang)
    find_program(CLANGXX_FOUND clang++)
    if(NOT CLANG_FOUND OR NOT CLANGXX_FOUND)
        message(FATAL_ERROR "clang compiler not found.")
        return()
    endif()
# Check for gcc compiler
elseif(COMPILER STREQUAL "gcc")
    find_program(GCC_FOUND gcc)
    find_program(GXX_FOUND g++)
    if(NOT GCC_FOUND OR NOT GXX_FOUND)
        message(FATAL_ERROR "gcc compiler not found.")
        return()
    endif()
# Check intel
elseif(COMPILER STREQUAL "intel")
#     find_program(INTEL_FOUND icpx)
#     if(NOT INTEL_FOUND)
#         message(FATAL_ERROR "intel compiler not found.")
#         return()
#     endif()
elseif(NOT COMPILER STREQUAL "gcc")
    message(FATAL_ERROR "Unrecognized compiler name: ${COMPILER}")
    return()
endif()

# Set compiler
if(COMPILER STREQUAL "clang")
    set(CMAKE_C_COMPILER "clang")
    set(CMAKE_CXX_COMPILER "clang++")
elseif(COMPILER STREQUAL "gcc")
    set(CMAKE_C_COMPILER "gcc")
    set(CMAKE_CXX_COMPILER "g++")
elseif(COMPILER STREQUAL "intel")
    set(CMAKE_CXX_COMPILER "icx-cl")
endif()

# Display selected compiler
message(STATUS "Using ${COMPILER} compiler")

# Create a 'bin' directory for the build
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# Set the path to the source directory
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)

# Include directory for header files
include_directories(include)

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
option(DEBUG_MODE "Enable debug mode. Passing -g and -O0 flags." OFF)
option(PERF "Enable performance analysis mode. Passing -gp -O0. Use either Valgrind or gprof for performance analysis." OFF)
option(PROD "Enable production mode. Passing  -Ofast -march=native -ffast-math -finline-functions -fomit-frame-pointer -DNDEBUG" ON)

# Add compile definitions based on options
if(DEBUG_MODE)
    message(STATUS "Debug mode is enabled.")
    
    # Add DEBUG_MODE definition to enable debug messages
    add_compile_definitions(DEBUG_MODE)
    
    # Add compiler flags for debugging with gdb. 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")

elseif(PERF) # Add compile definitions based on the PERF option
    message(STATUS "Performance analysis mode is enabled.")

    # Add NO_OUTPUT definition to disable output file creation. This prevents file output from modifying performance operations.
    add_compile_definitions(NO_OUTPUT)

    # Add compiler flags for performance analysis mode. This is NOT production mode. Use -DPROD=ON to enable production mode. 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -O0")

elseif(PROD) # Add compile definitions based on the PROD option
    message(STATUS "Production mode is enabled.")

    # Add optimization flags 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -ffast-math -finline-functions -fomit-frame-pointer -DNDEBUG")
else() # Add compile definitions based on the PROD option
    message(STATUS "No options set. Using -O0 compilation.")

    # Add optimization flags 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif()