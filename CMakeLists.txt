# Check if cmake is run from the project root directory
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(FATAL_ERROR "CMake should not be run from the project root directory.\n" 
        "Please run 'cmake ..' from the project build directory.\n"
        "If you previously attempted to run cmake from the project root directory,"
        "delete the auto-generated CMakeCache.txt file and CMakeFiles directory and try again.")
endif()

cmake_minimum_required(VERSION 3.23.5)

project(eSRIM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_C_COMPILER "gcc" CACHE STRING "Choose the C compiler: 'gcc', 'clang', or 'icx' (default is 'gcc')")
set(CMAKE_CXX_COMPILER "g++" CACHE STRING "Choose the C++ compiler: 'g++', 'clang++', or 'icpx' (default is 'g++')")

cmake_policy(SET CMP0135 NEW)

# gtest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/f8d7d77c06936315286eb55f8de22cd23c188571.zip
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

message(STATUS "Using ${CMAKE_C_COMPILER} C compiler")
message(STATUS "Using ${CMAKE_CXX_COMPILER} C++ compiler")

# Create a 'bin' directory for the build
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(SOURCE_DIR ${CMAKE_SOURCE_DIR}/src)
set(INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include)

set(SOURCE_FILES 
    ${SOURCE_DIR}/eSRIM.cpp
    ${SOURCE_DIR}/eSRIM_classes.cpp
    ${SOURCE_DIR}/utilities.cpp
)

set(MAIN_SOURCE_FILE 
    ${SOURCE_DIR}/main.cpp
)

add_executable(eSRIM ${MAIN_SOURCE_FILE} ${SOURCE_FILES})
include_directories(${INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -pthread")

if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    find_package(IntelSYCL REQUIRED)
    # Add your conditional code here that should only execute when using icpx
endif()

option(DEBUG_MODE "Enable debug mode. Passing -g and -O0 flags." OFF)
option(PERF "Enable performance analysis mode. Passing -gp -O0. Use either Valgrind or gprof for performance analysis." OFF)
option(PROD "Enable production mode. Passing  -Ofast -march=native -ffast-math -finline-functions -fomit-frame-pointer -DNDEBUG" ON)

if(DEBUG_MODE)
    message(STATUS "Debug mode is enabled.")
    add_compile_definitions(DEBUG_MODE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
elseif(PERF)
    message(STATUS "Performance analysis mode is enabled.")
    add_compile_definitions(NO_OUTPUT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
elseif(PROD)
    message(STATUS "Production mode is enabled.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -ffast-math -flto -finline-functions -fomit-frame-pointer -DNDEBUG")
    if(CMAKE_CXX_COMPILER STREQUAL "IntelLLVM")
        find_package(IntelSYCL REQUIRED)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost -fsycl -ipo -fp-model fast=2 -parallel -qopenmp")
    endif()
else()
    message(STATUS "No options set. Using -O0 compilation.")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif()

if(NOT (CMAKE_C_COMPILER STREQUAL "clang" AND CMAKE_CXX_COMPILER STREQUAL "clang++"))
    enable_testing()
    set(TEST_SOURCE_DIR ${CMAKE_SOURCE_DIR}/tests)

    # Add testing files
    set(TEST_SOURCE_FILES
        ${TEST_SOURCE_DIR}/mainFunctionTest.cpp
    )

    add_executable(eSRIMTests ${TEST_SOURCE_FILES} ${SOURCE_FILES})

    target_link_libraries(eSRIMTests gtest)
    add_test(NAME MainTest COMMAND eSRIMTests)

    # Include source and header directories
    target_include_directories(eSRIMTests PUBLIC ${INCLUDE_DIR})
    target_sources(eSRIMTests PUBLIC ${SOURCE_FILES})
endif()