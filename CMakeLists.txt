# Intel compiler lines
if (CMAKE_HOST_WIN32)
# need CMake 3.25.0+ for IntelLLVM support of target link properties on Windows
cmake_minimum_required(VERSION 3.25)
else()
# CMake 3.23.5 is the minimum recommended for IntelLLVM on Linux
cmake_minimum_required(VERSION 3.23.5)
endif()

project(eSRIM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

# Define compiler options
set(CMAKE_C_COMPILER "gcc" CACHE STRING "Choose the C compiler: 'gcc', 'clang', or 'icx' (default is 'gcc')")
set(CMAKE_CXX_COMPILER "g++" CACHE STRING "Choose the C++ compiler: 'g++', 'clang++', or 'icpx' (default is 'g++')")

cmake_policy(SET CMP0135 NEW)

# gtest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/f8d7d77c06936315286eb55f8de22cd23c188571.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Display selected compiler
message(STATUS "Using ${CMAKE_C_COMPILER} C compiler")
message(STATUS "Using ${CMAKE_CXX_COMPILER} C++ compiler")

# Create a 'bin' directory for the build
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/bin)
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
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build/bin)

# Intel-specific checks
if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
    find_package(IntelSYCL REQUIRED)
    # Add your conditional code here that should only execute when using icpx
endif()

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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")

elseif(PROD) # Add compile definitions based on the PROD option
    message(STATUS "Production mode is enabled.")

    # Add optimization flags 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -march=native -ffast-math -flto -finline-functions -fomit-frame-pointer -DNDEBUG")

    if(CMAKE_CXX_COMPILER STREQUAL "IntelLLVM")
        find_package(IntelSYCL REQUIRED)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost -fsycl")
    endif()

else() # Add compile definitions based on the PROD option
    message(STATUS "No options set. Using -O0 compilation.")

    # Add optimization flags 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
endif()

enable_testing()

# Define test executable
add_executable(mytests 
  tests/mainTest.cpp
  # Add more test files if needed
)

# Link Google Test and your project's libraries to the test executable
target_link_libraries(mytests gtest_main)

# Add tests to CMake
add_test(NAME MainTest COMMAND mytests)
# Add more tests as needed

add_custom_target(clean-cmake
    COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_BINARY_DIR}/CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/cmake_install.cmake
    COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/Makefile
    COMMENT "Cleaning CMake-generated files"
)
