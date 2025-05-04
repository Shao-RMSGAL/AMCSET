cmake_minimum_required(VERSION 3.23)

option(ENABLE_ASAN "Enable AddressSanitizer" OFF)

if(ENABLE_ASAN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address")
endif()

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/version.txt" PROJECT_VERSION_CONTENT)
string(STRIP "${PROJECT_VERSION_CONTENT}" PROJECT_VERSION)
project(amcset VERSION ${PROJECT_VERSION})

# Use the conanfile.py version number
if(NOT DEFINED PROJECT_VERSION)
    message(FATAL_ERROR "Version components not defined! Are you running through Conan?")
endif()
# For clangd lsp
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Require C++20
add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_20)

# Set compiler flags
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
set(only_clang "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang>")

target_compile_options(compiler_flags INTERFACE
  "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused;>>"
   "$<${only_clang}:$<BUILD_INTERFACE:-fbracket-depth=3500>>"
  "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
)

# Build shared libraries set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}") set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Pass version number to compiler
configure_file(${PROJECT_SOURCE_DIR}/common/include/config.h.in ${PROJECT_SOURCE_DIR}/common/include/config.h)

# Packages
find_package(Boost 1.85.0 REQUIRED)
find_package(GTest 1.16.0 REQUIRED)
find_package(glog 0.7.1 REQUIRED)
find_package(gflags 2.2.2 REQUIRED)
find_package(Qt6 REQUIRED COMPONENTS Core Gui Widgets 3DCore 3DExtras 3DRender 3DInput Graphs DataVisualization)
if (Qt6_FOUND)
    qt_standard_project_setup()
endif ()

# Configure subdirectories
add_subdirectory(common)
add_subdirectory(gui)
add_subdirectory(test)

if(WIN32)
    # set(Boost_DIR "${CMAKE_BINARY_DIR}/build/generators")
    # Target Windows 11 (_WIN32_WINNT = 0x0A00 for Windows 10 and later)
    add_definitions(-D_WIN32_WINNT=0x0A00)

    # Optionally define _CRT_SECURE_NO_WARNINGS to suppress MSVC warnings
    # add_definitions(-D_CRT_SECURE_NO_WARNINGS)

    # If you're using C++17 or newer, set the standard explicitly
    set(CMAKE_CXX_STANDARD 23)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Optionally set compile options for MSVC
    if(MSVC)
        add_compile_options(/W4)  # High warning level
    endif()
endif()

# Add package includes and link common
target_link_libraries(amcset-common PUBLIC ${Boost_LIBRARIES} glog::glog)

# Allow link access for common library to gui
target_link_libraries(amcset-gui PRIVATE 
    amcset-common
    Qt6::DataVisualization
    # Qt6::Graphs
    Qt6::Charts
    Qt6::Core
    Qt6::Widgets
    Qt6::Gui
    Qt6::3DCore
    Qt6::3DExtras
    Qt6::3DRender
    Qt6::3DInput
    gflags::gflags
)
