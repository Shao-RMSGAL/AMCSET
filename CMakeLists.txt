# Copyright 2025, Texas A&M University
#
# This file is part of AMCSET.
#
# AMCSET is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# AMCSET. If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.23)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
# set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address")

project(amcset VERSION 0.0.1)

option(ENABLE_ASAN "Enable AddressSanitizer" OFF)
if(ENABLE_ASAN)
    if(MSVC)
        add_compile_options(/fsanitize=address)
        add_link_options(/fsanitize=address)
    else()
        add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
        add_link_options(-fsanitize=address)
    endif()
endif()

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

# Performance: tune Release flags for maximum runtime performance
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

if(MSVC)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /GL /DNDEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
    target_compile_options(compiler_flags INTERFACE
        $<$<CONFIG:Release>:/O2;/GL>
    )
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -flto -funroll-loops -pipe -DNDEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")
    target_compile_options(compiler_flags INTERFACE
        $<$<CONFIG:Release>:-O3;-march=native;-flto;-funroll-loops;-pipe>
    )
endif()
target_compile_definitions(compiler_flags INTERFACE
    $<$<CONFIG:Release>:NDEBUG>
    $<$<CONFIG:Release>:ABSL_MIN_LOG_LEVEL=2>
    $<$<PLATFORM_ID:Windows>:_WIN32_WINNT=0x0A00>
    $<$<PLATFORM_ID:Windows>:_CRT_SECURE_NO_WARNINGS>
)
target_link_options(compiler_flags INTERFACE
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:/LTCG>
    $<$<AND:$<CONFIG:Release>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:-flto>
)

# Build shared libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Pass version number to compiler
configure_file(${PROJECT_SOURCE_DIR}/common/include/config.h.in ${PROJECT_SOURCE_DIR}/common/include/config.h)

# Packages
find_package(Boost 1.86.0 REQUIRED)
find_package(GTest 1.15.0 REQUIRED)
find_package(absl REQUIRED CONFIG)
find_package(gflags 2.2.2 REQUIRED)
# find_package(OpenSSL REQUIRED CONFIG COMPONENTS filesystem)
find_package(Qt6 6.11.1 COMPONENTS Core Gui Widgets Graphs 3DCore 3DExtras 3DRender 3DInput REQUIRED)
if (Qt6_FOUND)
    qt_standard_project_setup()
endif ()

# Configure subdirectories
add_subdirectory(common)
add_subdirectory(gui)
add_subdirectory(test)

if(WIN32)
    # Keep the project standard aligned with the Conan profile.
    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)

    # Optionally set compile options for MSVC
    if(MSVC)
        add_compile_options(/W4)  # High warning level
    endif()

endif()

# Add package includes and link common
# target_include_directories(amcset-common SYSTEM PUBLIC
#     ${Boost_INCLUDE_DIRS}
#     )
target_link_libraries(amcset-common PUBLIC
    ${Boost_LIBRARIES}
    absl::log
    absl::log_initialize
)


# Allow link access for common library to gui
target_include_directories(amcset-gui PRIVATE)
target_link_libraries(amcset-gui PRIVATE 
    amcset-common
    Qt6::Core
    Qt6::Widgets
    Qt6::Gui
    Qt6::Graphs
    Qt6::3DCore
    Qt6::3DExtras
    Qt6::3DRender
    Qt6::3DInput
    gflags::gflags
)