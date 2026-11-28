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

# Configure through the Conan presets (see qbuild.sh / qbuild.ps1). Targets:
#   amcset-common  static library with the simulation core (common/)
#   amcset         command-line front end (cli/), no Qt needed
#   amcset-gui     Qt 6 desktop application (gui/)
#   common-test    GoogleTest suite (test/)
#
# Options:
#   AMCSET_BUILD_GUI    build the Qt GUI (needs Qt 6.7+ with QtGraphs and Qt3D)
#   AMCSET_BUILD_CLI    build the command-line program
#   AMCSET_BUILD_TESTS  build the GoogleTest suite

# Keep in sync with conanfile.py's version.
project(amcset VERSION 1.0.1
        DESCRIPTION "Monte Carlo transport of ions and electrons and radiation damage in layered targets"
        HOMEPAGE_URL "https://github.com/Shao-RMSGAL/AMCSET"
        LANGUAGES C CXX)

option(AMCSET_BUILD_GUI "Build the Qt GUI" ON)
option(AMCSET_BUILD_CLI "Build the command-line program" ON)
option(AMCSET_BUILD_TESTS "Build the test suite" ON)

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
  # No fused multiply-adds where the source has none (MSVC's default): GCC
  # fuses wherever -march allows, differently in code that computes the same
  # expressions, which broke the bit-for-bit comparison of the CPU transport
  # with the GPU kernel compiled as C++ and exact equalities of event times.
  "$<${gcc_like_cxx}:-ffp-contract=off>"
)

# Performance: tune Release flags for maximum runtime performance.
# NOTE: the default -march=native makes Release binaries specific to the build machine's
# CPU, and these options are applied twice (via CMAKE_CXX_FLAGS_RELEASE and
# via compiler_flags); harmless but redundant.
# Set AMCSET_MARCH=x86-64-v2 for binaries that run on any recent x86-64 CPU
# (the Linux AppImage, tools/linux/make_appimage.sh).
set(AMCSET_MARCH "native" CACHE STRING "-march value of GCC/Clang Release builds")
# Portable Linux binaries (the CLI tarball of a release): the C++ runtime
# linked in, so that an older system's libstdc++ is enough.
option(AMCSET_STATIC_LIBSTDCXX "Link libstdc++ and libgcc statically (GCC)" OFF)
if(AMCSET_STATIC_LIBSTDCXX AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_link_options(-static-libstdc++ -static-libgcc)
endif()
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

if(MSVC)
    # /Zi + /DEBUG write a PDB next to each executable (no effect on the
    # generated code; /OPT:REF,ICF keep the Release linker defaults). The
    # crash handler and profilers need it for readable stack traces.
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /GL /DNDEBUG /Zi")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG /DEBUG /OPT:REF /OPT:ICF")
    target_compile_options(compiler_flags INTERFACE
        $<$<CONFIG:Release>:/O2;/GL>
    )
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=${AMCSET_MARCH} -flto -funroll-loops -pipe -DNDEBUG")
    # GCC: -flto=auto runs the link-time code generation on every core (plain
    # -flto runs it serially and warns about it).
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(AMCSET_LTO_LINK_FLAG "-flto=auto")
    else()
        set(AMCSET_LTO_LINK_FLAG "-flto")
    endif()
    # (-ffp-contract=off also at the link: the link-time code generation.)
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} ${AMCSET_LTO_LINK_FLAG} -ffp-contract=off")
    target_compile_options(compiler_flags INTERFACE
        $<$<CONFIG:Release>:-O3;-march=${AMCSET_MARCH};-flto;-funroll-loops;-pipe>
    )
endif()
# ABSL_MIN_LOG_LEVEL=2 compiles out LOG(INFO) and LOG(WARNING) in Release, so
# Release log files only contain errors.
target_compile_definitions(compiler_flags INTERFACE
    $<$<CONFIG:Release>:NDEBUG>
    $<$<CONFIG:Release>:ABSL_MIN_LOG_LEVEL=2>
    $<$<PLATFORM_ID:Windows>:_WIN32_WINNT=0x0A00>
    $<$<PLATFORM_ID:Windows>:_CRT_SECURE_NO_WARNINGS>
)
target_link_options(compiler_flags INTERFACE
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:MSVC>>:/LTCG>
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU>>:-flto=auto>
    $<$<AND:$<CONFIG:Release>,$<NOT:$<CXX_COMPILER_ID:MSVC,GNU>>>:-flto>
)

# Put every artifact directly in the build directory (e.g. build/Release).
# BUILD_SHARED_LIBS only affects PIC for amcset-common, which is always a
# STATIC library.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Pass version number to compiler (writes into the source tree; config.h is
# git-ignored).
configure_file(${PROJECT_SOURCE_DIR}/common/include/config.h.in ${PROJECT_SOURCE_DIR}/common/include/config.h)

# Packages
find_package(Boost 1.83 REQUIRED)
find_package(absl REQUIRED CONFIG)
find_package(Threads REQUIRED)
if(AMCSET_BUILD_TESTS)
    find_package(GTest REQUIRED)
endif()
if(AMCSET_BUILD_GUI)
    # Qt Widgets for the application, QRhi (Qt::GuiPrivate) for the GPU views.
    # Qt Shader Tools is optional: without it the pre-compiled shaders in
    # gui/shaders/baked are used.
    find_package(Qt6 6.7 COMPONENTS Core Gui Widgets REQUIRED)
    # The path views use QRhi, whose headers belong to GuiPrivate; the
    # executable is tied to the Qt version it was built with, which windeployqt
    # ships alongside it, so the private-module warning is silenced.
    set(QT_NO_PRIVATE_MODULE_WARNING ON)
    find_package(Qt6 6.7 QUIET COMPONENTS GuiPrivate ShaderTools)
    qt_standard_project_setup()
endif()

# Install layout (bin/, share/amcset/, share/doc/amcset/, ...) and the
# executables' Windows resources.
include(GNUInstallDirs)
include(${PROJECT_SOURCE_DIR}/packaging/windows/resources.cmake)

add_subdirectory(common)
if(AMCSET_BUILD_CLI)
    add_subdirectory(cli)
endif()
if(AMCSET_BUILD_GUI)
    add_subdirectory(gui)
endif()

# Profile-guided optimisation of the command-line program (MSVC; see
# docs/DEVELOPMENT.md "Compilers and flags"): configure with
# -DAMCSET_PGO=generate, build, run the training cases, then reconfigure with
# -DAMCSET_PGO=use and build again. The profile goes to AMCSET_PGO_DIR.
set(AMCSET_PGO "" CACHE STRING "Profile-guided optimisation of amcset: empty, generate or use")
set(AMCSET_PGO_DIR "${PROJECT_BINARY_DIR}/pgo" CACHE PATH "Directory of the PGO profile")
if(AMCSET_PGO AND AMCSET_BUILD_CLI)
    if(NOT MSVC)
        message(FATAL_ERROR "AMCSET_PGO is implemented for MSVC only")
    endif()
    file(MAKE_DIRECTORY ${AMCSET_PGO_DIR})
    if(AMCSET_PGO STREQUAL "generate")
        target_link_options(amcset PRIVATE /GENPROFILE:PGD=${AMCSET_PGO_DIR}/amcset.pgd)
    elseif(AMCSET_PGO STREQUAL "use")
        target_link_options(amcset PRIVATE /USEPROFILE:PGD=${AMCSET_PGO_DIR}/amcset.pgd)
    else()
        message(FATAL_ERROR "AMCSET_PGO must be empty, generate or use")
    endif()
endif()
if(AMCSET_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()

# Installation (desktop integration, documents, runtime libraries) and the
# packages built from it with CPack: docs/RELEASING.md.
add_subdirectory(packaging)
