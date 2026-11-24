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

# Stand-alone build of the determinism harness (no Conan, no Qt). See
# docs/TESTING.md. Requires Boost headers and abseil (with absl::log).
#
#   cmake -S tools/determinism_harness -B build/harness -DCMAKE_BUILD_TYPE=Release
#   cmake --build build/harness
#   ./build/harness/amcset-harness > before.txt   # run from the repo root

cmake_minimum_required(VERSION 3.23)
project(amcset_determinism_harness CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Root of the AMCSET checkout whose common/ sources are compiled in.
set(AMCSET_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE PATH
    "AMCSET source tree to build the harness against")

find_package(Boost REQUIRED)
find_package(absl REQUIRED CONFIG)
find_package(Threads REQUIRED)

# common/src/physics and common/src/trim are part of the core since v0.1.0.
file(GLOB_RECURSE AMCSET_COMMON_SOURCES ${AMCSET_SOURCE_DIR}/common/src/*.cpp)
set(amcset_VERSION_MAJOR 0)
set(amcset_VERSION_MINOR 1)
set(amcset_VERSION_PATCH 0)
configure_file(${AMCSET_SOURCE_DIR}/common/include/config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/generated/config.h)
add_executable(amcset-harness harness.cpp ${AMCSET_COMMON_SOURCES})
target_include_directories(amcset-harness PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/generated
    ${AMCSET_SOURCE_DIR}/common/include
    ${AMCSET_SOURCE_DIR}/common/generated)
target_link_libraries(amcset-harness PRIVATE
    Boost::boost absl::log absl::log_initialize Threads::Threads)
