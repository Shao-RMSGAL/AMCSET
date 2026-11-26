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


# Command-line program. The parsing/option logic lives in a small library so
# that the test suite can exercise it.
add_library(amcset-cli-lib STATIC
    ${PROJECT_SOURCE_DIR}/cli/src/cli.cpp
    ${PROJECT_SOURCE_DIR}/cli/include/cli.h
    ${PROJECT_SOURCE_DIR}/cli/src/interrupts.cpp
    ${PROJECT_SOURCE_DIR}/cli/include/interrupts.h
    )
target_include_directories(amcset-cli-lib PUBLIC ${PROJECT_SOURCE_DIR}/cli/include)
target_link_libraries(amcset-cli-lib PUBLIC amcset-common)

add_executable(amcset ${PROJECT_SOURCE_DIR}/cli/src/main.cpp)
target_link_libraries(amcset PRIVATE amcset-cli-lib)
# Version information, icon and manifest of the executable (Windows).
amcset_windows_resources(amcset "AMCSET command-line program")
install(TARGETS amcset RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
