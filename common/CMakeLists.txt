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

add_library(amcset-common STATIC
    ${PROJECT_SOURCE_DIR}/common/src/simulation.cpp
    ${PROJECT_SOURCE_DIR}/common/src/particle.cpp
    ${PROJECT_SOURCE_DIR}/common/src/electron.cpp
    ${PROJECT_SOURCE_DIR}/common/src/ion.cpp
    ${PROJECT_SOURCE_DIR}/common/src/quaternion.cpp
    ${PROJECT_SOURCE_DIR}/common/src/entropy_seeder.cpp
    ${PROJECT_SOURCE_DIR}/common/src/layer.cpp
    ${PROJECT_SOURCE_DIR}/common/src/volume.cpp
    ${PROJECT_SOURCE_DIR}/common/include/amcset_common.h
    ${PROJECT_SOURCE_DIR}/common/include/isotope_data.h
    ${PROJECT_SOURCE_DIR}/common/include/amcset_utilities.h
    ${PROJECT_SOURCE_DIR}/common/generated/isotope_database.h
    )

# state that common need PIC when the default is shared libraries
set_target_properties(amcset-common PROPERTIES
                        POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
                        )


target_link_libraries(amcset-common PUBLIC compiler_flags)

# define the symbol stating we are using the declspec(dllexport) when
# building on windows
target_compile_definitions(amcset-common PRIVATE "EXPORTING_AMCSET")

find_package(Python3 REQUIRED COMPONENTS Interpreter)

# Add the custom command to generate the isotope database
add_custom_command(
    OUTPUT ${PROJECT_SOURCE_DIR}/common/generated/isotope_database.h
    COMMAND
    ${Python3_EXECUTABLE}
    ${PROJECT_SOURCE_DIR}/common/scripts/generate_isotope_data.py
            ${PROJECT_SOURCE_DIR}/common/assets/isotope_data.txt
            ${PROJECT_SOURCE_DIR}/common/generated/isotope_database.h
    DEPENDS ${PROJECT_SOURCE_DIR}/common/assets/isotope_data.txt
            ${PROJECT_SOURCE_DIR}/common/scripts/generate_isotope_data.py
    COMMENT "Generating isotope database"
)


target_include_directories(amcset-common PUBLIC
    ${PROJECT_SOURCE_DIR}/common/include/
    ${PROJECT_SOURCE_DIR}/common/generated/
    )
