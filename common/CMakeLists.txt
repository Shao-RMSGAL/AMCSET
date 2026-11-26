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


# amcset-common: the simulation core shared by the GUI, the CLI and the tests.
add_library(amcset-common STATIC
    # Native AMCSET engine
    ${PROJECT_SOURCE_DIR}/common/src/simulation.cpp
    ${PROJECT_SOURCE_DIR}/common/src/particle.cpp
    ${PROJECT_SOURCE_DIR}/common/src/electron.cpp
    ${PROJECT_SOURCE_DIR}/common/src/ion.cpp
    ${PROJECT_SOURCE_DIR}/common/src/quaternion.cpp
    ${PROJECT_SOURCE_DIR}/common/src/entropy_seeder.cpp
    ${PROJECT_SOURCE_DIR}/common/src/layer.cpp
    ${PROJECT_SOURCE_DIR}/common/src/volume.cpp
    # Shared physics and data
    ${PROJECT_SOURCE_DIR}/common/src/element_data.cpp
    ${PROJECT_SOURCE_DIR}/common/src/run_config.cpp
    ${PROJECT_SOURCE_DIR}/common/src/data_files.cpp
    ${PROJECT_SOURCE_DIR}/common/src/diagnostics.cpp
    ${PROJECT_SOURCE_DIR}/common/src/physics/stopping.cpp
    ${PROJECT_SOURCE_DIR}/common/src/physics/scattering.cpp
    # SRIM mode (TRIM-compatible engine)
    ${PROJECT_SOURCE_DIR}/common/src/trim/trim_engine.cpp
    ${PROJECT_SOURCE_DIR}/common/src/trim/target_depth.cpp
    ${PROJECT_SOURCE_DIR}/common/src/trim/trim_output.cpp
    ${PROJECT_SOURCE_DIR}/common/src/trim/trim_export.cpp
    ${PROJECT_SOURCE_DIR}/common/src/trim/trim_in.cpp
    ${PROJECT_SOURCE_DIR}/common/src/trim/trim_dat.cpp
    ${PROJECT_SOURCE_DIR}/common/src/trim/single_precision.cpp
    ${PROJECT_SOURCE_DIR}/common/src/trim/trim_tables.cpp
    ${PROJECT_SOURCE_DIR}/common/src/trim/layer_crystal.cpp
    # Crystals (unit cells, orientation, cell files; not TRIM)
    ${PROJECT_SOURCE_DIR}/common/src/crystal/unit_cell.cpp
    ${PROJECT_SOURCE_DIR}/common/src/crystal/orientation.cpp
    ${PROJECT_SOURCE_DIR}/common/src/crystal/crystal.cpp
    ${PROJECT_SOURCE_DIR}/common/src/crystal/cell_io.cpp
    ${PROJECT_SOURCE_DIR}/common/src/crystal/lattice_search.cpp
    ${PROJECT_SOURCE_DIR}/common/src/crystal/materials.cpp
    ${PROJECT_SOURCE_DIR}/common/src/crystal/grain_map.cpp
    ${PROJECT_SOURCE_DIR}/common/src/crystal/grain_io.cpp
    # Output tables and file formats
    ${PROJECT_SOURCE_DIR}/common/src/output/units.cpp
    ${PROJECT_SOURCE_DIR}/common/src/output/table.cpp
    ${PROJECT_SOURCE_DIR}/common/src/output/stream.cpp
    ${PROJECT_SOURCE_DIR}/common/src/output/writers.cpp
    ${PROJECT_SOURCE_DIR}/common/src/output/parquet.cpp
    # Headers (for IDEs)
    ${PROJECT_SOURCE_DIR}/common/include/amcset_common.h
    ${PROJECT_SOURCE_DIR}/common/include/amcset_units.h
    ${PROJECT_SOURCE_DIR}/common/include/amcset_utilities.h
    ${PROJECT_SOURCE_DIR}/common/include/element_data.h
    ${PROJECT_SOURCE_DIR}/common/include/run_config.h
    ${PROJECT_SOURCE_DIR}/common/include/data_files.h
    ${PROJECT_SOURCE_DIR}/common/include/diagnostics.h
    ${PROJECT_SOURCE_DIR}/common/include/parallel_tasks.h
    ${PROJECT_SOURCE_DIR}/common/include/isotope_data.h
    ${PROJECT_SOURCE_DIR}/common/include/simulation_engine.h
    ${PROJECT_SOURCE_DIR}/common/include/tallies.h
    ${PROJECT_SOURCE_DIR}/common/include/physics/stopping.h
    ${PROJECT_SOURCE_DIR}/common/include/physics/scattering.h
    ${PROJECT_SOURCE_DIR}/common/include/physics/fast_math.h
    ${PROJECT_SOURCE_DIR}/common/include/trim/trim.h
    ${PROJECT_SOURCE_DIR}/common/include/trim/trim_in.h
    ${PROJECT_SOURCE_DIR}/common/include/trim/trim_dat.h
    ${PROJECT_SOURCE_DIR}/common/include/trim/single_precision.h
    ${PROJECT_SOURCE_DIR}/common/include/trim/trim_tables.h
    ${PROJECT_SOURCE_DIR}/common/include/output/units.h
    ${PROJECT_SOURCE_DIR}/common/include/output/table.h
    ${PROJECT_SOURCE_DIR}/common/include/output/text.h
    ${PROJECT_SOURCE_DIR}/common/include/output/stream.h
    ${PROJECT_SOURCE_DIR}/common/include/output/writers.h
    ${PROJECT_SOURCE_DIR}/common/include/trim/random.h
    ${PROJECT_SOURCE_DIR}/common/generated/isotope_database.h
    ${PROJECT_SOURCE_DIR}/common/generated/element_database.h
    )

set_target_properties(amcset-common PROPERTIES
                        POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
                        )

target_link_libraries(amcset-common PUBLIC
    compiler_flags
    Boost::boost
    absl::log
    absl::log_initialize
    absl::log_sink
    absl::log_sink_registry
    absl::log_entry
    absl::stacktrace
    absl::symbolize
    Threads::Threads
    )
# Crash reports: stack walking, symbols and minidumps (diagnostics.cpp).
if(WIN32)
    target_link_libraries(amcset-common PUBLIC dbghelp)
endif()

# Data files of the source tree (data/), found at run time by data_files.cpp
# when the program runs from the build folder.
target_compile_definitions(amcset-common PRIVATE AMCSET_SOURCE_DATA_DIR="${PROJECT_SOURCE_DIR}/data"
                                                  AMCSET_SOURCE_ASSETS_DIR="${PROJECT_SOURCE_DIR}/common/assets")
install(DIRECTORY ${PROJECT_SOURCE_DIR}/data/srim2013 DESTINATION ${CMAKE_INSTALL_DATADIR}/amcset)
# The native engine's electron tables (data_files.h: find_electron_table).
install(FILES ${PROJECT_SOURCE_DIR}/common/assets/electron_screening_parameters.csv
              ${PROJECT_SOURCE_DIR}/common/assets/mott_scattering_parameters.csv
        DESTINATION ${CMAKE_INSTALL_DATADIR}/amcset/electrons)

target_include_directories(amcset-common PUBLIC
    ${PROJECT_SOURCE_DIR}/common/include/
    ${PROJECT_SOURCE_DIR}/common/generated/
    )

# The generated databases are committed; they are regenerated when their
# source table or generator changes.
find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_custom_command(
        OUTPUT ${PROJECT_SOURCE_DIR}/common/generated/isotope_database.h
        COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/common/scripts/generate_isotope_data.py
                ${PROJECT_SOURCE_DIR}/common/assets/isotope_data.txt
                ${PROJECT_SOURCE_DIR}/common/generated/isotope_database.h
        DEPENDS ${PROJECT_SOURCE_DIR}/common/assets/isotope_data.txt
                ${PROJECT_SOURCE_DIR}/common/scripts/generate_isotope_data.py
        COMMENT "Generating isotope database"
    )
    add_custom_command(
        OUTPUT ${PROJECT_SOURCE_DIR}/common/generated/element_database.h
        COMMAND ${Python3_EXECUTABLE}
                ${PROJECT_SOURCE_DIR}/common/scripts/generate_element_data.py
                ${PROJECT_SOURCE_DIR}/common/assets/element_data.csv
                ${PROJECT_SOURCE_DIR}/common/generated/element_database.h
        DEPENDS ${PROJECT_SOURCE_DIR}/common/assets/element_data.csv
                ${PROJECT_SOURCE_DIR}/common/scripts/generate_element_data.py
        COMMENT "Generating element database"
    )
endif()
