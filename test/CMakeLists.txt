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


# common-test: unit, physics and regression tests of amcset-common.
add_executable(common-test
    ${PROJECT_SOURCE_DIR}/test/common/common_test.cpp
    ${PROJECT_SOURCE_DIR}/test/common/tallies_test.cpp
    ${PROJECT_SOURCE_DIR}/test/common/tallies_more_test.cpp
    ${PROJECT_SOURCE_DIR}/test/common/physics_test.cpp
    ${PROJECT_SOURCE_DIR}/test/common/trim_test.cpp
    ${PROJECT_SOURCE_DIR}/test/common/output_test.cpp
    ${PROJECT_SOURCE_DIR}/test/common/fast_math_test.cpp
    ${PROJECT_SOURCE_DIR}/test/common/diagnostics_test.cpp
    ${PROJECT_SOURCE_DIR}/test/common/crystal_test.cpp
    ${PROJECT_SOURCE_DIR}/test/common/compute_test.cpp
    )
target_link_libraries(common-test PRIVATE GTest::gtest_main amcset-common)

include(GoogleTest)
# Tests that read common/assets/*.csv expect the source tree as working
# directory. Optional SRIM-reference tests run when the environment variables
# AMCSET_SRIM_STOPPING_TABLE (and AMCSET_SRIM_DATA) are set; see TESTING.md.
gtest_discover_tests(common-test WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} DISCOVERY_TIMEOUT 60)

if(AMCSET_BUILD_CLI)
  # cli-test: run configuration and command-line parsing.
  add_executable(cli-test ${PROJECT_SOURCE_DIR}/test/cli/cli_test.cpp)
  target_link_libraries(cli-test PRIVATE GTest::gtest_main amcset-cli-lib)
  gtest_discover_tests(cli-test WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} DISCOVERY_TIMEOUT 60)

  # End-to-end smoke tests of the amcset executable.
  set(_cli_out ${CMAKE_CURRENT_BINARY_DIR}/cli-smoke)
  add_test(NAME cli.version COMMAND amcset version)
  set_tests_properties(cli.version PROPERTIES PASS_REGULAR_EXPRESSION "AMCSET [0-9]+\\.[0-9]+")
  add_test(NAME cli.template COMMAND amcset template)
  set_tests_properties(cli.template PROPERTIES PASS_REGULAR_EXPRESSION "\"layers\"")
  add_test(NAME cli.elements COMMAND amcset elements Fe)
  set_tests_properties(cli.elements PROPERTIES PASS_REGULAR_EXPRESSION "Fe")
  add_test(NAME cli.stopping COMMAND amcset stopping --ion H --target Fe --from "1 MeV" --to "1 MeV")
  add_test(NAME cli.template_save COMMAND amcset template srim -o ${_cli_out}/scenario.json)
  add_test(NAME cli.validate COMMAND amcset validate ${_cli_out}/scenario.json)
  set_tests_properties(cli.validate PROPERTIES DEPENDS cli.template_save)
  add_test(NAME cli.run_srim
           COMMAND amcset run ${_cli_out}/scenario.json -n 5 -j 2 -q -o ${_cli_out}/srim)
  set_tests_properties(cli.run_srim PROPERTIES DEPENDS cli.template_save)
  add_test(NAME cli.convert_trim_in
           COMMAND amcset convert ${_cli_out}/scenario.json ${_cli_out}/TRIM.IN)
  set_tests_properties(cli.convert_trim_in PROPERTIES DEPENDS cli.template_save)
  add_test(NAME cli.run_trim_in
           COMMAND amcset trim ${_cli_out}/TRIM.IN -n 5 -q -o ${_cli_out}/trim)
  set_tests_properties(cli.run_trim_in PROPERTIES DEPENDS cli.convert_trim_in)
  add_test(NAME cli.cell COMMAND amcset cell --structure bcc --sublattice Fe --a "2.8665 A" --beam "1 1 1")
  set_tests_properties(cli.cell PROPERTIES PASS_REGULAR_EXPRESSION "ion.angle = 54.73561 deg")
  add_test(NAME cli.cell_list COMMAND amcset cell --list)
  set_tests_properties(cli.cell_list PROPERTIES PASS_REGULAR_EXPRESSION "zincblende")
  # The grains of a polycrystal (written at configure time): mosaic grains in
  # registry form one group.
  file(WRITE ${_cli_out}/polycrystal.json [=[{"ion": {"element": "B", "energy": "30 keV"},
 "layers": [{"thickness": "600 A", "crystal": {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"],
   "polycrystal": {"grain_size": "150 A", "texture": "mosaic", "registry_below": "2 deg"}}}],
 "physics": {"transport": "crystal"}}
]=])
  add_test(NAME cli.grains COMMAND amcset grains ${_cli_out}/polycrystal.json)
  set_tests_properties(cli.grains PROPERTIES PASS_REGULAR_EXPRESSION "Registry +below 2 deg: 1 groups")
  add_test(NAME cli.devices COMMAND amcset devices)
  set_tests_properties(cli.devices PROPERTIES PASS_REGULAR_EXPRESSION "reference +0")
  # The GPU kernel on the CPU (reference backend) through the command line.
  add_test(NAME cli.run_reference
           COMMAND amcset run ${_cli_out}/scenario.json -n 32 -q --compute reference -o ${_cli_out}/reference)
  set_tests_properties(cli.run_reference PROPERTIES DEPENDS cli.template_save
                       PASS_REGULAR_EXPRESSION "Ran on +reference")
  add_test(NAME cli.bad_option COMMAND amcset run --no-such-option)
  set_tests_properties(cli.bad_option PROPERTIES WILL_FAIL TRUE)
  # Crash reports (diagnostics.h): a deliberate crash must produce a report
  # with the stack, the breadcrumb of the working thread and the rerun hint.
  # (The reports go next to the log in the build tree; no minidumps. On
  # Linux a fatal signal is raised again after the report, which CTest counts
  # as a crash, so the segv case runs on Windows only.)
  set(_crash_kinds thread-exception)
  if(WIN32)
    list(APPEND _crash_kinds segv)
  endif()
  foreach(_kind ${_crash_kinds})
    add_test(NAME cli.crash_report_${_kind}
             COMMAND amcset diag-crash ${_kind} --log ${_cli_out}/crash/diag.log)
    set_tests_properties(cli.crash_report_${_kind} PROPERTIES
      ENVIRONMENT "AMCSET_MINIDUMP=0"
      PASS_REGULAR_EXPRESSION "AMCSET crash report.*Stack of the crashing thread.*particle 41.*first_ion=32 -n 10")
  endforeach()
endif()
