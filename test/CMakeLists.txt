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

enable_testing()

# Common test
add_executable( common-test
    ${PROJECT_SOURCE_DIR}/test/common/common_test.cpp
    )
# Add tallies integration test to the common test binary
target_sources(common-test PRIVATE ${PROJECT_SOURCE_DIR}/test/common/tallies_test.cpp)
target_sources(common-test PRIVATE ${PROJECT_SOURCE_DIR}/test/common/tallies_more_test.cpp)

target_link_libraries(common-test PRIVATE
    GTest::gtest_main
    amcset-common
    )

target_include_directories( common-test PRIVATE
    ${PROJECT_SOURCE_DIR}/common/include
    )

include(GoogleTest)
gtest_discover_tests(common-test)
