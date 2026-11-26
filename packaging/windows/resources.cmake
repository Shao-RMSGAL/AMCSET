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

# amcset_windows_resources(<target> <description>): on Windows, gives the
# executable <target> its version information and icon (packaging/windows/
# version.rc.in, resources/icons/amcset.ico) and the application manifest
# (packaging/windows/amcset.manifest, merged by the linker). Nothing elsewhere.
function(amcset_windows_resources target description)
  if(NOT WIN32)
    return()
  endif()
  get_target_property(RC_NAME ${target} OUTPUT_NAME)
  if(NOT RC_NAME)
    set(RC_NAME ${target})
  endif()
  set(RC_DESCRIPTION "${description}")
  file(TO_CMAKE_PATH "${PROJECT_SOURCE_DIR}/resources/icons/amcset.ico" icon)
  set(RC_ICON_LINE "1 ICON \"${icon}\"")
  configure_file("${PROJECT_SOURCE_DIR}/packaging/windows/version.rc.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/${target}_version.rc" @ONLY)
  target_sources(${target} PRIVATE
    "${CMAKE_CURRENT_BINARY_DIR}/${target}_version.rc"
    "${PROJECT_SOURCE_DIR}/packaging/windows/amcset.manifest")
endfunction()
