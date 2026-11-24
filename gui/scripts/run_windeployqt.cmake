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

# Post-build step of amcset-gui on Windows: copy the Qt runtime next to the
# executable with windeployqt. A failure (for example windeployqt blocked by
# an application-control policy) is reported as a warning and does not fail
# the build; the GUI then needs the Qt DLLs on PATH (e.g. via conanrun.bat).
#
# Usage: cmake -DWINDEPLOYQT=<path> -DTARGET=<exe> -P run_windeployqt.cmake
if(NOT WINDEPLOYQT OR NOT EXISTS "${WINDEPLOYQT}")
  message(WARNING "windeployqt not found; Qt DLLs were not copied next to ${TARGET}")
  return()
endif()
execute_process(COMMAND "${WINDEPLOYQT}" --no-translations "${TARGET}"
                RESULT_VARIABLE result OUTPUT_QUIET ERROR_VARIABLE error)
if(NOT result EQUAL 0)
  message(WARNING "windeployqt failed (${result}): ${error}\n"
                  "Qt DLLs were not copied next to ${TARGET}; run the GUI with the "
                  "Qt runtime on PATH (build\\Release\\generators\\conanrun.bat).")
endif()
