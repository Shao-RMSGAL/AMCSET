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

# The AMCSET GUI: Qt Widgets, with the particle-path views drawn through QRhi
# (Vulkan, Direct3D, Metal or OpenGL).

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(GUI_SOURCES
    src/main.cpp
    src/app/main_window.cpp src/app/main_window.h
    src/app/guard.cpp src/app/guard.h
    src/app/unix_signals.cpp src/app/unix_signals.h
    src/model/scenario.cpp src/model/scenario.h
    src/model/target_models.cpp src/model/target_models.h
    src/model/quantity_text.cpp src/model/quantity_text.h
    src/model/output_location.cpp src/model/output_location.h
    src/run/run_controller.cpp src/run/run_controller.h
    src/run/track_store.cpp src/run/track_store.h
    src/run/datasets.cpp src/run/datasets.h
    src/setup/form_binder.h
    src/setup/ion_panel.cpp src/setup/ion_panel.h
    src/setup/target_panel.cpp src/setup/target_panel.h
    src/setup/calculation_panel.cpp src/setup/calculation_panel.h
    src/setup/crystal_panel.cpp src/setup/crystal_panel.h
    src/setup/grains_box.cpp src/setup/grains_box.h
    src/setup/output_panel.cpp src/setup/output_panel.h
    src/setup/output_spec_dialog.cpp src/setup/output_spec_dialog.h
    src/setup/periodic_table.cpp src/setup/periodic_table.h
    src/views/axis_ticks.cpp src/views/axis_ticks.h
    src/views/cell_view.cpp src/views/cell_view.h
    src/views/plot_widget.cpp src/views/plot_widget.h
    src/views/plots_panel.cpp src/views/plots_panel.h
    src/views/tables_panel.cpp src/views/tables_panel.h
    src/views/stats_panel.cpp src/views/stats_panel.h
    src/views/track_view.cpp src/views/track_view.h
    src/views/track_panel.cpp src/views/track_panel.h
    src/views/sample_scene.cpp src/views/sample_scene.h
)

add_executable(amcset-gui WIN32 ${GUI_SOURCES} resources/icons.qrc)
# Version information, icon and manifest of the executable (Windows).
amcset_windows_resources(amcset-gui "AMCSET")
install(TARGETS amcset-gui RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
target_include_directories(amcset-gui PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(amcset-gui PRIVATE compiler_flags amcset-common Qt6::Core Qt6::Gui Qt6::Widgets)

# QRhi lives in Qt's semi-public GUI headers (<rhi/qrhi.h>, Qt::GuiPrivate).
# Package managers that do not export Qt6::GuiPrivate still ship the headers
# under include/QtGui/<version>/QtGui.
if(TARGET Qt6::GuiPrivate)
  target_link_libraries(amcset-gui PRIVATE Qt6::GuiPrivate)
else()
  get_target_property(_qt_gui_includes Qt6::Gui INTERFACE_INCLUDE_DIRECTORIES)
  # Conan's CMakeDeps wraps the paths in generator expressions.
  string(REGEX REPLACE "\\$<[^:>]*:" "" _qt_gui_includes "${_qt_gui_includes}")
  string(REPLACE ">" "" _qt_gui_includes "${_qt_gui_includes}")
  if(DEFINED qt_PACKAGE_FOLDER_RELEASE)
    list(APPEND _qt_gui_includes "${qt_PACKAGE_FOLDER_RELEASE}/include" "${qt_PACKAGE_FOLDER_RELEASE}/include/QtGui")
  endif()
  set(_rhi_hints)
  foreach(_dir IN LISTS _qt_gui_includes)
    file(GLOB _found LIST_DIRECTORIES true "${_dir}/*/QtGui" "${_dir}/QtGui/*/QtGui")
    list(APPEND _rhi_hints "${_dir}/${Qt6_VERSION}/QtGui" "${_dir}/QtGui/${Qt6_VERSION}/QtGui" ${_found})
  endforeach()
  find_path(AMCSET_QT_RHI_INCLUDE_DIR rhi/qrhi.h HINTS ${_rhi_hints} NO_DEFAULT_PATH)
  if(NOT AMCSET_QT_RHI_INCLUDE_DIR)
    message(FATAL_ERROR "Qt's QRhi headers (rhi/qrhi.h) were not found; set AMCSET_QT_RHI_INCLUDE_DIR "
                        "to the folder include/QtGui/${Qt6_VERSION}/QtGui of the Qt installation.")
  endif()
  get_filename_component(_rhi_version_dir "${AMCSET_QT_RHI_INCLUDE_DIR}" DIRECTORY)
  target_include_directories(amcset-gui PRIVATE "${AMCSET_QT_RHI_INCLUDE_DIR}" "${_rhi_version_dir}")
endif()

# Shaders: the pre-compiled .qsb files in shaders/baked are used by default,
# because some machines block running qsb.exe (Device Guard). With
# AMCSET_COMPILE_SHADERS=ON they are compiled at build time by Qt Shader Tools.
# Regenerate the baked files with shaders/bake.sh after editing a shader.
option(AMCSET_COMPILE_SHADERS "Compile the GUI shaders with qsb at build time" OFF)
if(AMCSET_COMPILE_SHADERS AND COMMAND qt_add_shaders AND TARGET Qt6::qsb)
  qt_add_shaders(amcset-gui "amcset_shaders"
    PREFIX "/shaders"
    BASE shaders
    GLSL "100 es,120,150" HLSL 50 MSL 12
    FILES shaders/tracks.vert shaders/tracks.frag)
else()
  target_sources(amcset-gui PRIVATE shaders/shaders.qrc)
endif()

if(WIN32)
  # Copy the Qt runtime next to the executable (non-fatal; see
  # gui/scripts/run_windeployqt.cmake).
  get_target_property(QT_QMAKE_EXECUTABLE Qt6::qmake IMPORTED_LOCATION)
  get_filename_component(QT_BIN_DIR "${QT_QMAKE_EXECUTABLE}" DIRECTORY)
  find_program(WINDEPLOYQT_EXECUTABLE windeployqt HINTS "${QT_BIN_DIR}")
  add_custom_command(TARGET amcset-gui POST_BUILD
    COMMAND ${CMAKE_COMMAND} "-DWINDEPLOYQT=${WINDEPLOYQT_EXECUTABLE}"
            "-DTARGET=$<TARGET_FILE:amcset-gui>"
            -P "${PROJECT_SOURCE_DIR}/gui/scripts/run_windeployqt.cmake"
    COMMENT "Running windeployqt..."
    VERBATIM
  )
endif()
