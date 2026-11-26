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

# Installation and packages (docs/RELEASING.md).
#
# `cmake --install <build> --prefix <dir>` lays out:
#   bin/                      amcset, amcset-gui (Windows: and the Qt and C++ runtime DLLs)
#   share/amcset/srim2013/    the stopping tables (common/CMakeLists.txt)
#   share/doc/amcset/         README, CHANGELOG, LICENSE, the user documents
#   Linux only: share/applications, share/metainfo, share/icons/hicolor, share/man/man1
# `cpack` in the build folder packages that: ZIP and MSI (WiX, when found) on
# Windows; TGZ, DEB and RPM (rpmbuild, when found) on Linux.

set(AMCSET_APP_ID "io.github.Shao_RMSGAL.AMCSET")
set(AMCSET_PACKAGE_CONTACT "AMCSET maintainers <https://github.com/Shao-RMSGAL/AMCSET/issues>" CACHE STRING
    "Maintainer of the Linux packages (Debian's Maintainer field: name <e-mail>)")
string(TIMESTAMP AMCSET_RELEASE_DATE "%Y-%m-%d" UTC)
set(AMCSET_RELEASE_DATE "${AMCSET_RELEASE_DATE}" CACHE STRING "Release date in the AppStream metadata")

# --- Documents ----------------------------------------------------------------
set(docdir "${CMAKE_INSTALL_DATADIR}/doc/amcset")
install(FILES "${PROJECT_SOURCE_DIR}/README.md" "${PROJECT_SOURCE_DIR}/CHANGELOG.md" "${PROJECT_SOURCE_DIR}/LICENSE"
        DESTINATION "${docdir}")
install(FILES
          "${PROJECT_SOURCE_DIR}/docs/README.md" "${PROJECT_SOURCE_DIR}/docs/INSTALL.md"
          "${PROJECT_SOURCE_DIR}/docs/CLI.md" "${PROJECT_SOURCE_DIR}/docs/GUI.md"
          "${PROJECT_SOURCE_DIR}/docs/OUTPUT.md" "${PROJECT_SOURCE_DIR}/docs/SRIM_MODE.md"
          "${PROJECT_SOURCE_DIR}/docs/GLOSSARY.md" "${PROJECT_SOURCE_DIR}/docs/KNOWN_ISSUES.md"
        DESTINATION "${docdir}/docs")

# --- Linux desktop integration ---------------------------------------------------
if(UNIX AND NOT APPLE)
  # Manual pages, gzip-compressed as distributions install them.
  foreach(page amcset amcset-gui)
    configure_file("${PROJECT_SOURCE_DIR}/packaging/linux/${page}.1.in" "${CMAKE_CURRENT_BINARY_DIR}/${page}.1" @ONLY)
    file(ARCHIVE_CREATE OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${page}.1.gz" PATHS "${CMAKE_CURRENT_BINARY_DIR}/${page}.1"
         FORMAT raw COMPRESSION GZip)
  endforeach()
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/amcset.1.gz" DESTINATION "${CMAKE_INSTALL_MANDIR}/man1")
  if(AMCSET_BUILD_GUI)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/amcset-gui.1.gz" DESTINATION "${CMAKE_INSTALL_MANDIR}/man1")
    install(FILES "${PROJECT_SOURCE_DIR}/packaging/linux/${AMCSET_APP_ID}.desktop"
            DESTINATION "${CMAKE_INSTALL_DATADIR}/applications")
    configure_file("${PROJECT_SOURCE_DIR}/packaging/linux/${AMCSET_APP_ID}.metainfo.xml.in"
                   "${CMAKE_CURRENT_BINARY_DIR}/${AMCSET_APP_ID}.metainfo.xml" @ONLY)
    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${AMCSET_APP_ID}.metainfo.xml"
            DESTINATION "${CMAKE_INSTALL_DATADIR}/metainfo")
    foreach(size 16 22 24 32 48 64 128 256 512)
      install(FILES "${PROJECT_SOURCE_DIR}/resources/icons/amcset-${size}.png"
              DESTINATION "${CMAKE_INSTALL_DATADIR}/icons/hicolor/${size}x${size}/apps"
              RENAME "${AMCSET_APP_ID}.png")
    endforeach()
    install(FILES "${PROJECT_SOURCE_DIR}/resources/icons/amcset.svg"
            DESTINATION "${CMAKE_INSTALL_DATADIR}/icons/hicolor/scalable/apps"
            RENAME "${AMCSET_APP_ID}.svg")
  endif()
  # Debian's copyright file (machine-readable format).
  install(FILES "${PROJECT_SOURCE_DIR}/packaging/linux/copyright"
          DESTINATION "${docdir}")
endif()

# --- Windows runtime ---------------------------------------------------------------
if(WIN32)
  # The Visual C++ runtime next to the executables (app-local: no
  # redistributable and no administrator rights needed for the ZIP).
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_DESTINATION "${CMAKE_INSTALL_BINDIR}")
  set(CMAKE_INSTALL_UCRT_LIBRARIES FALSE)
  set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON) # (an older CMake may not know this compiler; see below)
  include(InstallRequiredSystemLibraries)
  if(NOT CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS)
    # Older CMake without this compiler's redistributable: take it from the
    # Visual Studio installation the compiler belongs to.
    get_filename_component(vc_dir "${CMAKE_CXX_COMPILER}" DIRECTORY)
    file(GLOB crt_dirs "${vc_dir}/../../../../../../Redist/MSVC/*/x64/Microsoft.VC*.CRT")
    list(SORT crt_dirs)
    list(POP_BACK crt_dirs crt)
    if(crt)
      file(GLOB crt_dlls "${crt}/*.dll")
      install(FILES ${crt_dlls} DESTINATION "${CMAKE_INSTALL_BINDIR}")
    else()
      message(WARNING "The Visual C++ runtime was not found; installed programs need the redistributable")
    endif()
  endif()
  if(AMCSET_BUILD_GUI)
    # The Qt runtime that windeployqt put next to amcset-gui.exe (gui/CMakeLists.txt):
    # only what the GUI uses.
    set(qt_dir "${PROJECT_BINARY_DIR}")
    install(FILES "${qt_dir}/Qt6Core.dll" "${qt_dir}/Qt6Gui.dll" "${qt_dir}/Qt6Widgets.dll"
            DESTINATION "${CMAKE_INSTALL_BINDIR}")
    install(FILES "${qt_dir}/d3dcompiler_47.dll" DESTINATION "${CMAKE_INSTALL_BINDIR}" OPTIONAL)
    install(FILES "${qt_dir}/platforms/qwindows.dll"
            DESTINATION "${CMAKE_INSTALL_BINDIR}/platforms")
    install(FILES "${qt_dir}/styles/qmodernwindowsstyle.dll"
            DESTINATION "${CMAKE_INSTALL_BINDIR}/styles" OPTIONAL)
    install(FILES "${qt_dir}/imageformats/qico.dll" "${qt_dir}/imageformats/qgif.dll"
            DESTINATION "${CMAKE_INSTALL_BINDIR}/imageformats" OPTIONAL)
  endif()
endif()

# --- Packages ------------------------------------------------------------------------
set(CPACK_PACKAGE_NAME "AMCSET")
set(CPACK_PACKAGE_VENDOR "Texas A&M University")
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Monte Carlo simulation of ions in matter and the radiation damage they cause")
set(CPACK_PACKAGE_DESCRIPTION
    "AMCSET simulates the transport of ions in layered amorphous and crystalline targets and the radiation damage they cause. Its SRIM mode reproduces SRIM-2013's TRIM statistically (TRIM.IN input, TRIM output files) and runs on all processor cores. A graphical program (amcset-gui) and a command-line program (amcset). The bundled SRIM-2013 stopping tables are derived from SRIM (J. F. Ziegler) and are free for non-commercial use only.")
set(CPACK_PACKAGE_HOMEPAGE_URL "${PROJECT_HOMEPAGE_URL}")
set(CPACK_PACKAGE_CONTACT "${AMCSET_PACKAGE_CONTACT}")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "AMCSET")
set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/resources/icons/amcset.ico")
set(CPACK_STRIP_FILES ON)
set(CPACK_SOURCE_IGNORE_FILES "/build.*/;/\\\\.git/;/output/;\\\\.conflict")
# The license as a .txt (WiX needs an extension it knows).
configure_file("${PROJECT_SOURCE_DIR}/LICENSE" "${CMAKE_CURRENT_BINARY_DIR}/LICENSE.txt" COPYONLY)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_BINARY_DIR}/LICENSE.txt")
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(AMD64|x86_64|amd64)$")
  set(amcset_arch x86_64)
else()
  set(amcset_arch "${CMAKE_SYSTEM_PROCESSOR}")
endif()

if(WIN32)
  set(CPACK_PACKAGE_FILE_NAME "AMCSET-${PROJECT_VERSION}-windows-x64")
  set(CPACK_GENERATOR ZIP)
  # MSI (WiX Toolset; v4 and later are .NET tools: `dotnet tool install --global wix`).
  find_program(AMCSET_WIX wix)
  find_program(AMCSET_WIX3_CANDLE candle PATHS "$ENV{WIX}/bin")
  if(AMCSET_WIX OR AMCSET_WIX3_CANDLE)
    list(APPEND CPACK_GENERATOR WIX)
  endif()
  if(AMCSET_WIX)
    set(CPACK_WIX_VERSION 4)
  endif()
  # Never change: Windows recognises later versions of the product by it.
  set(CPACK_WIX_UPGRADE_GUID "6F3C2B5E-8A41-4C7D-9E2B-1D5A7C3E9F40")
  set(CPACK_WIX_PRODUCT_ICON "${PROJECT_SOURCE_DIR}/resources/icons/amcset.ico")
  set(CPACK_WIX_PROGRAM_MENU_FOLDER "AMCSET")
  set(CPACK_WIX_ARCHITECTURE x64)
  set(CPACK_WIX_INSTALL_SCOPE perMachine) # Program Files, for every user
  set(CPACK_WIX_PROPERTY_ARPURLINFOABOUT "${PROJECT_HOMEPAGE_URL}")
  set(CPACK_WIX_PROPERTY_ARPHELPLINK "${PROJECT_HOMEPAGE_URL}/issues")
  set(CPACK_WIX_PROPERTY_ARPCONTACT "${PROJECT_HOMEPAGE_URL}/issues")
  # amcset on the (system) PATH, and a Start menu shortcut of the GUI.
  set(CPACK_WIX_PATCH_FILE "${PROJECT_SOURCE_DIR}/packaging/windows/wix_patch.xml")
  if(AMCSET_BUILD_GUI)
    set_property(INSTALL "${CMAKE_INSTALL_BINDIR}/amcset-gui.exe" PROPERTY CPACK_START_MENU_SHORTCUTS "AMCSET")
  endif()
elseif(UNIX AND NOT APPLE)
  set(CPACK_PACKAGE_FILE_NAME "amcset-${PROJECT_VERSION}-linux-${amcset_arch}")
  set(CPACK_GENERATOR TGZ)
  find_program(AMCSET_DPKG dpkg-deb)
  find_program(AMCSET_RPMBUILD rpmbuild)
  if(AMCSET_DPKG)
    list(APPEND CPACK_GENERATOR DEB)
  endif()
  if(AMCSET_RPMBUILD)
    list(APPEND CPACK_GENERATOR RPM)
  endif()
  # (DEB and RPM install to /usr; the tarball holds bin/ and share/ at its top.)
  # Debian
  set(CPACK_DEBIAN_PACKAGE_NAME "amcset")
  set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
  set(CPACK_DEBIAN_PACKAGE_SECTION "science")
  set(CPACK_DEBIAN_PACKAGE_PRIORITY "optional")
  set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  set(CPACK_DEBIAN_PACKAGE_HOMEPAGE "${PROJECT_HOMEPAGE_URL}")
  # RPM
  set(CPACK_RPM_PACKAGE_NAME "amcset")
  set(CPACK_RPM_FILE_NAME RPM-DEFAULT)
  set(CPACK_RPM_PACKAGE_RELEASE_DIST ON)
  set(CPACK_RPM_PACKAGE_LICENSE "GPL-3.0-or-later AND LicenseRef-SRIM-noncommercial")
  set(CPACK_RPM_PACKAGE_GROUP "Applications/Engineering")
  set(CPACK_RPM_PACKAGE_URL "${PROJECT_HOMEPAGE_URL}")
  set(CPACK_RPM_PACKAGE_AUTOREQ ON)
  set(CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST_ADDITION
      /usr/share/applications /usr/share/metainfo /usr/share/icons /usr/share/icons/hicolor
      /usr/share/man /usr/share/man/man1)
endif()

include(CPack)
