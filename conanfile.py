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

from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps


class amcsetRecipe(ConanFile):
    name = "amcset"
    version = "1.0.0"

    package_type = "application"

    license = "GPL-3.0"
    author = "Nathaniel Thomas, nathaniel@swbell.net"
    url = "https://github.com/Shao-RMSGAL/AMCSET"
    description = "AMCSET, a materials ion bombardment simulation program with capabilities for electron bombardment simulation."
    topics = ("simulation", "physics", "particles")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"
    # with_qt=False: Qt comes from elsewhere (the release workflow uses Qt's
    # own binaries; building Qt here takes hours), or no GUI is built.
    options = {"with_qt": [True, False]}
    default_options = {"with_qt": True}

    # Project sources
    # Everything CMakeLists.txt needs (common/ always; cli/, gui/ and test/
    # unless disabled with AMCSET_BUILD_CLI/GUI/TESTS=OFF).
    exports_sources = (
        "CMakeLists.txt",
        "cli/*",
        "gui/*",
        "common/*",
        "test/*",
    )

    def layout(self):
        cmake_layout(self)

    def generate(self):
        # Generate dependency information for CMake.
        deps = CMakeDeps(self)
        deps.generate()

        # Generate CMake toolchain.
        tc = CMakeToolchain(self)

        # Generate ConanPresets.json in the source directory.
        tc.user_presets_path = "ConanPresets.json"

        # Generate compile_commands.json.
        tc.variables["CMAKE_EXPORT_COMPILE_COMMANDS"] = True

        tc.generate()

    def requirements(self):
        self.requires("boost/[=1.89.0]", transitive_headers=True)
        self.requires("gtest/[=1.15.0]")
        self.requires("abseil/[=20260526.0]")
        if self.options.with_qt:
            self.requires("qt/[=6.11.1]")
        # libiconv 1.17 does not build with current MinGW/MSYS2 headers.
        self.requires("libiconv/[=1.18]", override=True)

    def build_requirements(self):
        self.tool_requires("cmake/[=3.29]")
        self.tool_requires("ninja/[=1.12]")

        if self.settings.os == "Windows":
            self.tool_requires("strawberryperl/[=5.32.1.1]")

    def configure(self):
        # ---- Boost ----

        self.options["boost"].without_thread = False
        self.options["boost"].bzip2 = False
        self.options["boost"].zlib = False
        self.options["boost"].without_test = True

        # ---- Qt ----

        if not self.options.with_qt:
            return
        # Use shared Qt libraries.
        self.options["qt"].shared = True
        self.options["qt"].with_pq = False

        # Enable Vulkan support.
        self.options["qt"].with_vulkan = True

        # The GUI needs only Qt Widgets and QRhi (part of Qt GUI); Qt Shader
        # Tools compiles its shaders (pre-compiled copies are used without
        # it). Qt 3D and Qt Graphs are no longer used; they stay enabled only
        # so that an existing Qt package is reused instead of rebuilt from
        # source (hours). Set both to False on a fresh setup.
        self.options["qt"].qt3d = True
        self.options["qt"].qtgraphs = True
        self.options["qt"].qtshadertools = True
        # Qt Quick Controls styles are unused and generate paths beyond Windows MAX_PATH.
        self.options["qt"].disabled_features = (
            "quickcontrols2_basic quickcontrols2_fusion quickcontrols2_fluentwinui3 "
            "quickcontrols2_imagine quickcontrols2_material quickcontrols2_stylekit "
            "quickcontrols2_universal quickcontrols2_windows"
        )

        # Wayland is a Linux build feature. The Windows xkbcommon package
        # does not expose the option required by Qt's Wayland validation.
        self.options["qt"].qtwayland = self.settings.os != "Windows"
        self.options["qt"].with_egl = True

    def build(self):
        cmake = CMake(self)

        cmake.configure()
        cmake.build()

        # Create a convenient compile_commands.json symlink
        # in the source directory. (qbuild.sh copies it instead.)
        import os

        compile_commands_path = os.path.join(
            self.build_folder,
            "compile_commands.json",
        )

        symlink_path = os.path.join(
            self.source_folder,
            "compile_commands.json",
        )

        # Remove an existing file/symlink.
        if os.path.lexists(symlink_path):
            os.remove(symlink_path)

        if os.path.exists(compile_commands_path):
            os.symlink(
                compile_commands_path,
                symlink_path,
            )
        else:
            self.output.warning(
                "compile_commands.json not found in build directory"
            )

    def package(self):
        cmake = CMake(self)
        cmake.install()
