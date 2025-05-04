from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps
from conan.tools.files import load
import os
import sys


def display_server():
    if "WAYLAND_DISPLAY" in os.environ:
        return "Wayland"
    elif "DISPLAY" in os.environ:
        return "X11"
    elif sys.platform.startswith("win"):
        return "Windows"
    else:
        raise Exception("Unable to determine display server")


class amcsetRecipe(ConanFile):
    name = "amcset"
    package_type = "application"
    # Optional metadata
    license = "GPL 3.0"
    author = "Nathaniel Thomas, nathaniel@swbell.net"
    url = "https://github.com/Shao-RMSGAL/AMCSET"
    description = "This is AMCSET, a materials ion bombardment simulation"
    "program with capabilities for electron bombardment simulation."

    topics = ("simulation", "physics", "particles")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"

    # Sources are located in the same place as this recipe, copy them to the
    # recipe
    exports_sources = "CMakeLists.txt", "gui/src", "common/src", "server/src"

    def set_version(self):
        self.version = load(self, os.path.join(self.recipe_folder, "version.txt"))

    def layout(self):
        cmake_layout(self)

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()
        tc = CMakeToolchain(self)
        tc.user_presets_path = "ConanPresets.json"
        tc.variables["CMAKE_EXPORT_COMPILE_COMMANDS"] = True
        tc.cache_variables["ENABLE_ASAN"] = "OFF"
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        #  self.conf["tools.system.package_manager:mode"] = "install"
        #  self.conf["tools.system.package_manager:sudo"] = True
        src_dir = self.source_folder
        build_dir = self.build_folder
        src_dir = self.source_folder
        build_dir = self.build_folder
        compile_commands_path = os.path.join(build_dir, "compile_commands.json")
        symlink_path = os.path.join(src_dir, "compile_commands.json")

        # Remove existing symlink if it exists
        if os.path.exists(symlink_path):
            os.remove(symlink_path)

        # Create new symlink
        if os.path.exists(compile_commands_path):
            os.symlink(compile_commands_path, symlink_path)
        else:
            self.output.warning("compile_commands.json not found in build directory")

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def requirements(self):
        self.requires("boost/[=1.85.0]")
        self.requires("gtest/[=1.16.0]")
        self.requires("gflags/[=2.2.2]")
        self.requires("glog/[=0.7.1]")
        self.requires("qt/[=6.7.3]")

    def build_requirements(self):
        self.tool_requires("cmake/[=3.29]")
        self.tool_requires("ninja/[=1.12]")
        if self.settings.os == "Windows":
            # Required for Qt on Windows
            self.tool_requires("strawberryperl/[=5.32.1.1]")

        # Boost settings

    def configure(self):
        boost = self.options["boost"]

        boost.shared = True  # or False, but be consistent
        boost.without_thread = False
        boost.bzip2 = False
        boost.zlib = False
        boost.transitive_headers = False

        # Qt settings
        qt = self.options["qt"]

        # Modules
        qt.qt3d = True
        qt.qtdatavis3d = True
        qt.qtgraphs = True
        qt.qtcharts = True
        qt.qtshadertools = True
        qt.qtdeclarative = True

        # Options
        qt.shared = True
        qt.gui = True
        qt.widgets = True

        #  qt.with_vulkan = True
        #  qt.with_mysql = True
        #  qt.with_glib = True
        #  qt.with_libjpeg = "libjpeg-turbo"
        #  qt.with_zstd = True
        #  qt.with_dbus = True
        #  qt.with_libalsa = True
        #  qt.with_gstreamer = True
        #  qt.with_pulseaudio = True
        #  qt.with_gssapi = False
        #  qt.opengl = "desktop"

        qt.with_fontconfig = True  # For font handling
        qt.with_harfbuzz = True  # For text rendering
        qt.with_freetype = True  # For font rendering

        # Enable required Qt modules
        if display_server() == "Wayland":
            self.options["qt"].qtwayland = True
        elif display_server() in ["X11", "Windows"]:
            self.options["qt"].qtwayland = False
