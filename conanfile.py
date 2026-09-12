from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps


class amcsetRecipe(ConanFile):
    name = "amcset"
    version = "0.0.1"

    package_type = "application"

    license = "GPL-3.0"
    author = "Nathaniel Thomas, nathaniel@swbell.net"
    url = "https://github.com/Shao-RMSGAL/AMCSET"
    description = "AMCSET, a materials ion bombardment simulation program with capabilities for electron bombardment simulation."
    topics = ("simulation", "physics", "particles")

    # Binary configuration
    settings = "os", "compiler", "build_type", "arch"

    # Project sources
    exports_sources = (
        "CMakeLists.txt",
        "gui/*",
        "common/*",
        "server/*",
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
        self.requires("boost/[=1.86.0]", transitive_headers=True)
        self.requires("gtest/[=1.15.0]")
        self.requires("gflags/[=2.2.2]")
        self.requires("glog/[=0.7.1]")
        self.requires("qt/[=6.11.1]")

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

        # Use shared Qt libraries.
        self.options["qt"].shared = True

        # Enable Vulkan support.
        self.options["qt"].with_vulkan = True

        # Enable Qt 3D and Graphs for 2D/3D plotting.
        self.options["qt"].qt3d = True
        self.options["qt"].qtgraphs = True
        self.options["qt"].qtshadertools = True

        # Enable Wayland support.
        #
        # This is a build-time feature, not something that should
        # depend on whether DISPLAY/WAYLAND_DISPLAY happens to exist
        # when Conan is executed.
        self.options["qt"].qtwayland = True
        self.options["qt"].with_egl = True

    def build(self):
        cmake = CMake(self)

        cmake.configure()
        cmake.build()

        # Create a convenient compile_commands.json symlink
        # in the source directory.
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
