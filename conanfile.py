from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMake, cmake_layout, CMakeDeps


class amcsetRecipe(ConanFile):
    name = "amcset"
    version = "1.0"
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

    def layout(self):
        cmake_layout(self)

    def generate(self):
        deps = CMakeDeps(self)
        deps.generate()
        tc = CMakeToolchain(self)
        tc.user_presets_path = 'ConanPresets.json'
        tc.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        self.conf["tools.system.package_manager:mode"] = "install"
        self.conf["tools.system.package_manager:sudo"] = True

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def requirements(self):
        self.requires("boost/1.85.0")
        self.requires("gtest/1.14.0")
        self.requires("glog/0.7.1")

    def build_requirements(self):
        self.tool_requires("cmake/[>3.23]")
        self.tool_requires("ninja/[>1.12]")
