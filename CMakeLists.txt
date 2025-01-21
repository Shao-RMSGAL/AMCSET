cmake_minimum_required(VERSION 3.23)
project(amcset VERSION 1.0)

# For clangd lsp
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Require C++20
add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_23)

# Set compiler flags
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
target_compile_options(compiler_flags INTERFACE
  "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused>>"
  "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
)

# Build shared libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Pass version number to compiler
configure_file(common/include/config.h.in config.h)


# Packages
find_package(Boost 1.85.0 REQUIRED)
find_package(GTest 1.14.0 REQUIRED)
find_package(glog 0.7.1 REQUIRED)
find_package(gflags 2.2.2 REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(Qt6 COMPONENTS Core Gui Widgets REQUIRED)
if (Qt6_FOUND)
    qt_standard_project_setup()
endif ()

# Configure subdirectories
add_subdirectory(common)
add_subdirectory(server)
add_subdirectory(test)

# Add package includes and link common
target_include_directories(amcset-common SYSTEM PUBLIC 
    ${Boost_INCLUDE_DIRS}
    )

# target_include_directories(amcset-gui PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(amcset-common PUBLIC ${Boost_LIBRARY_DIRS} glog::glog)

# Allow link access for common library to server
target_link_libraries(amcset-server PRIVATE amcset-common 
    Qt6::Core
    Qt6::Widgets
    Qt6::Gui
    glog::glog
    gflags::gflags
)
