cmake_minimum_required(VERSION 3.23)

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
# set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
# set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fsanitize=address")

project(amcset VERSION 0.0.1)

# Use the conanfile.py version number
if(NOT DEFINED PROJECT_VERSION)
    message(FATAL_ERROR "Version components not defined! Are you running through Conan?")
endif()
# For clangd lsp
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Require C++23
add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_20)

# Set compiler flags
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
set(only_clang "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang>")

target_compile_options(compiler_flags INTERFACE
  "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused;>>"
   "$<${only_clang}:$<BUILD_INTERFACE:-fbracket-depth=3500>>"
  "$<${msvc_cxx}:$<BUILD_INTERFACE:-W3>>"
)

# Build shared libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Pass version number to compiler
configure_file(${PROJECT_SOURCE_DIR}/common/include/config.h.in ${PROJECT_SOURCE_DIR}/common/include/config.h)

# Packages
find_package(Boost 1.85.0 REQUIRED)
find_package(GTest 1.14.0 REQUIRED)
find_package(glog 0.7.1 REQUIRED)
find_package(gflags 2.2.2 REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(Qt6 COMPONENTS Core Gui Widgets 3DCore 3DExtras 3DRender 3DInput REQUIRED)
if (Qt6_FOUND)
    qt_standard_project_setup()
endif ()

# Configure subdirectories
add_subdirectory(common)
add_subdirectory(gui)
add_subdirectory(test)

# Add package includes and link common
target_include_directories(amcset-common SYSTEM PUBLIC
    ${Boost_INCLUDE_DIRS}
    )

target_include_directories(amcset-gui PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(amcset-common PUBLIC ${Boost_LIBRARY_DIRS} glog::glog)

# Allow link access for common library to gui
target_link_libraries(amcset-gui PRIVATE 
    amcset-common
    Qt6::Core
    Qt6::Widgets
    Qt6::Gui
    Qt6::3DCore
    Qt6::3DExtras
    Qt6::3DRender
    Qt6::3DInput
    gflags::gflags
)
