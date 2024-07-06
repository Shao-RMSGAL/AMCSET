cmake_minimum_required(VERSION 3.23)
project(amcset VERSION 1.0)

# For clangd lsp
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Require C++20
add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_20)

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

# Configure subdirectories
add_subdirectory(common)
add_subdirectory(server)
add_subdirectory(gui)
add_subdirectory(test)

# Add package includes and link common
target_include_directories(amcset-common SYSTEM PUBLIC 
    ${Boost_INCLUDE_DIRS}
    glog::glog
    )

target_include_directories(amcset-gui PRIVATE ${Boost_INCLUDE_DIRS})
target_link_libraries(amcset-common PUBLIC ${Boost_LIBRARY_DIRS})

# Allow link access for common library to server
target_link_libraries(amcset-server PRIVATE amcset-common)
target_link_libraries(amcset-gui PRIVATE amcset-common)

# Install directories (When running cmake --install)
install(TARGETS amcset-gui
    BUNDLE  DESTINATION .
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

qt_generate_deploy_app_script(
    TARGET amcset-gui
    OUTPUT_SCRIPT deploy_script
    NO_UNSUPPORTED_PLATFORM_ERROR
)
install(SCRIPT ${deploy_script})
