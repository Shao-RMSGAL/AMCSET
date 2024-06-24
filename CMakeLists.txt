cmake_minimum_required(VERSION 3.23)

project(amcset VERSION 1.0)

# Require C++ 20
add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_20)

# Set compiler flags
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
target_compile_options(compiler_flags INTERFACE
  "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused>>"
  )

# Build shared libraries
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Pass version number to compiler
configure_file(common/include/config.h.in config.h)

# Configure subdirectories
add_subdirectory(server)
add_subdirectory(gui)
add_subdirectory(common)

# Packages
find_package(Boost 1.85.0)

# Add package includes and link common
target_include_directories(amcset_common PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(amcset_common PUBLIC ${BOOST_LIBRARY_DIRS})

# Allow link access for common library to gui and server
target_link_libraries(amcset_server PRIVATE amcset_common)
target_link_libraries(amcset_gui PRIVATE amcset_common)
