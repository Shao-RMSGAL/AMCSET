cmake_minimum_required(VERSION 3.23)

project(amcset VERSION 1.0)

add_library(compiler_flags INTERFACE)
target_compile_features(compiler_flags INTERFACE cxx_std_20)

set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")

target_compile_options(compiler_flags INTERFACE
  "$<${gcc_like_cxx}:$<BUILD_INTERFACE:-Wall;-Wextra;-Wshadow;-Wformat=2;-Wunused>>"
  )

option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

configure_file(common/include/config.h.in config.h)

add_subdirectory(server)
add_subdirectory(gui)
add_subdirectory(common)

target_link_libraries(amcset_server PUBLIC amcset_common)
target_link_libraries(amcset_gui PUBLIC amcset_common)
