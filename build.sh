#!/bin/bash
#
# This script executes the commands needed to configure the build environment
# when the repository is downloaded for the first time. Do not run multiple
# times, unless you have deleted the contents of the /build directory or
# restarted your terminal

# Attempt to run conan from sudo. Run as 
# sudo -E ./build.sh
if [ "$EUID" -eq 0 ]; then
    sudo -E $HOME/.local/bin/conan install . --build=missing 
    sudo -E $HOME/.local/bin/conan install . \
        -s build_type=Debug --build=missing
    exit 1
else
    conan install . --build=missing
    conan install . -s build_type=Debug --build=missing
fi
source ./build/Debug/generators/conanbuild.sh

$HOME/Qt/Tools/CMake/bin/cmake -G Ninja \
    --preset debug \
    -DQt6_DIR=$HOME/Qt/6.7.2/gcc_64/lib/cmake/Qt6 \
    -DCMAKE_PREFIX_PATH=$HOME/Qt/6.7.2/gcc_64

cmake --build --preset debug

# /home/nathaniel/Qt/6.7.2/gcc_64/bin/qt-cmake --build --preset debug
# cmake --preset debug

# Commands for release build
# cmake --preset release
# cmake --build --preset release
