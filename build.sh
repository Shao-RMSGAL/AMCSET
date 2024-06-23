#!/bin/bash
#
# This script executes the commands needed to configure the build environment
# when the repository is downloaded for the first time. Do not run multiple
# times, unless you have deleted the contents of the /build directory or
# restarted your terminal

conan install .
conan install . -s build_type=Debug

source ./build/Debug/generators/conanbuild.sh

cmake --preset debug
cmake --build --preset debug
