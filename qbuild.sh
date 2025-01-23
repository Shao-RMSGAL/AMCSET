# A script for quickly building the AMCSET program (Tested on Linux)
#
# This project uses conan as a package manager, which is _supposed_
# to solve portability issues of dependancies across various platform.
# Unfortunately, conan is somewhat unreliable in my experience, 
# and builds will occasionally fail. Your best bet is use to a 
# standard, clean Ubuntu distribution to build this project, 
# but if you're willing to deal with a bit of troubleshooting,
# this project should build anywhere that the project dependancies
# are supported. This includes Windows.
#
# If builds fail, try deleting conan cache with 'conan remove '*'' or 
# `conan cache clean`. Then try to reinstall. You may find that rebooting
# can also solve some problems. Also try to make adjustments to your
# conan profile if needed. Avoid modifying the conanfile.py, as this defines
# the state of dependancies, and messing with them may mess with the intended
# behavior of the code.
#
# Instructions (Linux):
# 1. Have the following installed:
#   - Conan (https://conan.io/)
# 2. Run `conan profile detect`, then modify ~/.conan2/profiles/default to 
#   use C++23. (Set compiler.cppstd=gnu23)
#   Here is an example:
#
#   [settings]
#       arch=x86_64
#       build_type=Debug
#       compiler=gcc
#       compiler.cppstd=gnu23
#       compiler.libcxx=libstdc++11
#       compiler.version=14
#       os=Linux
#   # Optionally choose your favorite generator, for example, Ninja
#   [conf]
#       tools.cmake.cmaketoolchain:generator=Ninja
# 3. Run ./qbuild -c -b -r from the project root directory
#   -c tells conan to install all required packages.
#   -b builds the code
#   -r runs the resulting binary once built.

# Building

# Control whether to run in debug or release mode. Debug seems to be broken.
mode="Debug"
mode_lower="debug"
project_dir="$HOME/Code/C++/AMCSET"
src_dir="$project_dir"
cmake_cmd="cmake"
build_dir="$project_dir/build/$mode"

conan_flags="--build=missing --settings=build_type=$mode"
cmake_flags="--preset conan-$mode_lower"
run_flags="-v=3"

test_dir="$project_dir/build/$mode/test"
export TESTDIR="$build_dir/log"

check_status() {
    if [ $? -ne 0 ]; then 
        echo "An error occured"
        exit 1
    fi
}

build() {
    source ./build/$mode/generators/conanbuild.sh
    $cmake_cmd $cmake_flags
    check_status
    $cmake_cmd --build $cmake_flags
    check_status
    cd -
}

test() {
    cd $test_dir
    ctest -V
    check_status
}

run() {
    $build_dir/amcset-server $run_flags
    check_status
}

conan_run() {
    conan install . --settings=build_type=Debug --build=missing
    check_status
    conan install . --settings=build_type=Release --build=missing
    check_status
}

if [ "$#" -eq 0 ]; then
    build
else
    for arg in "$@"
    do
        if [ "$arg" == "--test" ] ||  [ "$arg" == "-t" ]; then 
            test
        fi 

        if [ "$arg" == "--run" ] || [ "$arg" == "-r" ]; then
            run
        fi

        if [ "$arg" == "--conan" ] || [ "$arg" == "-c" ]; then
            conan_run
        fi

        if [ "$arg" == "--build" ] || [ "$arg" == "-b" ]; then
            build
        fi
    done
fi
