# A script for quickly building the AMCSET program (Tested on Linux)
# Instructions:
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
# 3. Run ./qbuild -c -b from the project root directory

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
    $build_dir/amcset-server
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
