# A script for quickly building the AMCSET program (Tested on Linux)
# Instructions:
# 1. Have the following installed:
#   - Conan (https://conan.io/)
#   - Cmake (https://cmake.org/)
# 2. Run `conan profile detect`, then modify ~/.conan2/profiles/default to 
#   use C++23. (Set compiler.cppstd=gnu23)
# 3. Run ./qbuild -c -b from the project root directory

# Building

project_dir="$HOME/Code/C++/AMCSET"
build_dir="$project_dir/build/Debug"
src_dir="$project_dir"
cmake_cmd="cmake"
cmake_flags="--preset conan-debug"
conan_flags="--build=missing --settings=build_type=Debug"
test_dir="$project_dir/build/Debug/test"
export TESTDIR="$build_dir/log"

check_status() {
    if [ $? -ne 0 ]; then 
        echo "An error occured"
        exit 1
    fi
}

build() {
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
    conan install . $conan_flags 
    check_status
    conan install .
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
