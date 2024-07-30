# Building

project_dir="$HOME/Code/C++/AMCSET"
build_dir="$project_dir/build/Debug"
src_dir="$project_dir"
cmake_cmd="cmake"
cmake_flags="--preset conan-debug"
conan_flags="-s build_type=Debug --build=missing"
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
    $build_dir/amcset-server &
    $build_dir/amcset-gui
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
