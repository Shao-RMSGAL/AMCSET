# Building

project_dir="$HOME/Code/C++/AMCSET"
build_dir="$project_dir/build/Debug"
src_dir="$project_dir"
cmake_cmd="$HOME/Qt/Tools/CMake/bin/cmake"
test_dir="$project_dir/build/Debug/test"
cmake_flags="-DCMAKE_EXPORT_COMPILE_COMMANDS=1"
export TESTDIR="$build_dir/log"

build() {
    cd $build_dir
    $cmake_cmd -S $src_dir -B $build_dir $cmake_flags
    $cmake_cmd --build $build_dir --target all 
    cd -
}

test() {
    cd $test_dir
    ctest -V
}

run() {
    $build_dir/amcset-server &
    $build_dir/amcset-gui
}

conan_run() {
    conan install . --build=missing
    conan install . -s build_type=Debug --build=missing
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
