# Building

Once you have cloned the repository, simply run the build.sh script to build
the program.

```bash
./build.sh
```

# Running on Windows

1. Install Python

[https://www.python.org/downloads/](https://www.python.org/downloads/)

2. Install CMake (Download __Previous Release (3.31.7)__, the current
version will not play nice with Conan)

[https://cmake.org/download/](https://cmake.org/download/)

2. Install Conan

```powershell
pip install conan
```

3. Conan profile detect

```powershell
conan profile detect
```

4. Change C++ version

Open `C:\Users\natha\.conan2\profiles\default` and change the C++ version to 
17:

```
[settings]
arch=x86_64
build_type=Release
compiler=msvc
compiler.cppstd=14
compiler.runtime=dynamic
compiler.version=194
os=Windows
```

to

```
[settings]
arch=x86_64
build_type=Release
compiler=msvc
compiler.cppstd=17
compiler.runtime=dynamic
compiler.version=194
os=Windows
```
4. Use conan to install dependencies for debug and release mode
    (this will take several minutes)

```powershell
conan install . --build=missing --settings=build_type=Debug
conan install . --build=missing --settings=build_type=Release
```

5. Build the source code

```powershell
.\build\generators\conanbuild.bat
cmake --preset conan-release
```
