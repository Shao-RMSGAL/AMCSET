# A very simple script to build a release-ready version of AMCSET on Windows.
# You may also run each line individually if preferred. Copy and paste into 
# a Powershell terminal (not a command prompt). Refer to the README.md file
# for instructions on configuring a build environment.
#
# Install dependencies (First run can take ~0.5 hr)
conan install . --build=missing -s build_type=Debug
conan install . --build=missing

# Configure build environment (Fast)
cmake --preset conan-release

# Build AMCSET (May take a minute or so)
cmake --build --preset conan-release

# Copy Qt .dll files into the build directory. (Fast)
C:\Users\natha\.conan2\p\b\qt181826b5d51c2\p\bin\windeployqt.exe .\build\Release\amcset-gui.exe

# Run AMCSET (optional)
.\build\Release\amcset-gui.exe
