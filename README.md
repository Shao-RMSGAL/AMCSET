# Building

Once you have cloned the repository, simply run the build.sh script to build
the program.

```bash
./build.sh
```

# Running on Windows

Run the PowerShell build script from the project root. The installer uses
`winget` for host applications, MSYS2 for the MinGW-w64 compiler, and creates
a dedicated Conan profile that uses GCC rather than MSVC.

```powershell
.\qbuild.ps1 -i -c -b
```

After the environment is installed, the equivalent operations are:

```powershell
.\qbuild.ps1 -c  # Install Conan dependencies
.\qbuild.ps1 -b  # Build Release
.\qbuild.ps1 -t  # Run tests
.\qbuild.ps1 -r  # Run the GUI
```

`-i` requires Microsoft's `winget` App Installer. It installs Python, CMake,
MSYS2, the MinGW-w64 toolchain, and Conan when they are missing.
