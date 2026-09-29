<#
.SYNOPSIS
    A script for quickly building the AMCSET program (Windows equivalent of qbuild)

.DESCRIPTION
    This project uses conan as a package manager, which is _supposed_
    to solve portability issues of dependencies across various platforms.
    Unfortunately, conan is somewhat unreliable in experience, and builds
    will occasionally fail. This includes Windows, though you may need to
    deal with a bit of troubleshooting depending on your toolchain.

    If builds fail, try deleting the conan cache with 'conan remove "*"' or
    'conan cache clean'. Then try to reinstall. You may find that rebooting
    can also solve some problems. Also try adjusting your conan profile if
    needed. Avoid modifying conanfile.py, as this defines the state of
    dependencies, and messing with it may break the intended behavior of
    the code.

.NOTES
    Instructions (Windows):
    1. Have the following installed:
       - Conan (https://conan.io/)
       - CMake
       - A supported compiler toolchain (e.g. MSVC / Visual Studio Build Tools)
    2. Run `conan profile detect`, then modify the generated default profile
       (usually at $env:USERPROFILE\.conan2\profiles\default) to use C++23.
       (Set compiler.cppstd=23)
       Example:

       [settings]
           arch=x86_64
           build_type=Debug
           compiler=msvc
           compiler.cppstd=23
           compiler.runtime=dynamic
           compiler.version=193
           os=Windows
       # Optionally choose your favorite generator, for example, Ninja
       [conf]
           tools.cmake.cmaketoolchain:generator=Ninja

    3. Run .\qbuild.ps1 -Conan -Build -Run from the project root directory
       -Conan tells conan to install all required packages.
       -Build builds the code.
       -Run runs the resulting binary once built.

.PARAMETER Conan
    Install dependencies via conan.

.PARAMETER Build
    Configure and build the project.

.PARAMETER Run
    Run the resulting binary.

.PARAMETER Test
    Run ctest in the test directory.
#>

[CmdletBinding()]
param(
    [Alias("c")][switch]$Conan,
    [Alias("b")][switch]$Build,
    [Alias("r")][switch]$Run,
    [Alias("t")][switch]$Test
)

$ErrorActionPreference = "Continue"
if ($PSVersionTable.PSVersion.Major -ge 7) {
    $PSNativeCommandUseErrorActionPreference = $false
}

# Control whether to run in debug or release mode. Debug seems to be broken.
$mode = "Release"
$modeLower = "release"

$ProjectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$SrcDir = $ProjectDir
$CmakeCmd = "cmake"
$BuildDir = Join-Path $ProjectDir "build\$mode"
$CmakeConfigurePreset = "conan-default"
$CmakeBuildPreset = "conan-release"

$CmakeFlags = @("--preset", "conan-$modeLower")
$RunFlags = @("--minloglevel=3", "--gui=true")

# $TestDir = Join-Path $BuildDir "test"
$TestDir = "build\test"
$env:TESTDIR = Join-Path $BuildDir "log"

function Check-Status {
    param([int]$ExitCode)
    if ($ExitCode -ne 0) {
        Write-Host "An error occurred"
        exit 1
    }
}

function Copy-CompileCommands {
    $src = Join-Path $BuildDir "compile_commands.json"
    if (Test-Path $src) {
        Copy-Item -Path $src -Destination (Join-Path $ProjectDir "compile_commands.json") -Force
        Check-Status $LASTEXITCODE
    } else {
        Write-Host "Warning: compile_commands.json not found in $BuildDir"
    }
}

function Test-ConanBuildFiles {
    $generatorDir = Join-Path $ProjectDir "build\generators"
    $toolchainPath = Join-Path $generatorDir "conan_toolchain.cmake"
    $presetPath = Join-Path $generatorDir "CMakePresets.json"
    $legacyReleasePreset = Join-Path $ProjectDir "build\Release\generators\CMakePresets.json"
    $legacyDebugPreset = Join-Path $ProjectDir "build\Debug\generators\CMakePresets.json"

    if (-not (Test-Path $toolchainPath) -or
        -not (Test-Path $presetPath) -or
        (Test-Path $legacyReleasePreset) -or
        (Test-Path $legacyDebugPreset)) {
        return $false
    }

    $preset = Get-Content -Raw $presetPath
    if ($preset -match "/opt/conan/|/workspace/|/home/|/usr/bin/") {
        return $false
    }

    return $true
}

function Invoke-Build {
    if (-not (Test-ConanBuildFiles)) {
        Write-Host "Conan build files are missing or stale; installing dependencies..."
        Invoke-Conan
    }

    # Enable AddressSanitizer only for debug builds; keep Release fast
    $asanFlag = "OFF"
    if ($modeLower -eq "debug") {
        $asanFlag = "ON"
    }

    & $CmakeCmd --fresh --preset $CmakeConfigurePreset `
        "-DENABLE_ASAN=$asanFlag" `
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON" `
        "-DCMAKE_BUILD_TYPE=$mode" 2>&1 | Out-Host
    Check-Status $LASTEXITCODE

    & $CmakeCmd --build --preset $CmakeBuildPreset 2>&1 | Out-Host
    Check-Status $LASTEXITCODE

    Copy-CompileCommands
}

function Invoke-Test {
    Push-Location $TestDir
    try {
        ctest -V
        Check-Status $LASTEXITCODE
    } finally {
        Pop-Location
    }
}

function Invoke-Run {
    $exe = Join-Path $BuildDir "amcset-gui.exe"
    & $exe @RunFlags
    Check-Status $LASTEXITCODE
}

function Invoke-Conan {
    Push-Location $ProjectDir
    try {
        Remove-Item -Recurse -Force (Join-Path $ProjectDir "build") `
            -ErrorAction SilentlyContinue

        & conan install . `
            "--settings=build_type=$mode" `
            "--settings=compiler.cppstd=20" `
            "--build=missing" 2>&1 | Out-Host
        Check-Status $LASTEXITCODE
    } finally {
        Pop-Location
    }
}

# Note: On Windows, compiler.cstd is generally not set the same way as on
# Linux/gcc; MSVC's C standard is tied to the compiler version. Adjust the
# conan settings above (and your conan profile) to match your toolchain.

if (-not ($Conan -or $Build -or $Run -or $Test)) {
    Invoke-Build
} else {
    if ($Conan) { Invoke-Conan }
    if ($Build) { Invoke-Build }
    if ($Test)  { Invoke-Test }
    if ($Run)   { Invoke-Run }
}