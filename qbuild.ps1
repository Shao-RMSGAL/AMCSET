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
    2. Run `conan profile detect`. The -Conan step overrides
       compiler.cppstd=20 and build_type=Release on the command line, so the
       profile (usually at $env:USERPROFILE\.conan2\profiles\default) mainly
       needs the right compiler. Example:

       [settings]
           arch=x86_64
           build_type=Release
           compiler=msvc
           compiler.cppstd=20
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

# Only Release is supported (see qbuild.sh).
$mode = "Release"
$modeLower = "release"

$ProjectDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$CmakeCmd = $null
$BuildDir = Join-Path $ProjectDir "build\$mode"
$CmakeConfigurePreset = "conan-release"
$CmakeBuildPreset = "conan-release"

$RunFlags = @("--minloglevel=3", "--gui=true")

# NOTE: This differs from qbuild.sh (build/Release/test) and from $BuildDir
# above; the CI workflow (.github/workflows/tests.yml) also uses build/test.
# Unverified which one the Windows toolchain actually produces -- see
# docs/KNOWN_ISSUES.md. (Was a path relative to the current directory; it is
# now anchored at the repository root, which is identical when the script is
# run from the root as documented.)
$TestDir = Join-Path $ProjectDir "build\test"

function Check-Status {
    param([int]$ExitCode)
    if ($ExitCode -ne 0) {
        Write-Host "An error occurred"
        exit 1
    }
}

function Test-ConanBuildFiles {
    $generatorDir = Join-Path $BuildDir "generators"
    $toolchainPath = Join-Path $generatorDir "conan_toolchain.cmake"
    $presetPath = Join-Path $generatorDir "CMakePresets.json"
    $buildEnvPath = Join-Path $generatorDir "conanbuild.bat"

    if (-not (Test-Path $toolchainPath) -or
        -not (Test-Path $presetPath) -or
        -not (Test-Path $buildEnvPath)) {
        return $false
    }

    $presets = Get-Content -Raw $presetPath | ConvertFrom-Json
    $configurePreset = $presets.configurePresets |
        Where-Object { $_.name -eq "conan-$modeLower" } |
        Select-Object -First 1
    if (-not $configurePreset -or
        $configurePreset.toolchainFile -match "/opt/conan/|/workspace/|/home/|/usr/bin/") {
        return $false
    }

    return $true
}

function Invoke-Build {
    Push-Location $ProjectDir
    try {
        if (-not (Test-ConanBuildFiles)) {
            Write-Host "Conan build files are missing or stale; installing dependencies..."
            Invoke-Conan
        }

        $presetPath = Join-Path $BuildDir "generators\CMakePresets.json"
        $presets = Get-Content -Raw $presetPath | ConvertFrom-Json
        $configurePreset = $presets.configurePresets |
            Where-Object { $_.name -eq "conan-$modeLower" } |
            Select-Object -First 1
        $CmakeCmd = $configurePreset.cmakeExecutable
        if (-not $CmakeCmd -or -not (Test-Path $CmakeCmd)) {
            $cmakeCommand = Get-Command cmake -ErrorAction SilentlyContinue
            if (-not $cmakeCommand) {
                throw "CMake executable was not found in the Conan preset or PATH."
            }
            $CmakeCmd = $cmakeCommand.Source
        }

        # Enable AddressSanitizer only for debug builds; keep Release fast
        $asanFlag = "OFF"
        if ($modeLower -eq "debug") {
            $asanFlag = "ON"
        }

        $buildEnvPath = Join-Path $BuildDir "generators\conanbuild.bat"
        $buildCommand = 'call "{0}" && "{1}" --fresh --preset {2} -DENABLE_ASAN={3} -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DCMAKE_BUILD_TYPE={4} && "{1}" --build --preset {5}' -f `
            $buildEnvPath, $CmakeCmd, $CmakeConfigurePreset, $asanFlag, $mode, $CmakeBuildPreset
        & $env:ComSpec /d /s /c $buildCommand 2>&1 | Out-Host
        Check-Status $LASTEXITCODE

        Copy-CompileCommands
    } finally {
        Pop-Location
    }
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
            "--conf=tools.cmake.cmaketoolchain:generator=Ninja" `
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