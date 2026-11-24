@echo off
REM Configure, build and test AMCSET on Windows (MSVC + Conan).
REM
REM Double-click it, or run it from a command prompt:
REM     tools\windows\build_and_test.bat            build and test
REM     tools\windows\build_and_test.bat --conan    (re)install the Conan dependencies first
REM
REM The Conan dependencies (Qt, Boost, abseil, GoogleTest) are installed
REM with "conan install" only when build\Release\generators does not exist yet or
REM --conan is given; building Qt from source takes hours the first time.
REM Output goes to build\build_and_test.log; build\build_and_test.status ends
REM with SUCCESS, BUILD FAILED or TESTS FAILED.
setlocal
cd /d "%~dp0\..\.."
if not exist build mkdir build
set "LOG=%CD%\build\build_and_test.log"
set "STATUS=%CD%\build\build_and_test.status"
echo RUNNING> "%STATUS%"
echo AMCSET build started %date% %time%> "%LOG%"

if "%~1"=="--conan" goto conan
if exist build\Release\generators\conanbuild.bat goto configure
:conan
echo Installing Conan dependencies...
conan install . --build=missing -s build_type=Release -s compiler.cppstd=20 >> "%LOG%" 2>&1 || goto buildfail

:configure
call build\Release\generators\conanbuild.bat >> "%LOG%" 2>&1
echo Configuring...
cmake --preset conan-release -DAMCSET_BUILD_GUI=ON -DAMCSET_BUILD_CLI=ON -DAMCSET_BUILD_TESTS=ON >> "%LOG%" 2>&1 || goto buildfail
echo Building...
cmake --build --preset conan-release -- -k 0 >> "%LOG%" 2>&1 || goto buildfail
call build\Release\generators\conanrun.bat >> "%LOG%" 2>&1
REM Windows may refuse to start a newly built unsigned program (Device Guard /
REM Smart App Control); cmd then prints the reason here.
echo Checking that amcset.exe starts...>> "%LOG%"
build\Release\amcset.exe version >> "%LOG%" 2>&1 || echo amcset.exe did not start or failed: see above. Windows Security, Protection history lists programs it blocked.>> "%LOG%"
echo Testing...
ctest --preset conan-release --output-on-failure >> "%LOG%" 2>&1 || goto testfail
echo SUCCESS>> "%STATUS%"
echo Build and tests succeeded. See %LOG%
goto end

:buildfail
echo BUILD FAILED>> "%STATUS%"
echo Build failed. See %LOG%
goto end

:testfail
echo TESTS FAILED>> "%STATUS%"
echo Some tests failed. See %LOG%

:end
echo AMCSET build finished %date% %time%>> "%LOG%"
REM Keep the window open for two minutes (any key closes it).
if "%AMCSET_NO_PAUSE%"=="" timeout /t 120
endlocal
