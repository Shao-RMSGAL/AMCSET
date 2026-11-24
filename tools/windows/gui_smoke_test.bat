@echo off
REM Smoke test of the GUI: start amcset-gui, run the default scenario (or the
REM file given as argument), save images of all views and quit. Each run gets
REM its own folder build\gui_smoke\<n> with the images, log.txt, status.txt and
REM exit_code.txt, which amcset-gui writes when it has closed its window
REM (0 = the run finished; status.txt also ends with the code when cmd gets
REM to write it).
setlocal
cd /d "%~dp0\..\.."
set "OUT=%CD%\build\gui_smoke\%RANDOM%"
mkdir "%OUT%"
echo %OUT%> "%CD%\build\gui_smoke_last.txt"
echo RUNNING> "%OUT%\status.txt"
call build\Release\generators\conanrun.bat >NUL 2>&1
set QT_FORCE_STDERR_LOGGING=1
build\Release\amcset-gui.exe --screenshots "%OUT%" --quit --timeout 600 %* > "%OUT%\log.txt" 2>&1
set "RC=%ERRORLEVEL%"
REM The build folder may be synchronised (Nextcloud, OneDrive), which can hold
REM status.txt for a moment: retry the write.
for /L %%i in (1,1,20) do (
  (echo %RC%>> "%OUT%\status.txt") 2>NUL && goto written
  ping -n 2 127.0.0.1 >NUL
)
:written
endlocal
