@echo off
REM Speed of AMCSET's SRIM mode against TRIM.exe (see benchmark_vs_trim.ps1).
REM     tools\windows\benchmark_vs_trim.bat [SRIM folder]
REM Close other programs first; results in build\benchmark\benchmark.csv.
setlocal
cd /d "%~dp0\..\.."
call build\Release\generators\conanrun.bat >NUL 2>&1
if "%~1"=="" (
  powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0benchmark_vs_trim.ps1" > build\benchmark_log.txt 2>&1
) else (
  powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0benchmark_vs_trim.ps1" -Srim "%~1" > build\benchmark_log.txt 2>&1
)
echo %ERRORLEVEL%> build\benchmark_exit.txt
endlocal
