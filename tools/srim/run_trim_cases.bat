@echo off
REM Runs run_trim_cases.ps1 (see README.md). Keep this folder inside the SRIM-2013 installation.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_trim_cases.ps1" %*
pause
