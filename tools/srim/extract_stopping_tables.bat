@echo off
REM Runs extract_stopping_tables.ps1 (see README.md). Keep this folder inside the SRIM-2013 installation.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0extract_stopping_tables.ps1" %*
pause
