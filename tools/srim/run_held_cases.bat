@echo off
REM Runs the held TRIM.DAT recoil cases (codes 6/7); answer TRIM's message box in each window.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_held_cases.ps1"
pause
