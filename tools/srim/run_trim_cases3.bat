@echo off
REM Runs TRIM reference batch 3 (trim_cases3 -> results\trim3), see run_trim_cases3.ps1.
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0run_trim_cases3.ps1"
pause
