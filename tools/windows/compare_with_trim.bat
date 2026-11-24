@echo off
REM Compare AMCSET with TRIM reference runs (tools\trim_compare\batch.py).
REM
REM     tools\windows\compare_with_trim.bat [SRIM folder] [TRIM cases folder]
REM
REM The SRIM folder (default %USERPROFILE%\Applications\SRIM-2013-Pro) holds the
REM TRIM reference runs; AMCSET uses its default stopping, the SRIM-2013 tables in
REM data\srim2013 (set STOPPING_TABLE to compare with another file). The cases folder
REM (default <SRIM>\_amcset\results\trim3) holds one finished TRIM run per folder.
REM AMCSET runs every finished case with several seeds into
REM build\trim_validation\<cases folder name>\ and compare.py writes cmp\*.json;
REM re-running skips finished runs. The statistics go to analysis.txt/.json.
setlocal
cd /d "%~dp0\..\.."
set "SRIM=%~1"
if "%SRIM%"=="" set "SRIM=%USERPROFILE%\Applications\SRIM-2013-Pro"
set "CASES=%~2"
if "%CASES%"=="" set "CASES=%SRIM%\_amcset\results\trim3"
for %%I in ("%CASES%") do set "NAME=%%~nxI"
set "OUT=%CD%\build\trim_validation\%NAME%"
if not exist "%OUT%" mkdir "%OUT%"
call build\Release\generators\conanrun.bat >NUL 2>&1
echo Comparing %CASES% into %OUT%
set "TABLE_ARG="
if not "%STOPPING_TABLE%"=="" set TABLE_ARG=--stopping-table "%STOPPING_TABLE%"
python3 tools\trim_compare\batch.py "%CASES%" "%OUT%" --amcset build\Release\amcset.exe %TABLE_ARG% ^
  --seeds 101,102,103,104,105,106,107,108 --replicate-seeds 101,102,103,104,105,106,107,108,109,110,111,112 > "%OUT%\batch_log.txt" 2>&1
echo %ERRORLEVEL%> "%OUT%\batch_exit.txt"
REM Statistics against every TRIM run of the same input (also the older reference
REM runs in <SRIM>\_amcset\results\trim, if present).
set "MORE="
if exist "%SRIM%\_amcset\results\trim" if /I not "%CASES%"=="%SRIM%\_amcset\results\trim" set "MORE=%SRIM%\_amcset\results\trim"
python3 tools\trim_compare\analyze.py "%OUT%" "%CASES%" %MORE% -v --json "%OUT%\analysis.json" > "%OUT%\analysis.txt" 2>&1
echo Done. See %OUT%
endlocal
