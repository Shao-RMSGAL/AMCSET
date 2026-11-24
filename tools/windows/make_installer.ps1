# Builds build\installer\AMCSET-Setup.exe: a single self-extracting installer
# (Windows IExpress, no extra tools) from the Release build in build\Release.
# Run tools\windows\build_and_test.bat first. Usage:
#   powershell -ExecutionPolicy Bypass -File tools\windows\make_installer.ps1
$ErrorActionPreference = 'Stop'
$root    = Resolve-Path (Join-Path $PSScriptRoot '..\..')
$rel     = Join-Path $root 'build\Release'
$out     = Join-Path $root 'build\installer'
$stage   = Join-Path $out 'stage'
$app     = Join-Path $stage 'app'
$pack    = Join-Path $out 'pack'
$version = ([regex]::Match((Get-Content (Join-Path $root 'CMakeLists.txt') -Raw),
            'VERSION\s+(\d+\.\d+\.\d+)')).Groups[1].Value
if (-not $version) { $version = '0.1.0' }

foreach ($d in $stage, $pack) { if (Test-Path $d) { Remove-Item $d -Recurse -Force } }
New-Item -ItemType Directory -Force $app, $pack | Out-Null

# Application, the Qt libraries it links, and the Qt plugins it needs.
foreach ($f in 'amcset.exe', 'amcset-gui.exe', 'Qt6Core.dll', 'Qt6Gui.dll', 'Qt6Widgets.dll', 'd3dcompiler_47.dll') {
  Copy-Item (Join-Path $rel $f) $app
}
foreach ($p in 'platforms\qwindows.dll', 'styles\qmodernwindowsstyle.dll',
               'imageformats\qgif.dll', 'imageformats\qico.dll') {
  $dst = Join-Path $app (Split-Path $p -Parent)
  New-Item -ItemType Directory -Force $dst | Out-Null
  Copy-Item (Join-Path $rel $p) $dst
}
# App-local Visual C++ runtime, so no redistributable (and no admin) is needed.
$crt = Get-ChildItem (Join-Path ${env:ProgramFiles} 'Microsoft Visual Studio\*\*\VC\Redist\MSVC\*\x64\Microsoft.VC*.CRT') -Directory |
       Sort-Object FullName | Select-Object -Last 1
foreach ($f in 'vcruntime140.dll', 'vcruntime140_1.dll', 'msvcp140.dll', 'msvcp140_1.dll', 'msvcp140_2.dll') {
  Copy-Item (Join-Path $crt.FullName $f) $app
}
# SRIM-2013 stopping tables (found at run time in data\ next to the exe).
Copy-Item (Join-Path $root 'data') (Join-Path $app 'data') -Recurse
foreach ($f in 'README.md', 'LICENSE', 'CHANGELOG.md') {
  if (Test-Path (Join-Path $root $f)) { Copy-Item (Join-Path $root $f) $app }
}
$docs = Join-Path $app 'docs'; New-Item -ItemType Directory -Force $docs | Out-Null
foreach ($f in 'CLI.md', 'GUI.md', 'SRIM_MODE.md', 'OUTPUT.md', 'GLOSSARY.md') {
  Copy-Item (Join-Path $root "docs\$f") $docs
}
Copy-Item (Join-Path $PSScriptRoot 'installer\uninstall.ps1') $app

# payload.zip (forward-slash entries) + installer script.
Add-Type -AssemblyName System.IO.Compression.FileSystem
$zip = Join-Path $pack 'payload.zip'
[System.IO.Compression.ZipFile]::CreateFromDirectory($app, $zip,
  [System.IO.Compression.CompressionLevel]::Optimal, $false)
Copy-Item (Join-Path $PSScriptRoot 'installer\install.ps1') $pack
Set-Content (Join-Path $pack 'version.txt') $version -Encoding ASCII

# IExpress directive file.
$exe = Join-Path $out 'AMCSET-Setup.exe'
if (Test-Path $exe) { Remove-Item $exe -Force }
$sed = @"
[Version]
Class=IEXPRESS
SEDVersion=3
[Options]
PackagePurpose=InstallApp
ShowInstallProgramWindow=0
HideExtractAnimation=1
UseLongFileName=1
InsideCompressed=0
CAB_FixedSize=0
CAB_ResvCodeSigning=0
RebootMode=N
InstallPrompt=%InstallPrompt%
DisplayLicense=
FinishMessage=
TargetName=%TargetName%
FriendlyName=%FriendlyName%
AppLaunched=%AppLaunched%
PostInstallCmd=<None>
AdminQuietInstCmd=%AppLaunched%
UserQuietInstCmd=%AppLaunched%
SourceFiles=SourceFiles
[Strings]
InstallPrompt=Install AMCSET $version for the current user?
TargetName=$exe
FriendlyName=AMCSET $version setup
AppLaunched=powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File install.ps1
FILE0="install.ps1"
FILE1="payload.zip"
FILE2="version.txt"
[SourceFiles]
SourceFiles0=$pack\
[SourceFiles0]
%FILE0%=
%FILE1%=
%FILE2%=
"@
$sedPath = Join-Path $out 'AMCSET.sed'
Set-Content $sedPath $sed -Encoding ASCII
$p = Start-Process "$env:windir\System32\iexpress.exe" -ArgumentList '/N', '/Q', $sedPath -Wait -PassThru -NoNewWindow
if (-not (Test-Path $exe)) { throw "IExpress did not create $exe (exit code $($p.ExitCode))" }
"{0}  {1:N1} MB" -f $exe, ((Get-Item $exe).Length / 1MB)
