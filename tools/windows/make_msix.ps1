# Builds build\packages\AMCSET-<version>-windows-x64.msix from a Release build
# (tools\windows\build_and_test.bat, or the release workflow): installs the
# build into a staging folder, adds the logos, the package manifest
# (packaging\windows\AppxManifest.xml.in) and the resource index, and packs it
# with the Windows SDK's makeappx (which validates the manifest).
#
#   powershell -ExecutionPolicy Bypass -File tools\windows\make_msix.ps1 `
#       [-BuildDir build\Release] [-Publisher "CN=..."] [-Out build\packages]
#
# An MSIX installs only when it is signed by a certificate the computer trusts,
# and -Publisher must then be exactly the subject of that certificate (for
# Azure Artifact Signing: the subject shown on the certificate profile). Sign
# the result with tools\windows\sign.ps1. Without a certificate the package
# can still be checked (makeappx unpack, the App Certification Kit).
param(
  [string] $BuildDir = 'build\Release',
  [string] $Publisher = 'CN=AMCSET Development',
  [string] $Out = 'build\packages',
  [string] $CMake = 'cmake'
)
$ErrorActionPreference = 'Stop'
$root = Resolve-Path (Join-Path $PSScriptRoot '..\..')
$version = ([regex]::Match((Get-Content (Join-Path $root 'CMakeLists.txt') -Raw), 'VERSION\s+(\d+\.\d+\.\d+)')).Groups[1].Value
if (-not $version) { throw 'No version in CMakeLists.txt' }

$sdk = Get-ChildItem "${env:ProgramFiles(x86)}\Windows Kits\10\bin\10.*\x64\makeappx.exe" -ErrorAction SilentlyContinue |
       Sort-Object { [version]($_.Directory.Parent.Name) } | Select-Object -Last 1
if (-not $sdk) { throw 'makeappx.exe not found: install the Windows SDK.' }
$makeappx = $sdk.FullName
$makepri = Join-Path $sdk.DirectoryName 'makepri.exe'

$stage = Join-Path $root 'build\msix\stage'
if (Test-Path $stage) { Remove-Item $stage -Recurse -Force }
New-Item -ItemType Directory -Force $stage | Out-Null

# The installed files (bin\, share\) under AMCSET\.
& $CMake --install (Join-Path $root $BuildDir) --prefix (Join-Path $stage 'AMCSET')
if ($LASTEXITCODE -ne 0) { throw "cmake --install failed ($LASTEXITCODE)" }

# Logos: the scale and target-size variants (resolved through resources.pri)
# and plain copies (used if no index is built).
$assets = Join-Path $stage 'Assets'
New-Item -ItemType Directory -Force $assets | Out-Null
Copy-Item (Join-Path $root 'resources\icons\msix\*.png') $assets
foreach ($n in 'Square44x44Logo', 'Square150x150Logo', 'StoreLogo') {
  Copy-Item (Join-Path $assets "$n.scale-100.png") (Join-Path $assets "$n.png")
}

$manifest = (Get-Content (Join-Path $root 'packaging\windows\AppxManifest.xml.in') -Raw)
$manifest = $manifest.Replace('@VERSION@', $version).Replace('@PUBLISHER@', [Security.SecurityElement]::Escape($Publisher))
[IO.File]::WriteAllText((Join-Path $stage 'AppxManifest.xml'), $manifest, (New-Object Text.UTF8Encoding $false))

# The resource index of the logo variants.
if (Test-Path $makepri) {
  $cfg = Join-Path $root 'build\msix\priconfig.xml'
  & $makepri createconfig /cf $cfg /dq en-US /pv 10.0.0 /o | Out-Null
  # One index in the package (no resource packs, which only bundles use).
  [xml]$x = Get-Content $cfg
  $pk = $x.SelectSingleNode('//packaging')
  if ($pk) { [void]$pk.ParentNode.RemoveChild($pk); $x.Save($cfg) }
  & $makepri new /pr $stage /cf $cfg /mn (Join-Path $stage 'AppxManifest.xml') /of (Join-Path $stage 'resources.pri') /o | Out-Null
  if ($LASTEXITCODE -ne 0) { throw "makepri failed ($LASTEXITCODE)" }
}

New-Item -ItemType Directory -Force (Join-Path $root $Out) | Out-Null
$msix = Join-Path (Join-Path $root $Out) "AMCSET-$version-windows-x64.msix"
& $makeappx pack /o /h SHA256 /d $stage /p $msix
if ($LASTEXITCODE -ne 0) { throw "makeappx failed ($LASTEXITCODE)" }
"{0}  {1:N1} MB" -f $msix, ((Get-Item $msix).Length / 1MB)
