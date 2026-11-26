# Signs AMCSET's Windows files (exe, dll, msi, msix) with signtool, adds an
# RFC 3161 time stamp, and verifies each signature. Files that already carry a
# valid signature (e.g. Microsoft's C++ runtime) are left alone. See
# docs/RELEASING.md "Code signing" for setting up a certificate.
#
# Azure Artifact Signing (Microsoft's signing service):
#   powershell -File tools\windows\sign.ps1 -Metadata metadata.json `
#       -Dlib <...>\bin\x64\Azure.CodeSigning.Dlib.dll  <files or folders>
#   metadata.json: {"Endpoint": "https://eus.codesigning.azure.net",
#                   "CodeSigningAccountName": "...", "CertificateProfileName": "..."}
#   (signs as the identity of `az login` or the environment's Azure credentials;
#   needs the .NET 8 runtime and signtool 10.0.22621.755 or later.)
# A certificate in the Windows certificate store (OV/EV, also on a token):
#   powershell -File tools\windows\sign.ps1 -Thumbprint <SHA-1 thumbprint> <files or folders>
# A .pfx file (password in the environment variable AMCSET_PFX_PASSWORD):
#   powershell -File tools\windows\sign.ps1 -Pfx cert.pfx <files or folders>
# Folders are searched recursively for *.exe, *.dll, *.msi, *.msix.
param(
  [Parameter(Position = 0, ValueFromRemainingArguments = $true)] [string[]] $Paths,
  [string] $Metadata,
  [string] $Dlib,
  [string] $Thumbprint,
  [string] $Pfx,
  [string] $TimestampUrl,
  [string] $SignTool,
  [switch] $Force
)
$ErrorActionPreference = 'Stop'

if (-not $Paths) { throw 'Name the files or folders to sign.' }
$modes = @([bool]$Metadata, [bool]$Thumbprint, [bool]$Pfx) | Where-Object { $_ }
if ($modes.Count -ne 1) { throw 'Give exactly one of -Metadata (with -Dlib), -Thumbprint or -Pfx.' }

if (-not $SignTool) {
  # The newest Windows SDK's x64 signtool.
  $SignTool = Get-ChildItem "${env:ProgramFiles(x86)}\Windows Kits\10\bin\10.*\x64\signtool.exe" -ErrorAction SilentlyContinue |
              Sort-Object { [version]($_.Directory.Parent.Name) } | Select-Object -Last 1 -ExpandProperty FullName
  if (-not $SignTool) { throw 'signtool.exe not found: install the Windows SDK or give -SignTool.' }
}

$args_sign = @('sign', '/fd', 'SHA256', '/td', 'SHA256')
if ($Metadata) {
  if (-not $Dlib -or -not (Test-Path $Dlib)) { throw '-Metadata needs -Dlib <path to Azure.CodeSigning.Dlib.dll> (NuGet package Microsoft.ArtifactSigning.Client).' }
  if (-not $TimestampUrl) { $TimestampUrl = 'http://timestamp.acs.microsoft.com' }
  $args_sign += @('/tr', $TimestampUrl, '/dlib', (Resolve-Path $Dlib).Path, '/dmdf', (Resolve-Path $Metadata).Path)
} elseif ($Thumbprint) {
  if (-not $TimestampUrl) { throw '-Thumbprint needs -TimestampUrl: the RFC 3161 time-stamp server of the certificate authority.' }
  $args_sign += @('/tr', $TimestampUrl, '/sha1', $Thumbprint)
} else {
  if (-not $TimestampUrl) { throw '-Pfx needs -TimestampUrl: the RFC 3161 time-stamp server of the certificate authority.' }
  $args_sign += @('/tr', $TimestampUrl, '/f', (Resolve-Path $Pfx).Path)
  if ($env:AMCSET_PFX_PASSWORD) { $args_sign += @('/p', $env:AMCSET_PFX_PASSWORD) }
}

$files = foreach ($p in $Paths) {
  if (Test-Path $p -PathType Container) {
    Get-ChildItem $p -Recurse -File -Include *.exe, *.dll, *.msi, *.msix | Select-Object -ExpandProperty FullName
  } elseif (Test-Path $p -PathType Leaf) {
    (Resolve-Path $p).Path
  } else { throw "Not found: $p" }
}
$todo = @()
foreach ($f in $files) {
  $s = Get-AuthenticodeSignature $f
  if ($s.Status -eq 'Valid' -and -not $Force) { Write-Host "already signed: $f ($($s.SignerCertificate.Subject))"; continue }
  $todo += $f
}
if (-not $todo) { Write-Host 'Nothing to sign.'; exit 0 }

# signtool takes many files per call; keep the command lines short.
for ($i = 0; $i -lt $todo.Count; $i += 20) {
  $chunk = $todo[$i..([Math]::Min($i + 19, $todo.Count - 1))]
  & $SignTool @args_sign @chunk
  if ($LASTEXITCODE -ne 0) { throw "signtool sign failed ($LASTEXITCODE)" }
}
foreach ($f in $todo) {
  & $SignTool verify /pa /q $f
  if ($LASTEXITCODE -ne 0) { throw "signature does not verify: $f" }
}
Write-Host "Signed and verified $($todo.Count) file(s)."
