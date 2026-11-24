# Installs and registers a GitHub Actions self-hosted runner on this machine
# for the "Build and test" workflow (.github/workflows/tests.yml).
#
# 1. On GitHub: repository > Settings > Actions > Runners > New self-hosted
#    runner > Windows. Copy the token from the "Configure" step (it is valid
#    for about an hour).
# 2. powershell -ExecutionPolicy Bypass -File tools\windows\setup_runner.ps1 -Token <token> [-Autostart]
#
# The runner runs as the current user (not as a service), so it sees the
# user's Conan cache (Qt, built once) and Python. -Autostart adds a scheduled
# task that starts it at logon; without it, start C:\actions-runner\run.cmd.
param(
  [Parameter(Mandatory = $true)][string]$Token,
  [string]$Repo = 'https://github.com/Shao-RMSGAL/AMCSET',
  [string]$Dir = 'C:\actions-runner',
  [string]$Name = $env:COMPUTERNAME,
  [switch]$Autostart
)
$ErrorActionPreference = 'Stop'
[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12

New-Item -ItemType Directory -Force $Dir | Out-Null
if (-not (Test-Path (Join-Path $Dir 'config.cmd'))) {
  $rel = Invoke-RestMethod 'https://api.github.com/repos/actions/runner/releases/latest' -Headers @{ 'User-Agent' = 'amcset-setup' }
  $ver = $rel.tag_name.TrimStart('v')
  $asset = "actions-runner-win-x64-$ver.zip"
  $zip = Join-Path $env:TEMP $asset
  Write-Host "Downloading $asset ..."
  Invoke-WebRequest "https://github.com/actions/runner/releases/download/v$ver/$asset" -OutFile $zip -UseBasicParsing
  $m = [regex]::Match($rel.body, "$([regex]::Escape($asset))\s*<!--\s*hash\s*:\s*([0-9a-f]{64})")
  if ($m.Success) {
    $h = (Get-FileHash $zip -Algorithm SHA256).Hash.ToLower()
    if ($h -ne $m.Groups[1].Value) { throw "SHA-256 mismatch for $asset" }
    Write-Host 'SHA-256 verified.'
  } else { Write-Warning 'No published hash found in the release notes; not verified.' }
  Expand-Archive $zip $Dir -Force
  Remove-Item $zip
}
Push-Location $Dir
try {
  & .\config.cmd --unattended --replace --url $Repo --token $Token --name $Name `
    --labels 'windows,x64' --work '_work'
  if ($LASTEXITCODE -ne 0) { throw "config.cmd failed ($LASTEXITCODE)" }
} finally { Pop-Location }

if ($Autostart) {
  $action  = New-ScheduledTaskAction -Execute (Join-Path $Dir 'run.cmd') -WorkingDirectory $Dir
  $trigger = New-ScheduledTaskTrigger -AtLogOn -User $env:USERNAME
  $set     = New-ScheduledTaskSettingsSet -ExecutionTimeLimit ([TimeSpan]::Zero) -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries
  Register-ScheduledTask -TaskName 'GitHub Actions runner (AMCSET)' -Action $action -Trigger $trigger -Settings $set -Force | Out-Null
  Write-Host 'Scheduled task added: the runner starts at logon.'
}
Write-Host "Runner registered. Start it now with: $Dir\run.cmd"
