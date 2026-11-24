# Removes AMCSET (installed by AMCSET-Setup.exe). Lives in the install folder;
# started from Settings > Installed apps or the Start Menu entry.
# -Quiet skips the confirmation and the final message.
param([switch]$Quiet)
Add-Type -AssemblyName System.Windows.Forms
$app = Split-Path -Parent $MyInvocation.MyCommand.Path
if (-not $Quiet) {
  $a = [System.Windows.Forms.MessageBox]::Show(
    "Remove AMCSET from`n$app ?`n`nYour own scenario and result files are not touched.",
    'AMCSET uninstall', 'YesNo', 'Question')
  if ($a -ne 'Yes') { exit 0 }
}
Get-Process amcset-gui, amcset -ErrorAction SilentlyContinue |
  Where-Object { $_.Path -like "$app\*" } | Stop-Process -Force
Remove-Item (Join-Path $env:APPDATA 'Microsoft\Windows\Start Menu\Programs\AMCSET') -Recurse -Force -ErrorAction SilentlyContinue
Remove-Item (Join-Path ([Environment]::GetFolderPath('Desktop')) 'AMCSET.lnk') -Force -ErrorAction SilentlyContinue
Remove-Item 'HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\AMCSET' -Recurse -Force -ErrorAction SilentlyContinue
# This script lives in the folder it deletes: a helper cmd removes the folder
# after this PowerShell process has exited.
Start-Process cmd.exe -ArgumentList "/c ping -n 3 127.0.0.1 >nul & rmdir /s /q `"$app`"" -WindowStyle Hidden
if (-not $Quiet) {
  [void][System.Windows.Forms.MessageBox]::Show('AMCSET has been removed.', 'AMCSET uninstall', 'OK', 'Information')
}
