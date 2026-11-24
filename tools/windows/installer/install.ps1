# Installer payload script: runs inside the self-extracting AMCSET-Setup.exe
# (see tools\windows\make_installer.ps1). Per-user install, no administrator
# rights needed.
$ErrorActionPreference = 'Stop'
$here    = Split-Path -Parent $MyInvocation.MyCommand.Path
$appName = 'AMCSET'
$version = (Get-Content (Join-Path $here 'version.txt') -Raw).Trim()
$target  = Join-Path $env:LOCALAPPDATA "Programs\$appName"
$startDir = Join-Path $env:APPDATA "Microsoft\Windows\Start Menu\Programs\$appName"
Add-Type -AssemblyName System.Windows.Forms

function Show-Message($text, $icon = 'Information') {
  [void][System.Windows.Forms.MessageBox]::Show($text, "$appName setup", 'OK', $icon)
}

try {
  Write-Host "Installing $appName $version to $target ..."
  if (Get-Process amcset-gui -ErrorAction SilentlyContinue |
      Where-Object { $_.Path -like "$target\*" }) {
    Show-Message "Please close AMCSET and run the installer again." 'Warning'
    exit 1
  }
  if (Test-Path $target) { Remove-Item $target -Recurse -Force }
  New-Item -ItemType Directory -Force $target | Out-Null
  Expand-Archive -Path (Join-Path $here 'payload.zip') -DestinationPath $target -Force

  # Shortcuts
  $shell = New-Object -ComObject WScript.Shell
  New-Item -ItemType Directory -Force $startDir | Out-Null
  $lnk = $shell.CreateShortcut((Join-Path $startDir "$appName.lnk"))
  $lnk.TargetPath = Join-Path $target 'amcset-gui.exe'
  $lnk.WorkingDirectory = [Environment]::GetFolderPath('MyDocuments')
  $lnk.Description = 'AMCSET ion transport and radiation damage simulator'
  $lnk.Save()
  $lnk = $shell.CreateShortcut((Join-Path $startDir "$appName command prompt.lnk"))
  $lnk.TargetPath = "$env:windir\System32\cmd.exe"
  $lnk.Arguments = "/k set `"PATH=$target;%PATH%`" && cd /d `"%USERPROFILE%`" && amcset --help"
  $lnk.Description = 'Command prompt with the amcset command on the PATH'
  $lnk.Save()
  $lnk = $shell.CreateShortcut((Join-Path $startDir "Uninstall $appName.lnk"))
  $lnk.TargetPath = "$env:windir\System32\WindowsPowerShell\v1.0\powershell.exe"
  $lnk.Arguments = "-NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File `"$target\uninstall.ps1`""
  $lnk.IconLocation = (Join-Path $target 'amcset-gui.exe')
  $lnk.Description = "Remove $appName from this computer"
  $lnk.Save()
  $lnk = $shell.CreateShortcut((Join-Path ([Environment]::GetFolderPath('Desktop')) "$appName.lnk"))
  $lnk.TargetPath = Join-Path $target 'amcset-gui.exe'
  $lnk.WorkingDirectory = [Environment]::GetFolderPath('MyDocuments')
  $lnk.Save()

  # "Apps > Installed apps" entry (current user)
  $key = "HKCU:\Software\Microsoft\Windows\CurrentVersion\Uninstall\$appName"
  New-Item -Path $key -Force | Out-Null
  $props = @{
    DisplayName     = "$appName $version"
    DisplayVersion  = $version
    Publisher       = 'Texas A&M University'
    InstallLocation = $target
    DisplayIcon     = (Join-Path $target 'amcset-gui.exe')
    UninstallString = "powershell.exe -NoProfile -ExecutionPolicy Bypass -WindowStyle Hidden -File `"$target\uninstall.ps1`""
  }
  foreach ($k in $props.Keys) { New-ItemProperty -Path $key -Name $k -Value $props[$k] -Force | Out-Null }
  New-ItemProperty -Path $key -Name NoModify -Value 1 -PropertyType DWord -Force | Out-Null
  New-ItemProperty -Path $key -Name NoRepair -Value 1 -PropertyType DWord -Force | Out-Null

  $answer = [System.Windows.Forms.MessageBox]::Show(
    "$appName $version is installed in`n$target`n`nStart it now?", "$appName setup", 'YesNo', 'Information')
  if ($answer -eq 'Yes') { Start-Process (Join-Path $target 'amcset-gui.exe') }
} catch {
  Show-Message ("Installation failed:`n" + $_.Exception.Message) 'Error'
  exit 1
}
