# Extract SRIM-2013 electronic/nuclear stopping tables with SR Module (SRModule.exe).
# Written for AMCSET's "SRIM mode". Reads sr_jobs.csv (one ion/target pair per row) and
# sr_energies_keV.txt (energy list), runs SR Module in several private copies of the
# "SR Module" folder in parallel, and stores each output table in results\stopping\<id>.txt.
# Re-running skips tables that already exist. Close the window to abort.
param([int]$Workers = 0)
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$srim = Split-Path -Parent $root
if ($Workers -le 0) { $Workers = [Math]::Max(1, [Environment]::ProcessorCount - 2) }
$work  = Join-Path $root 'work_sr'
$out   = Join-Path $root 'results\stopping'
$claim = Join-Path $root 'work_sr\claims'
New-Item -ItemType Directory -Force -Path $work, $out, $claim | Out-Null
$jobs = Import-Csv (Join-Path $root 'sr_jobs.csv')
$energies = Get-Content (Join-Path $root 'sr_energies_keV.txt')
Write-Host ("SR Module extraction: {0} tables, {1} workers" -f $jobs.Count, $Workers)
for ($w = 0; $w -lt $Workers; $w++) {
  $d = Join-Path $work ("w" + $w)
  if (-not (Test-Path (Join-Path $d 'SRModule.exe'))) { Copy-Item -Recurse -Force (Join-Path $srim 'SR Module') $d }
}
$worker = {
  param($wdir, $jobs, $energies, $out, $claim)
  foreach ($j in $jobs) {
    $dest = Join-Path $out ($j.id + '.txt')
    if (Test-Path $dest) { continue }
    try { New-Item -ItemType File -Path (Join-Path $claim ($j.id + '.claim')) -ErrorAction Stop | Out-Null } catch { continue }
    $lines = @('---Stopping/Range Input Data (Number-format: Period = Decimal Point)',
               '---Output File Name', '"sr_out.txt"',
               '---Ion(Z), Ion Mass(u)', ("{0}   {1}" -f $j.z1, $j.m1),
               '---Target Data: (Solid=0,Gas=1), Density(g/cm3), Compound Corr.', ("{0}    {1}    1" -f $j.phase, $j.density),
               '---Number of Target Elements', ' 1 ',
               '---Target Elements: (Z), Target name, Stoich, Target Mass(u)', ("{0}   `"{1}`"   1   {2}" -f $j.z2, $j.name, $j.m2),
               '---Output Stopping Units (1-8)', ' 7',
               '---Ion Energy : E-Min(keV), E-Max(keV)', ' 0 0') + $energies + @(' 0')
    Set-Content -Path (Join-Path $wdir 'SR.IN') -Value $lines -Encoding ASCII
    $f = Join-Path $wdir 'sr_out.txt'
    if (Test-Path $f) { Remove-Item -Force $f }
    $p = Start-Process -FilePath (Join-Path $wdir 'SRModule.exe') -WorkingDirectory $wdir -PassThru
    if (-not $p.WaitForExit(180000)) { try { $p.Kill() } catch {} ; Add-Content (Join-Path $out ('_error_' + $j.id + '.txt')) 'timeout'; continue }
    if (Test-Path $f) { Move-Item -Force $f $dest } else { Add-Content (Join-Path $out ('_error_' + $j.id + '.txt')) 'no output' }
  }
}
$running = @()
for ($w = 0; $w -lt $Workers; $w++) {
  $running += Start-Job -ScriptBlock $worker -ArgumentList (Join-Path $work ("w" + $w)), $jobs, $energies, $out, $claim
}
while (($running | Where-Object { $_.State -eq 'Running' }).Count -gt 0) {
  $n = (Get-ChildItem $out -Filter 'sr_*.txt').Count
  Write-Host ("[{0}] {1}/{2} tables" -f (Get-Date -Format T), $n, $jobs.Count)
  Set-Content -Path (Join-Path $out '_progress.txt') -Value ("{0} {1}/{2}" -f (Get-Date -Format o), $n, $jobs.Count)
  Start-Sleep -Seconds 20
}
$running | Receive-Job | Out-Null
Set-Content -Path (Join-Path $out '_done.txt') -Value (Get-Date -Format o)
Write-Host "Done."
