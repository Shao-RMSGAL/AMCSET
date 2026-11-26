# Run the TRIM reference cases of batch 3 (folder trim_cases3, results in results\trim3). Based on version 2.
# Version 2 also collects the summary files TRIM writes next to TRIM.exe in batch mode
# (RANGE.txt, VACANCY.txt, NOVAC.txt, IONIZ.txt, PHONON.txt, E2RECOIL.txt, LATERAL.txt, TDATA.txt).
# A case is complete when results\trim\<case>\RANGE.txt and _done.txt exist. Cases finished by
# version 1 (run_trim_cases.ps1) without those files are salvaged from the version-1 worker folder
# when possible, otherwise re-run. At most -MaxProcs TRIM.exe processes run at once (all versions).
# -CaseSet / -OutSet choose another batch (e.g. trim_cases4 -> results\trim4, gen_cases4.py).
param([int]$Workers = 0, [int]$MaxProcs = 0, [int]$TimeoutHours = 12, [string]$Filter = '*',
      [string]$CaseSet = 'trim_cases3', [string]$OutSet = 'trim3')
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$srim = Split-Path -Parent $root
if ($Workers -le 0) { $Workers = [Math]::Max(1, [Environment]::ProcessorCount - 2) }
if ($MaxProcs -le 0) { $MaxProcs = $Workers }
$work   = Join-Path $root 'work_trim2'
$work1  = Join-Path $root 'work_trim'
$out    = Join-Path $root ('results\' + $OutSet)
$claim  = Join-Path $work $(if ($OutSet -eq 'trim3') { 'claims3' } else { 'claims_' + $OutSet })
New-Item -ItemType Directory -Force -Path $work, $out, $claim | Out-Null
$cases = Get-ChildItem (Join-Path $root $CaseSet) -Directory -Filter $Filter | Sort-Object Name
Write-Host ("TRIM reference cases: {0}, {1} workers, at most {2} TRIM processes" -f $cases.Count, $Workers, $MaxProcs)
$skip = @('_amcset', 'SRIM Tutorials', 'SRIMHelp', 'SRIMPICS', 'Demo', 'SRIM-Setup', 'SR Module')
for ($w = 0; $w -lt $Workers; $w++) {
  $d = Join-Path $work ("w" + $w)
  if (-not (Test-Path (Join-Path $d 'TRIM.exe'))) {
    New-Item -ItemType Directory -Force -Path $d | Out-Null
    Get-ChildItem $srim | Where-Object { $skip -notcontains $_.Name } | ForEach-Object { Copy-Item -Recurse -Force $_.FullName $d }
    New-Item -ItemType Directory -Force -Path (Join-Path $d 'SRIM Outputs'), (Join-Path $d 'SRIM Restore') | Out-Null
  }
}
$worker = {
  param($wdir, $casedirs, $out, $claim, $claim1, $work1, $timeoutMs, $maxProcs)
  $rootFiles = @('RANGE.txt', 'VACANCY.txt', 'NOVAC.txt', 'IONIZ.txt', 'PHONON.txt', 'E2RECOIL.txt', 'LATERAL.txt', 'TDATA.txt')
  function Test-Complete($dest) { (Test-Path (Join-Path $dest '_done.txt')) -and (Test-Path (Join-Path $dest 'RANGE.txt')) }
  function Try-Salvage($c, $dest) {
    $done = Join-Path $dest '_done.txt'
    if (-not (Test-Path $done) -or (Test-Path (Join-Path $dest 'RANGE.txt'))) { return $false }
    $m = [regex]::Match((Get-Content $done -Raw), 'worker=(w\d+)')
    if (-not $m.Success) { return $false }
    $w1 = Join-Path $work1 $m.Groups[1].Value
    $tin = Join-Path $w1 'TRIM.IN'
    if (-not (Test-Path $tin) -or -not (Test-Path (Join-Path $w1 'RANGE.txt'))) { return $false }
    if ((Get-FileHash $tin).Hash -ne (Get-FileHash (Join-Path $c 'TRIM.IN')).Hash) { return $false }
    if ((Get-Item (Join-Path $w1 'RANGE.txt')).LastWriteTime -lt (Get-Item $done).LastWriteTime.AddMinutes(-10)) { return $false }
    foreach ($f in $rootFiles) { $p = Join-Path $w1 $f; if (Test-Path $p) { Copy-Item -Force $p $dest } }
    Add-Content -Path $done -Value 'salvaged_root_files=True'
    return $true
  }
  while ($true) {
    $pending = 0; $ran = $false
    foreach ($c in $casedirs) {
      $name = Split-Path -Leaf $c
      $dest = Join-Path $out $name
      if (Test-Complete $dest) { continue }
      if (Try-Salvage $c $dest) { continue }
      $pending++
      # still being run by the version-1 runner?
      if ((Test-Path (Join-Path $claim1 ($name + '.claim'))) -and -not (Test-Path (Join-Path $dest '_done.txt'))) { continue }
      if (Test-Path (Join-Path $claim ($name + '.claim'))) { continue }
      while (@(Get-Process -Name TRIM -ErrorAction SilentlyContinue).Count -ge $maxProcs) { Start-Sleep -Seconds 15 }
      try { New-Item -ItemType File -Path (Join-Path $claim ($name + '.claim')) -ErrorAction Stop | Out-Null } catch { continue }
      foreach ($sub in 'SRIM Outputs', 'SRIM Restore') { Get-ChildItem (Join-Path $wdir $sub) -File -ErrorAction SilentlyContinue | Remove-Item -Force }
      foreach ($f in @('TRIM.DAT', 'TRIM.dat') + $rootFiles) { $p = Join-Path $wdir $f; if (Test-Path $p) { Remove-Item -Force $p } }
      Copy-Item -Force (Join-Path $c '*') $wdir
      Set-Content -Path (Join-Path $wdir 'TRIMAUTO') -Value @('1', '') -Encoding ASCII
      $t0 = Get-Date
      $p = Start-Process -FilePath (Join-Path $wdir 'TRIM.exe') -WorkingDirectory $wdir -PassThru
      $finished = $p.WaitForExit($timeoutMs)
      if (-not $finished) { try { $p.Kill() } catch {} }
      $secs = ((Get-Date) - $t0).TotalSeconds
      New-Item -ItemType Directory -Force -Path $dest | Out-Null
      Get-ChildItem $dest -File | Remove-Item -Force
      foreach ($sub in 'SRIM Outputs', 'SRIM Restore') { Get-ChildItem (Join-Path $wdir $sub) -File -ErrorAction SilentlyContinue | Copy-Item -Destination $dest -Force }
      foreach ($f in $rootFiles) { $q = Join-Path $wdir $f; if (Test-Path $q) { Copy-Item -Force $q $dest } }
      Copy-Item -Force (Join-Path $wdir 'TRIM.IN') $dest
      if (-not (Test-Path (Join-Path $dest 'RANGE.txt'))) { Set-Content -Path (Join-Path $dest 'RANGE.txt') -Value 'missing: TRIM wrote no RANGE.txt' }
      Set-Content -Path (Join-Path $dest '_done.txt') -Value ("finished={0} seconds={1:F1} worker={2} end={3} runner=3" -f $finished, $secs, (Split-Path -Leaf $wdir), (Get-Date -Format o))
      $ran = $true
    }
    if ($pending -eq 0) { break }
    if (-not $ran) { Start-Sleep -Seconds 30 }
  }
}
$dirs = $cases | ForEach-Object { $_.FullName }
$running = @()
for ($w = 0; $w -lt $Workers; $w++) {
  $running += Start-Job -ScriptBlock $worker -ArgumentList (Join-Path $work ("w" + $w)), $dirs, $out, $claim, (Join-Path $work1 'claims'), $work1, ($TimeoutHours * 3600000), $MaxProcs
  Start-Sleep -Seconds 3
}
while (($running | Where-Object { $_.State -eq 'Running' }).Count -gt 0) {
  $n = (Get-ChildItem $out -Directory | Where-Object { (Test-Path (Join-Path $_.FullName '_done.txt')) -and (Test-Path (Join-Path $_.FullName 'RANGE.txt')) }).Count
  $np = @(Get-Process -Name TRIM -ErrorAction SilentlyContinue).Count
  Write-Host ("[{0}] {1}/{2} cases complete, {3} TRIM processes" -f (Get-Date -Format T), $n, $cases.Count, $np)
  Set-Content -Path (Join-Path $out '_progress3.txt') -Value ("{0} {1}/{2} trim_processes={3}" -f (Get-Date -Format o), $n, $cases.Count, $np)
  Start-Sleep -Seconds 30
}
$running | Receive-Job | Out-File (Join-Path $out '_runner3_log.txt')
Set-Content -Path (Join-Path $out '_all_done3.txt') -Value (Get-Date -Format o)
Write-Host "Done."
