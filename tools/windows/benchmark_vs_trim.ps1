# Speed of AMCSET's SRIM mode against SRIM-2013's TRIM.exe on the same inputs.
# Run through benchmark_vs_trim.bat (which puts the runtime libraries on PATH).
#
# For each case the TRIM.IN of a TRIM reference run is copied with the ion
# count below, then run by TRIM.exe alone (in a private copy of the SRIM
# folder under build\benchmark\trim_work, TRIMAUTO = 1) and by amcset.exe with
# one thread and with all threads. Nothing else should run meanwhile.
# Results: build\benchmark\benchmark.csv (wall-clock seconds).
param([string]$Srim = (Join-Path $env:USERPROFILE 'Applications\SRIM-2013-Pro'),
      [string]$Only = '')
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..\..')).Path
$out = Join-Path $root 'build\benchmark'
$work = Join-Path $out 'trim_work'
$amcset = Join-Path $root 'build\Release\amcset.exe'
New-Item -ItemType Directory -Force -Path $out | Out-Null
if (-not (Test-Path (Join-Path $work 'TRIM.exe'))) {
  New-Item -ItemType Directory -Force -Path $work | Out-Null
  Get-ChildItem $Srim -File | Copy-Item -Destination $work
  Copy-Item -Recurse -Force (Join-Path $Srim 'Data') $work
}
foreach ($d in 'SRIM Outputs', 'SRIM Restore') { New-Item -ItemType Directory -Force -Path (Join-Path $work $d) | Out-Null }

# name, ions (chosen so that TRIM needs one to a few minutes)
$cases = @(
  @('008_fe_fe_1000keV_kp', 2000), @('005_fe_fe_1000keV_fc', 200),
  @('011_fe_fe_100keV_fc', 500), @('013_fe_fe_10keV_fc', 2000),
  @('015_fe_fe_1keV_fc', 5000), @('044_ar_cu_1keV_fc', 5000),
  @('027_h_fe_1000keV_kp', 5000), @('036_he_fe_1000keV_kp', 2000),
  @('035_he_fe_1000keV_fc', 1000), @('052_fe_316ss_2000keV_fc', 100),
  @('053_fe_316ss_2000keV_kp', 500), @('046_co_co_1000keV_fc', 200))
$threads = [Environment]::ProcessorCount
$csv = Join-Path $out 'benchmark.csv'
if (-not (Test-Path $csv)) {
  Set-Content -Path $csv -Value 'case,ions,trim_s,amcset_1_thread_s,amcset_all_threads_s,threads,speedup_1_thread,speedup_all_threads'
}
$done = (Get-Content $csv | Select-Object -Skip 1 | ForEach-Object { ($_ -split ',')[0] })
foreach ($c in $cases) {
  $name = $c[0]; $ions = $c[1]
  if ($Only -and $name -notmatch $Only) { continue }
  if ($done -contains $name) { continue }
  $src = $null
  foreach ($d in 'trim', 'trim3') {
    $p = Join-Path $Srim "_amcset\results\$d\$name\TRIM.IN"
    if (Test-Path $p) { $src = $p; break }
  }
  if (-not $src) { Write-Host "skip $name (no TRIM.IN)"; continue }
  $lines = Get-Content $src
  $t = ($lines[2].Trim() -split '\s+')
  $t[4] = "$ions"; $t[6] = "$ions"
  $lines[2] = '     ' + ($t -join '   ')
  $case_dir = Join-Path $out "inputs\$name"
  New-Item -ItemType Directory -Force -Path $case_dir | Out-Null
  Set-Content -Path (Join-Path $case_dir 'TRIM.IN') -Value $lines -Encoding ASCII
  Set-Content -Path (Join-Path $work 'TRIM.IN') -Value $lines -Encoding ASCII
  Set-Content -Path (Join-Path $work 'TRIMAUTO') -Value @('1', '') -Encoding ASCII
  Get-ChildItem (Join-Path $work 'SRIM Outputs') -File | Remove-Item -Force
  Write-Host "$name ($ions ions): TRIM..."
  $sw = [Diagnostics.Stopwatch]::StartNew()
  Start-Process -FilePath (Join-Path $work 'TRIM.exe') -WorkingDirectory $work -Wait
  $trim_s = $sw.Elapsed.TotalSeconds
  $times = @()
  foreach ($j in 1, 0) {
    $o = Join-Path $case_dir ("amcset_j$j")
    $sw = [Diagnostics.Stopwatch]::StartNew()
    $args_ = @('trim', ('"' + (Join-Path $case_dir 'TRIM.IN') + '"'), '-q', '-s', '1', '-j', "$j",
               '--formats', 'trim', '-o', ('"' + $o + '"'))
    $p = Start-Process -FilePath $amcset -ArgumentList $args_ -NoNewWindow -Wait -PassThru `
           -RedirectStandardOutput (Join-Path $case_dir "amcset_j$j.log") -RedirectStandardError (Join-Path $case_dir "amcset_j$j.err")
    if ($p.ExitCode -ne 0) { throw "amcset failed on $name" }
    $times += $sw.Elapsed.TotalSeconds
  }
  $line = '{0},{1},{2:F1},{3:F2},{4:F2},{5},{6:F1},{7:F1}' -f $name, $ions, $trim_s, $times[0], $times[1], $threads,
          ($trim_s / $times[0]), ($trim_s / $times[1])
  Add-Content -Path $csv -Value $line
  Write-Host $line
}
Write-Host "Done: $csv"
