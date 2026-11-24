# Runs the TRIM.DAT recoil cases (cascade codes 6/7) that the batch runner held back:
# TRIM shows a message box for these codes even in batch mode, which must be answered.
# Each case runs in its own copy of the installation (work_trim2\h<k>); outputs go to
# results\trim3\<case> like the batch runner's.
param([string[]]$Cases = @('375_dat_fe_fe_pka_fc', '397_dat_fecr_pka_fc', '418_dat_fe_fe_pka_kp'))
$ErrorActionPreference = 'Stop'
$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$work = Join-Path $root 'work_trim2'
$out  = Join-Path $root 'results\trim3'
$rootFiles = @('RANGE.txt', 'VACANCY.txt', 'NOVAC.txt', 'IONIZ.txt', 'PHONON.txt', 'E2RECOIL.txt', 'LATERAL.txt', 'TDATA.txt')
$procs = @()
for ($k = 0; $k -lt $Cases.Count; $k++) {
  $wdir = Join-Path $work ('h' + ($k + 1))
  foreach ($sub in 'SRIM Outputs', 'SRIM Restore') { Get-ChildItem (Join-Path $wdir $sub) -File -ErrorAction SilentlyContinue | Remove-Item -Force }
  foreach ($f in @('TRIM.DAT', 'TRIM.dat') + $rootFiles) { $p = Join-Path $wdir $f; if (Test-Path $p) { Remove-Item -Force $p } }
  Copy-Item -Force (Join-Path (Join-Path $root 'trim_cases3') ($Cases[$k] + '\*')) $wdir
  Set-Content -Path (Join-Path $wdir 'TRIMAUTO') -Value @('1', '') -Encoding ASCII
  Write-Host ("Starting {0} in {1}" -f $Cases[$k], $wdir)
  $procs += [pscustomobject]@{ Case = $Cases[$k]; Dir = $wdir; Start = Get-Date; Proc = (Start-Process -FilePath (Join-Path $wdir 'TRIM.exe') -WorkingDirectory $wdir -PassThru) }
}
foreach ($p in $procs) {
  $p.Proc.WaitForExit()
  $dest = Join-Path $out $p.Case
  New-Item -ItemType Directory -Force -Path $dest | Out-Null
  Get-ChildItem $dest -File | Remove-Item -Force
  foreach ($sub in 'SRIM Outputs', 'SRIM Restore') { Get-ChildItem (Join-Path $p.Dir $sub) -File -ErrorAction SilentlyContinue | Copy-Item -Destination $dest -Force }
  foreach ($f in $rootFiles) { $q = Join-Path $p.Dir $f; if (Test-Path $q) { Copy-Item -Force $q $dest } }
  Copy-Item -Force (Join-Path $p.Dir 'TRIM.IN') $dest
  Copy-Item -Force (Join-Path $p.Dir 'TRIM.DAT') $dest -ErrorAction SilentlyContinue
  if (-not (Test-Path (Join-Path $dest 'RANGE.txt'))) { Set-Content -Path (Join-Path $dest 'RANGE.txt') -Value 'missing: TRIM wrote no RANGE.txt' }
  Set-Content -Path (Join-Path $dest '_done.txt') -Value ("finished=True seconds={0:F1} worker={1} end={2} runner=held" -f ((Get-Date) - $p.Start).TotalSeconds, (Split-Path -Leaf $p.Dir), (Get-Date -Format o))
  Write-Host ("Finished {0}" -f $p.Case)
}
Write-Host "Done."
