# tools/srim — SRIM-2013 data extraction and TRIM reference runs

SRIM-2013 is free for non-commercial use. AMCSET bundles one product of
these scripts, the stopping tables in `data/srim2013/` (used by default;
see its README); TRIM's own outputs are not distributed. These scripts
generate, from a local SRIM-2013 installation (Windows), the stopping tables
and the TRIM reference results used to validate SRIM mode
(`docs/VALIDATION.md` §2). You only need them to regenerate the tables or to
add reference cases.

Copy this folder into the SRIM installation as a subfolder (for example
`C:\SRIM-2013\_amcset`); the PowerShell scripts look for `TRIM.exe` and
`SR Module` in the parent folder.

| Step | Script | Result |
|---|---|---|
| 1 | `python gen_cases.py <out>` (env `ATOMDATA=<SRIM>\Data\ATOMDATA`) | `sr_jobs.csv` (projectile/target pairs), `sr_energies_keV.txt`, `trim_cases\NNN_name\TRIM.IN` (reference cases), `atomdata.json` |
| 2 | `extract_stopping_tables.bat` (`.ps1 -Workers N`) | `results\stopping\sr_ZZZ_ZZZ_s.txt`, one SR Module table per pair (parallel copies of SR Module; re-running skips finished tables) |
| 3 | `python pack_stopping.py results/stopping results/srim2013_stopping.amcset.gz` | one gzip'ed table file; AMCSET reads it uncompressed: `gunzip` it and replace `data/srim2013/srim2013_stopping.amcset` (or pass it as `--stopping-table` / `stopping.table`) |
| 4 | `run_trim_cases.bat` (`.ps1 -Workers N -MaxProcs N -Filter '0*'`) | `results\trim\<case>\` with every TRIM output file (TRIM is run in batch mode via TRIMAUTO in private copies of the installation) |
| 5 | `python gen_cases3.py atomdata.json trim_cases3` | batch 3 of the reference cases (seed replicates, layers/alloys, Ni, gases, TRIM.DAT inputs, options) |
| 6 | `run_trim_cases3.bat` (`.ps1`, same options) | `results\trim3\<case>\` |
| 7 | `run_held_cases.bat` (`.ps1 -Cases ...`) | the TRIM.DAT recoil cases (codes 6/7) that step 6 skips: TRIM shows a message box for them even in batch mode, which the user must answer; results in `results\trim3\<case>\` |
| 8 | `python gen_cases4.py atomdata.json trim_cases4`, then `run_trim_cases3.ps1 -CaseSet trim_cases4 -OutSet trim4` | batch 4: Ni ions in Fe, Ni, Mo and Si at 1-11.5 MeV with Co and Cu controls (the Ni stopping of SRIM-2013, `docs/VALIDATION.md` §3, §9); results in `results\trim4\<case>\` |

On the workstation the installation is
`C:\Users\natha\Applications\SRIM-2013-Pro\_amcset` and the reference
results are in `results\trim` and `results\trim3` there (161 TRIM runs of
75 inputs, the set behind `docs/VALIDATION.md` §2.2).

Then compare with AMCSET, one case:

```
amcset trim results\trim\<case>\TRIM.IN -o amcset_out\<case>
python <repo>\tools\trim_compare\compare.py results\trim\<case> amcset_out\<case>
```

or every case with several seeds and statistics:
`<repo>\tools\windows\compare_with_trim.bat` (`docs/DEVELOPMENT.md` §6).

Notes:

- SR Module and TRIM show a window per run; the scripts start several at a
  time. Close the console window to stop; finished outputs are kept and
  re-running resumes.
- `pack_stopping.py` reads `sr_jobs.csv` from the current directory.
- A full table set (92 projectiles × 92 targets, solid and gas phases) has
  ~17000 tables and takes several hours to extract. TRIM uses gas-phase
  stopping for every element of a gas layer (e.g. the carbon of CO2), so the
  gas tables cover all elements; an older table set with gas tables for the
  natural gases only is completed by re-running steps 1-3 (finished tables
  are skipped).
- Do not commit TRIM outputs or other table files; `.gitignore` excludes
  `*.amcset`, `*.amcset.gz` and `srim_tables/` except the bundled
  `data/srim2013/srim2013_stopping.amcset`. Replacing that file changes
  every SRIM-mode result: rerun the TRIM comparison and say so in the
  CHANGELOG.
