# tools/ — validation, data generation, build helpers

| Folder | What | Runs on |
|---|---|---|
| `trim_compare/` | `trimout.py` (readers for TRIM's output files), `compare.py` (one TRIM run vs one AMCSET run), `batch.py` (run AMCSET with several seeds on every finished TRIM case), `analyze.py` (statistics over all runs of each input: z per quantity, χ² of profiles) | any Python 3 (standard library only) |
| `windows/` | `build_and_test.bat`, `compare_with_trim.bat` (batch + analyze over the reference set), `benchmark_vs_trim.bat`/`.ps1` (controlled speed test), `gui_smoke_test.bat`, `make_installer.ps1` (single-file installer, `installer/` holds its scripts) | the Windows workstation |
| `srim/` | SRIM-2013 stopping-table extraction (SR Module), packing into `srim2013_stopping.amcset`, TRIM reference-case generation (`gen_cases.py`, `gen_cases3.py`) and batch runs of TRIM.exe | Windows, inside the SRIM installation; see `srim/README.md` |
| `damage_models/` | `pka_damage_functions.py` (PKAs of one energy in pure targets: TRIM full-cascade vacancies and replacements, NRT, Kinchin-Pease, arc-dpa, rpa per PKA), `arc_dpa_vs_srim_postprocessing.py` (the dpa table against Mitsi et al. 2023's SRIM-2013 arc-dpa); docs/SRIM_MODE.md §6.3 | any Python 3 (standard library only) |
| `crystal/` | `channeling_checks.py`: the crystal transport's checks (directions and parameters for B in Si, the channeling dip of He in W, the normalisation of the stopping division, random-orientation averages, the cut-off at Ed against NRT). `calibration.py`: the transport against measurements (the dip of 1 MeV He in Ge ⟨110⟩; B in Si at 0.5, 15 and 30 keV for several stopping divisions), with the measured values and their sources in the script; runs are cached in `<out>/runs`. `polycrystal_checks.py`: grains of several sizes, textures and boundary widths against the single crystal and the amorphous target, and grains from a grain file. All write `results.txt`, the source of the numbers in docs/CLI.md "Crystals" | any Python 3 (standard library only); about 2, 6 and 1 minutes on the workstation |
| `determinism_harness/` | Prints every collision point of fixed-seed native-engine histories, to prove a refactor is bit-identical | any; `docs/TESTING.md` |

## The validation loop (SRIM mode vs TRIM)

```bat
:: all reference cases, AMCSET 8 seeds each (12 for cases with TRIM replicates), then statistics
tools\windows\compare_with_trim.bat [SRIM folder] [TRIM cases folder]
:: default cases: %USERPROFILE%\Applications\SRIM-2013-Pro\_amcset\results\trim3
:: output: build\trim_validation\<cases folder name>\{<case>\<seed>\, cmp\, progress.txt, analysis.txt, analysis.json}
```

Re-running skips finished runs: move or rename the output folder to start
a fresh comparison (earlier rounds are kept as `trim3_v1`, `trim3_v2`, …).
`analysis.txt` starts with the per-quantity summary (units, |z| < 2, |z| < 3,
χ²/n, largest |z|) followed by every unit; units are TRIM runs of the same
input. A unit with |z| > 3 in several correlated quantities is usually one
run-to-run fluctuation; a real difference persists with more seeds and in
related cases.

One case by hand:

```bat
amcset trim <case>\TRIM.IN -s 101 -q --formats trim -o out\101
python tools\trim_compare\compare.py <case> out\101
```

## Rules

- These scripts are part of the evidence chain for `docs/SRIM_MODE.md`;
  keep `compare.py`/`analyze.py` definitions stable, or rerun and update the
  documented numbers when you change them.
- Python tools use only the standard library (they run with the plain
  Windows Python).
- TRIM reference data and SR Module output live outside the repository
  (`<SRIM>\_amcset\results`); never copy them into it. The one exception is
  the packed stopping file in `data/srim2013/`.
- Anything that starts TRIM.exe or SR Module opens windows and runs for hours:
  ask the user first.
