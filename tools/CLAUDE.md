# tools/ — validation, data generation, build helpers

| Folder | What | Runs on |
|---|---|---|
| `trim_compare/` | `trimout.py` (readers for TRIM's output files), `compare.py` (one TRIM run vs one AMCSET run), `batch.py` (run AMCSET with several seeds on every finished TRIM case), `analyze.py` (statistics over all runs of each input: z per quantity, χ² of profiles) | any Python 3 (standard library only) |
| `windows/` | `build_and_test.bat`, `compare_with_trim.bat` (batch + analyze over the reference set), `benchmark_vs_trim.bat`/`.ps1` (controlled speed test), `gui_smoke_test.bat`, `setup_runner.ps1` (registers the self-hosted CI runner); for releases `make_msix.ps1` (MSIX from a build, Windows SDK) and `sign.ps1` (signtool with Azure Artifact Signing or a certificate; docs/RELEASING.md) | the Windows workstation; the release workflow |
| `linux/` | `make_appimage.sh` (single-file Linux AppImage: builds in `~/amcset-appimage`, or packages `AMCSET_BUILD_DIR`; lays the files out with `cmake --install` and bundles Qt with linuxdeploy). Desktop entry, AppStream data and icons are in `packaging/linux` and `resources/icons` | Linux or WSL; the release workflow |
| `srim/` | SRIM-2013 stopping-table extraction (SR Module), packing into `srim2013_stopping.amcset`, TRIM reference-case generation (`gen_cases.py`, `gen_cases3.py`, `gen_cases4.py`: Ni ions) and batch runs of TRIM.exe | Windows, inside the SRIM installation; see `srim/README.md` |
| `damage_models/` | `pka_damage_functions.py` (PKAs of one energy in pure targets: TRIM full-cascade vacancies and replacements, NRT, Kinchin-Pease, arc-dpa, rpa per PKA), `arc_dpa_vs_srim_postprocessing.py` (the dpa table against Mitsi et al. 2023's SRIM-2013 arc-dpa); docs/VALIDATION.md §5 | any Python 3 (standard library only) |
| `validation/` | `srim_table_scan.py` (every ion's SRIM-2013 electronic stopping against its neighbours' at the same velocity, in nine targets: finds Ni and Au), `stopping_reference.py` (electronic stopping against the measurements of the IAEA database and NIST PSTAR/ASTAR; reads those downloads, which are not in the repository), `literature_cases.py` (AMCSET against measurements, the values and sources in the script and in `range_data.py`, `sputter_data.py`: Ni ranges and peak depths with and without `physics.stopping_corrections.nickel_ions`, 350 ranges in Si, SiO2, C, B, SiC, sputtering yields, published SRIM damage calculations, As channeling tails in Si, 50 keV Sb in Si, ionization yields of Ge and Si recoils); runs cached in `<out>/runs`, structured results in `<out>/results.json`; docs/VALIDATION.md §2.4, §3, §4, §6, §7. `make_figures.py` draws the 35 figures of the review `docs/validation/` from all these outputs (and the TRIM comparison, damage-model and crystal outputs) into `docs/validation/figures/` | any Python 3 (standard library only; `make_figures.py` needs matplotlib); a minute (stopping) and about half an hour (every literature case) |
| `crystal/` | `channeling_checks.py`: the crystal transport's checks (directions and parameters for B in Si, the channeling dip of He in W, the normalisation of the stopping division, random-orientation averages, the cut-off at Ed against NRT). `calibration.py`: the transport against measurements (the dip of 1 MeV He in Ge ⟨110⟩; B in Si at 0.5, 15 and 30 keV for several stopping divisions), with the measured values and their sources in the script; runs are cached in `<out>/runs`. `polycrystal_checks.py`: grains of several sizes, textures and boundary widths against the single crystal and the amorphous target, and grains from a grain file. `range_study.py`: why ranges in crystals are longer (seven ion-crystal pairs along axes, planes, random directions and orientations, static lattice, tilt scans, 0.5 keV B against measured profiles). All write `results.txt`, the source of the numbers in docs/CLI.md "Crystals" | any Python 3 (standard library only); about 2, 6, 1 and 5 minutes on the workstation |
| `gpu_compare/` | `compare_backends.py`: AMCSET on the CPU against GPU backends (`--backends opencl:mixed,hip:mixed,...`), seven cases, independent seeds, z per quantity of summary.json against Student's t (`--resummarize` recomputes a summary); `speed.py`: wall time and compute report per case and backend (`--live BACKEND`: the first N ions of one run, for the live-plot animation); `make_figures.py`: the plots and animations of docs/GPU.md into `docs/gpu/figures/` (the SVG diagrams there are drawn by hand); docs/GPU.md "CPU against GPU", "Speed", "Figures" | any Python 3 (standard library only; `make_figures.py` needs matplotlib and numpy); minutes to an hour (uses the GPU: ask first, it may be shared) |
| `perf/` | `bench.py` (speed and bit-identity of a SRIM-mode change over 13 reference cases: `run`, `compare`, `table`, `ab`), `fast_vs_exact.py` (statistical check of `physics.fast_math`), `fuzz_inputs.py` (input robustness), `winsample.cpp` (sampling profiler for Windows); docs/DEVELOPMENT.md §5, §7 | any Python 3 (standard library only); winsample: Windows |
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
  Windows Python). The exceptions are `validation/make_figures.py` and
  `gpu_compare/make_figures.py`, which need matplotlib (and numpy); they
  only draw figures and compute no evidence.
- After rerunning a validation script, redraw the review's figures
  (`python tools/validation/make_figures.py`) and check the numbers quoted
  in `docs/validation/*.md` against the new output.
- TRIM reference data and SR Module output live outside the repository
  (`<SRIM>\_amcset\results`); never copy them into it. The one exception is
  the packed stopping file in `data/srim2013/`.
- Anything that starts TRIM.exe or SR Module opens windows and runs for hours:
  ask the user first.
