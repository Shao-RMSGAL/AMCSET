# Validation of AMCSET

Every comparison of AMCSET's results with an independent reference in one
place: SRIM-2013's TRIM.exe (which SRIM mode reproduces), measurements
(stopping, ranges, channeling, ionization), and molecular dynamics (damage).
Each section gives the cases, the parameters, the script that reruns them and
where its output is, then what agrees, what deviates and why. Numbers are
quoted from the scripts' output files named in each section.

**The review with figures**: [docs/validation/](validation/README.md)
explains these comparisons chapter by chapter with 35 figures (drawn by
`tools/validation/make_figures.py` from the same output files). This
document stays the reference for the tables and parameters.

## 1. Summary

| What | Reference | Result | Section |
|---|---|---|---|
| SRIM mode against TRIM.exe | 161 SRIM-2013 TRIM runs of 75 inputs; 11 runs of Ni, Co and Cu ions (batch 4) | 1031 of 1032 compared quantities within 3σ; open: replacements in Xe gas (+4.6 %, z = 5.0). Ni ions: 166 of 169 quantities and 21 of 22 profiles within 3σ, the rest in one run (§2.3) | §2 |
| Published SRIM damage calculations | Stoller et al. 2013, Weber and Zhang 2019 | vacancies per ion 0.93-1.01 of theirs, Stoller's Fe cases and damage energies to 0-2 % | §2.4 |
| SRIM-2013 electronic stopping tables | the tables' own Z1 trend | Ni ions up to 39 % below their neighbours Co and Cu at 0.02-0.3 MeV/u, in every target; Au ions up to 26 % below Pt and Hg; no other ion | §3.1 |
| Electronic stopping against measurements | IAEA database: 38 268 measurements, 901 ion-element pairs; NIST PSTAR/ASTAR | SRIM tables: median 0.996, 81 % within 10 %, 94 % within 20 %; Ni ions 0.898 at 0.025-0.25 MeV/u (1.042 with `nickel_ions`), the data sets disagreeing among themselves; H and He in elements within 3 % of NIST above 0.1 MeV (water 4 %) | §3.2, §3.3 |
| Ranges of Ni ions | measured peak depths, 2-11.5 MeV, in Ni, Fe and two steels (Chauhan et al. 2026, Pena et al. 2025, Song et al. 2026) | SRIM-2013 stopping: peaks 1.2-1.9 × too deep; with `nickel_ions` (§9) six of ten within 9 % (3.5-9 MeV), 5 MeV in Ni and 11.5 MeV in Fe 1.21 ×, 2 MeV 1.2-1.8 ×; TRIM.exe itself has the long ranges (§2.3) | §4.1, §2.3 |
| Ranges of 3-83 ions in Si, SiO2, C, B, SiC | Behar et al. 1995, Oosterhoff 1988, Zhang et al. 2009 (RBS, SIMS, nuclear reactions; 350 points) | Si and SiO2: median 0.99-1.06; heavy ions in C, B, SiC 12-25 % too short (the known SRIM deficit); straggling too small in light targets | §4.3 |
| Range of 50 keV Sb in Si | Titze et al. 2021 (SIMS) | AMCSET 33.5 nm, measured 27.5 ± 1.3 nm (+22 %); straggling 9.8 against 9.7 nm | §4.2 |
| Sputtering yields | Yamamura et al. compilations: Ar on six elements, H, D, He on Fe and W (158 points) | median 1.05 above 0.5 keV, 88 % within a factor 2; per pair 0.68 (Ar on Au) to 1.88 (He on W); no yield near the light-ion thresholds | §4.5 |
| Damage functions | arc-dpa/rpa MD fits, SRIM post-processing, 6 K resistivity | arc-dpa reproduces the published post-processing to 0.43 % rms; TRIM's vacancies overestimate surviving defects 4-19 × | §5 |
| Crystal transport | channeling dip of He in Ge; B profiles in Si; As channeling tails against tilt (Tasch et al. 1994) | dip and peak positions within the data; the stopping division is not fixed by them; As tails 9-20 × too high at 100 nm (perfect crystal) | §6 |
| Ionization of keV recoils | Ge (Bonhomme 2022), Si (Chavarria 2016, Izraelevitch 2017) | Ge 25-32 % below the data, Si 28 % to 3 × above; not fixed by a change of the stopping | §7 |
| Native engine (electrons) | — | not validated | §8 |

The deviations and the corrections that can be switched on: §9. What is not
covered yet: §10.

## How the cases are rerun

| Section | Command (Windows, after `build\Release\generators\conanrun.bat`) | Output |
|---|---|---|
| §2 | `tools\windows\compare_with_trim.bat` (about an hour on all cores) | `build\trim_validation\trim3\analysis.txt` |
| §3.1 | `python tools\validation\srim_table_scan.py --exe build\Release\amcset.exe --out build\validation\table_scan --threshold 0.15` | `results.txt`, `ratios.csv` there |
| §3.2, §3.3 | `python tools\validation\stopping_reference.py --exe build\Release\amcset.exe --iaea StoppingPower.csv --nist DIR --out build\validation\stopping` (a minute) | `results.txt`, `iaea_ratios.csv` there |
| §2.4, §4, §6 (As), §7 | `python tools\validation\literature_cases.py --exe build\Release\amcset.exe --out build\validation\literature` (4 seeds; about half an hour) | `results.txt`, `ranges.csv`, `sputtering.csv` there; runs cached in `runs\` |
| §5 | `tools\damage_models\pka_damage_functions.py`, `arc_dpa_vs_srim_postprocessing.py` | `build\damage_models\` |
| §6 | `tools\crystal\calibration.py`, `range_study.py`, `channeling_checks.py` | `build\crystal_calibration\`, `build\range_study\` |

The scripts use only Python's standard library and the `amcset` executable;
the measured values and their sources are written into the scripts next to
the scenarios that reproduce them (the longer tables in
`tools/validation/range_data.py` and `sputter_data.py`). The stopping
comparison reads two downloads that are not in the repository: the IAEA
database's full download (`StoppingPower.csv`, https://nds.iaea.org/stopping/,
release 2026-09 here) and NIST's PSTAR and ASTAR tables for Si, Fe, Cu, W, Au
and water (default energies, saved as `pstar_<material>.txt` and
`astar_<material>.txt` with NIST's columns).

## 2. Against SRIM-2013 TRIM.exe

Reference runs: SRIM-2013.00 TRIM in batch mode (TRIM.IN, `TRIM.exe`),
1000-20000 ions per case, TRIM defaults except as noted (Ed 40 eV, Eb 0,
Es 4.34 eV for Fe). AMCSET: the same TRIM.IN via `amcset trim`, SRIM
stopping tables, the same number of ions. Comparison tool:
`tools/trim_compare/compare.py` (vacancies are integrated from the depth
tables, because TRIM rounds the totals in VACANCY.txt to integers).

Every number from a Monte Carlo run carries a statistical error, and TRIM's
reference values carry the same error as AMCSET's. To tell a real
discrepancy from noise, each case was also run with four further seeds; the
spread of the five AMCSET runs estimates the standard error σ of one run
(TRIM's or AMCSET's), and the mean of the five runs is compared with TRIM:
z = (mean − TRIM) / (σ √(1 + 1/5)). |z| < 2 is consistent with statistics.

Statistics of the comparison (78 cases, `physics.transport: trim2013`,
November 2026; thin-target cases left out of the summary). With σ estimated
from five runs, z follows a t distribution with 4 degrees of freedom, for
which 88 % of the values lie within ±2 and 96 % within ±3 when the two
programs agree in the mean:

| Quantity | Cases | |z| < 2 | |z| < 3 | Median |Δ| % | Largest |z| |
|---|---|---|---|---|---|
| Mean range | 76 | 65 (86 %) | 72 (95 %) | 0.14 | 4.1 (H → Cr/Fe 1 MeV KP) |
| Range straggling | 76 | 67 (88 %) | 76 (100 %) | 1.32 | 2.4 (He → Ti 1 MeV FC) |
| Vacancies | 76 | 67 (88 %) | 75 (99 %) | 0.38 | -4.7 (H → Fe 3 MeV KP) |
| Replacements | 34 | 29 (85 %) | 31 (91 %) | 0.85 | 5.8 (Fe → Cr/Fe 2 MeV FC) |
| Sputtering yield (≥ 0.01/ion) | 20 | 13 (65 %) | 16 (80 %) | 3.50 | -4.8 (Fe → Fe 2 MeV FC) |
| Ionization by ion | 76 | 70 (92 %) | 74 (97 %) | 0.04 | 5.6 (H → Ti 1 MeV KP) |
| Ionization by recoils | 76 | 67 (88 %) | 72 (95 %) | 0.62 | 4.0 (Fe → Fe 2 MeV KP) |
| Phonons by ion | 76 | 65 (86 %) | 69 (91 %) | 0.14 | 7.3 (Fe → Fe 10 keV FC) |
| Phonons by recoils | 76 | 69 (91 %) | 75 (99 %) | 0.28 | -4.4 (H → Fe 3 MeV KP) |

- **Fe → Fe** (the priority case), 1 keV - 5 MeV: ranges, vacancies (quick
  and full cascades), replacements, the energy partition and the sputtering
  yields agree with TRIM within its statistical scatter. Median deviations:
  range 0.55 %, vacancies 0.13 %.
- **Other transition metals, H and He**: the same; the largest vacancy
  deviation is H → Fe 3 MeV KP (−1.6 %, z −4.7).
- **Sputtering**: yields from 1 keV to 500 keV agree within 1-4 % (Fe → Fe
  1 keV 1.95 vs 1.98 atoms/ion, 10 keV 4.92 vs 5.10; Ar → Cu 1 keV 3.59 vs
  3.55). Yields of MeV ions rest on a few rare
  near-surface cascades per thousand ions; their z values overstate the
  significance, because the five-run σ of such clustered counts is itself
  uncertain.
- Backscattering of MeV ions (0.01-0.1 % of the ions) is counting noise.
- Thin-target cases ("thin") transmit > 99.9 % of the ions; their range
  statistics rest on a handful of stopped ions.

AMCSET's values are the means of the five runs, z as above.

| Case | Ions | Range TRIM (Å) | AMCSET | Δ % (z) | Vacancies/ion TRIM | AMCSET | Δ % (z) | Replacements Δ % (z) | Sputtered/ion TRIM | AMCSET | Δ % (z) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Fe → Fe 1 MeV FC | 2000 | 3545 | 3550 | +0.1 (+0.1) | 8371 | 8362 | -0.1 (-0.3) | -0.2 (-0.4) | 1.81 | 2.13 | +17.7 (+1.5) |
| Fe → Fe 1 MeV FC 60° | 2000 | 1874 | 1901 | +1.4 (+0.9) | 8083 | 8098 | +0.2 (+0.9) | -0.3 (-0.9) | 7.94 | 8.16 | +2.7 (+0.5) |
| Fe → Fe 1 MeV FC Ed 25/Eb 3 | 2000 | 3479 | 3556 | +2.2 (+1.5) | 13585 | 13499 | -0.6 (-1.3) | -0.9 (-1.6) | 2.15 | 2.10 | -2.5 (-0.4) |
| Fe → Fe 1 MeV KP | 10000 | 3530 | 3544 | +0.4 (+1.7) | 4086 | 4079 | -0.2 (-1.5) | – | – | – | – |
| Fe → Fe 1 MeV KP 60° | 10000 | 1892 | 1903 | +0.6 (+1.1) | 4011 | 4002 | -0.2 (-1.8) | – | – | – | – |
| Fe → Fe 1 MeV KP Ed 25/Eb 3 | 10000 | 3560 | 3541 | -0.5 (-4.0) | 6531 | 6538 | +0.1 (+1.0) | – | – | – | – |
| Fe → Fe 100 keV FC | 5000 | 351 | 345 | -1.7 (-1.6) | 1313 | 1315 | +0.1 (+1.3) | +0.0 (+0.1) | 4.73 | 4.89 | +3.5 (+2.6) |
| Fe → Fe 100 keV KP | 10000 | 346 | 345 | -0.6 (-1.0) | 658 | 659 | +0.1 (+2.1) | – | – | – | – |
| Fe → Fe 10 keV FC | 5000 | 53.1 | 54.2 | +2.1 (+3.0) | 147 | 146 | -0.2 (-1.1) | -0.4 (-0.6) | 5.10 | 4.92 | -3.5 (-2.0) |
| Fe → Fe 10 keV KP | 20000 | 54.2 | 53.8 | -0.6 (-1.7) | 75.2 | 75.2 | -0.0 (-0.0) | – | – | – | – |
| Fe → Fe 1 keV FC | 5000 | 14.1 | 14.1 | +0.4 (+0.5) | 15.1 | 15.1 | +0.1 (+0.5) | +2.4 (+2.0) | 1.98 | 1.95 | -1.8 (-1.4) |
| Fe → Fe 1 keV KP | 20000 | 13.5 | 13.5 | +0.1 (+0.2) | 7.72 | 7.71 | -0.1 (-1.0) | – | – | – | – |
| Fe → Fe 1 keV mono | 5000 | 14.0 | 14.0 | +0.2 (+0.3) | 15.0 | 15.1 | +0.4 (+1.2) | – | 1.91 | 1.93 | +1.3 (+1.4) |
| Fe → Fe 2 MeV FC | 2000 | 6726 | 6806 | +1.2 (+1.3) | 12164 | 12023 | -1.2 (-2.1) | -0.7 (-1.1) | 1.64 | 1.33 | -19.0 (-4.8) |
| Fe → Fe 2 MeV KP | 10000 | 6799 | 6787 | -0.2 (-1.4) | 5831 | 5835 | +0.1 (+0.7) | – | – | – | – |
| Fe → Fe 5 MeV FC | 1000 | 13860 | 13834 | -0.2 (-0.5) | 16340 | 16296 | -0.3 (-0.8) | -0.5 (-1.3) | 0.543 | 0.706 | +30.0 (+1.0) |
| Fe → Fe 5 MeV KP | 10000 | 13820 | 13806 | -0.1 (-0.8) | 7890 | 7893 | +0.0 (+0.3) | – | – | – | – |
| Fe → Fe 500 keV FC | 2000 | 1727 | 1759 | +1.9 (+1.2) | 5209 | 5196 | -0.2 (-0.9) | -0.5 (-1.2) | 2.97 | 3.01 | +1.5 (+0.2) |
| Fe → Fe 500 keV KP | 10000 | 1765 | 1749 | -0.9 (-2.5) | 2556 | 2558 | +0.1 (+1.2) | – | – | – | – |
| Fe → Fe-9Cr 2 MeV FC | 1000 | 6843 | 6824 | -0.3 (-0.6) | 11944 | 11988 | +0.4 (+0.6) | +0.0 (+0.1) | 1.05 | 1.29 | +23.1 (+3.1) |
| Fe → Fe-9Cr 2 MeV KP | 10000 | 6817 | 6812 | -0.1 (-0.6) | 5807 | 5801 | -0.1 (-1.0) | – | – | – | – |
| H → Fe 1 MeV FC | 5000 | 64940 | 64837 | -0.2 (-2.2) | 18.7 | 19.4 | +4.1 (+1.9) | +5.6 (+2.1) | – | – | – |
| H → Fe 1 MeV KP | 20000 | 65210 | 65261 | +0.1 (+3.3) | 10.2 | 10.1 | -1.4 (-1.4) | – | – | – | – |
| H → Fe 1 MeV thin KP | 20000 | 18150 | 16965 | -6.5 (-0.3) | 1.48 | 1.55 | +4.8 (+1.2) | – | – | – | – |
| H → Fe 100 keV FC | 5000 | 4553 | 4558 | +0.1 (+0.5) | 9.49 | 9.30 | -2.0 (-2.3) | -0.8 (-0.6) | 0.0002 | 0.00048 | (few atoms) |
| H → Fe 100 keV KP | 20000 | 4546 | 4539 | -0.1 (-0.9) | 4.94 | 4.98 | +0.8 (+1.1) | – | – | – | – |
| H → Fe 10 keV FC | 5000 | 745 | 743 | -0.3 (-0.3) | 3.82 | 3.79 | -0.9 (-0.6) | +0.9 (+0.3) | 0.0066 | 0.00432 | (few atoms) |
| H → Fe 10 keV KP | 20000 | 752 | 751 | -0.2 (-0.6) | 2.06 | 2.06 | -0.0 (-0.0) | – | – | – | – |
| H → Fe 3 MeV FC | 5000 | 3.533e+05 | 3.532e+05 | -0.0 (-0.4) | 32.5 | 32.6 | +0.3 (+0.2) | +0.2 (+0.1) | – | – | – |
| H → Fe 3 MeV KP | 20000 | 3.546e+05 | 3.547e+05 | +0.0 (+1.7) | 17.2 | 16.9 | -1.6 (-4.7) | – | – | – | – |
| He → Fe 1 MeV FC | 5000 | 16210 | 16204 | -0.0 (-0.2) | 146 | 146 | +0.2 (+0.2) | +0.0 (+0.0) | 0.0028 | 0.00144 | (few atoms) |
| He → Fe 1 MeV KP | 20000 | 16200 | 16205 | +0.0 (+0.8) | 74.9 | 75.0 | +0.1 (+0.2) | – | – | – | – |
| He → Fe 100 keV FC | 5000 | 3013 | 3004 | -0.3 (-0.7) | 101 | 101 | -0.3 (-0.4) | +0.8 (+0.8) | 0.022 | 0.0113 | -48.5 (-1.7) |
| He → Fe 100 keV KP | 20000 | 3017 | 3016 | -0.0 (-0.2) | 52.1 | 52.0 | -0.2 (-1.2) | – | – | – | – |
| He → Fe 10 keV FC | 5000 | 458 | 455 | -0.5 (-0.7) | 36.4 | 36.6 | +0.6 (+1.4) | +3.6 (+4.8) | 0.121 | 0.118 | -2.2 (-0.3) |
| He → Fe 10 keV KP | 20000 | 459 | 457 | -0.4 (-0.8) | 19.3 | 19.3 | -0.0 (-0.2) | – | – | – | – |
| He → Fe 2 MeV thin KP | 20000 | 14060 | 13689 | -2.6 (-1.1) | 11.7 | 12.1 | +2.9 (+3.3) | – | – | – | – |
| He → Fe 3 MeV FC | 5000 | 49760 | 49796 | +0.1 (+0.9) | 173 | 172 | -0.6 (-0.5) | -1.0 (-0.6) | 0.0006 | 0.0006 | (few atoms) |
| He → Fe 3 MeV KP | 20000 | 49820 | 49838 | +0.0 (+1.2) | 88.8 | 88.6 | -0.3 (-0.5) | – | – | – | – |
| Ar → Cu 1 keV FC | 5000 | 14.3 | 14.5 | +0.9 (+1.1) | 14.4 | 14.4 | -0.4 (-1.1) | +1.4 (+0.9) | 3.55 | 3.59 | +1.0 (+1.5) |
| Ar → Cu 1 keV mono | 5000 | 14.2 | 14.3 | +0.4 (+0.5) | 14.4 | 14.4 | +0.1 (+0.3) | – | 3.53 | 3.58 | +1.4 (+0.8) |
| Co → Co 1 MeV FC | 1000 | 3213 | 3200 | -0.4 (-0.3) | 9807 | 9783 | -0.2 (-0.6) | -0.3 (-0.6) | 2.87 | 2.83 | -1.3 (-0.6) |
| Co → Co 1 MeV KP | 5000 | 3209 | 3210 | +0.0 (+0.1) | 4388 | 4388 | +0.0 (+0.0) | – | – | – | – |
| Cr → Cr 1 MeV FC | 1000 | 3855 | 3871 | +0.4 (+0.4) | 6883 | 6856 | -0.4 (-0.7) | -0.7 (-0.9) | 1.32 | 1.48 | +11.8 (+4.0) |
| Cr → Cr 1 MeV KP | 5000 | 3865 | 3876 | +0.3 (+0.7) | 3566 | 3565 | -0.0 (-0.2) | – | – | – | – |
| Cu → Cu 1 MeV FC | 1000 | 3196 | 3241 | +1.4 (+1.1) | 11570 | 11544 | -0.2 (-0.6) | -0.9 (-1.4) | 4.85 | 4.07 | -16.0 (-2.4) |
| Cu → Cu 1 MeV KP | 5000 | 3219 | 3236 | +0.5 (+0.8) | 4936 | 4932 | -0.1 (-0.5) | – | – | – | – |
| Fe → 316 SS 2 MeV FC | 1000 | 6722 | 6646 | -1.1 (-2.6) | 12177 | 12341 | +1.3 (+2.3) | +2.2 (+3.0) | 1.59 | 1.27 | -20.1 (-4.0) |
| Fe → Cr/Fe 2 MeV FC | 1000 | 6669 | 6673 | +0.1 (+0.0) | 11109 | 11333 | +2.0 (+2.4) | +6.1 (+5.8) | 1.16 | 1.10 | -4.8 (-0.4) |
| H → Co 1 MeV KP | 10000 | 61730 | 61719 | -0.0 (-0.4) | 10.1 | 10.2 | +0.5 (+0.2) | – | – | – | – |
| H → Cr 1 MeV KP | 10000 | 69300 | 69280 | -0.0 (-0.7) | 9.33 | 9.46 | +1.4 (+1.3) | – | – | – | – |
| H → Cr/Fe 1 MeV KP | 20000 | 65350 | 65380 | +0.0 (+4.1) | 10.2 | 10.1 | -1.1 (-2.6) | – | – | – | – |
| H → Cu 1 MeV KP | 10000 | 68300 | 68236 | -0.1 (-3.0) | 11.1 | 11.2 | +1.5 (+1.6) | – | – | – | – |
| H → Mn 1 MeV KP | 10000 | 70100 | 70093 | -0.0 (-0.3) | 9.97 | 10.0 | +0.5 (+0.3) | – | – | – | – |
| H → Mo 1 MeV KP | 10000 | 67260 | 67233 | -0.0 (-0.9) | 9.40 | 9.45 | +0.6 (+0.6) | – | – | – | – |
| H → Nb 1 MeV KP | 10000 | 75400 | 75467 | +0.1 (+2.2) | 8.93 | 8.69 | -2.6 (-2.5) | – | – | – | – |
| H → Ni 1 MeV KP | 10000 | 63830 | 63788 | -0.1 (-2.8) | 11.5 | 11.8 | +2.5 (+2.9) | – | – | – | – |
| H → Ta 1 MeV KP | 10000 | 63930 | 63855 | -0.1 (-1.9) | 11.7 | 11.8 | +1.3 (+1.4) | – | – | – | – |
| H → Ti 1 MeV KP | 10000 | 1.06e+05 | 1.06e+05 | -0.0 (-2.1) | 8.54 | 8.42 | -1.4 (-1.0) | – | – | – | – |
| H → V 1 MeV KP | 10000 | 83150 | 83148 | -0.0 (-0.2) | 8.64 | 8.70 | +0.8 (+0.5) | – | – | – | – |
| H → W 1 MeV KP | 10000 | 55020 | 54961 | -0.1 (-1.8) | 11.5 | 11.6 | +0.9 (+1.1) | – | – | – | – |
| H → W 5 keV KP | 20000 | 311 | 310 | -0.4 (-0.8) | 0.335 | 0.332 | -0.9 (-0.9) | – | – | – | – |
| H → Zr 1 MeV KP | 10000 | 1.058e+05 | 1.057e+05 | -0.1 (-1.2) | 9.10 | 9.06 | -0.4 (-0.3) | – | – | – | – |
| He → Co 1 MeV FC | 2000 | 15850 | 15839 | -0.1 (-0.2) | 164 | 165 | +0.3 (+0.2) | +0.2 (+0.1) | 0.001 | 0.0008 | (few atoms) |
| He → Co 1 MeV KP | 10000 | 15830 | 15840 | +0.1 (+1.2) | 80.9 | 80.8 | -0.2 (-0.3) | – | – | – | – |
| He → Cr 1 MeV FC | 2000 | 16260 | 16255 | -0.0 (-0.1) | 129 | 132 | +2.1 (+0.9) | +1.5 (+0.8) | 0.0015 | 0.0008 | (few atoms) |
| He → Cr 1 MeV KP | 10000 | 16270 | 16252 | -0.1 (-1.8) | 69.9 | 70.3 | +0.6 (+0.9) | – | – | – | – |
| He → Cu 1 MeV FC | 2000 | 19130 | 19216 | +0.5 (+1.2) | 205 | 200 | -2.4 (-1.4) | -3.8 (-1.4) | – | – | – |
| He → Mn 1 MeV FC | 2000 | 17440 | 17429 | -0.1 (-0.2) | 134 | 134 | +0.6 (+0.3) | +0.7 (+0.3) | 0.002 | 0.0037 | (few atoms) |
| He → Mn 1 MeV KP | 10000 | 17420 | 17417 | -0.0 (-0.4) | 69.4 | 69.7 | +0.6 (+0.7) | – | – | – | – |
| He → Mo 1 MeV FC | 2000 | 15480 | 15463 | -0.1 (-0.3) | 123 | 124 | +1.0 (+0.4) | +2.1 (+0.6) | 0.0005 | 0.0006 | (few atoms) |
| He → Mo 1 MeV KP | 10000 | 15440 | 15453 | +0.1 (+0.8) | 64.4 | 64.4 | -0.0 (-0.0) | – | – | – | – |
| He → Nb 1 MeV FC | 2000 | 15700 | 15754 | +0.3 (+0.9) | 102 | 102 | -0.2 (-0.1) | -0.1 (-0.0) | – | – | – |
| He → Nb 1 MeV KP | 10000 | 15750 | 15743 | -0.0 (-0.5) | 55.8 | 55.8 | +0.0 (+0.1) | – | – | – | – |
| He → Ni 1 MeV FC | 2000 | 16750 | 16745 | -0.0 (-0.1) | 187 | 185 | -1.3 (-0.7) | -1.1 (-0.5) | – | – | – |
| He → Ta 1 MeV FC | 2000 | 16480 | 16450 | -0.2 (-0.3) | 149 | 150 | +0.7 (+0.3) | +1.4 (+0.4) | 0.0005 | 0.0009 | (few atoms) |
| He → Ti 1 MeV FC | 2000 | 23380 | 23370 | -0.0 (-0.2) | 102 | 102 | +0.7 (+0.3) | -0.4 (-0.1) | – | – | – |
| He → V 1 MeV FC | 2000 | 17730 | 17776 | +0.3 (+1.0) | 103 | 103 | -0.1 (-0.1) | +1.0 (+0.3) | 0.0005 | 0.0005 | (few atoms) |

### 2.1 Differences found after the first round, and their causes

The second round (§2.2) traced the differences that were left in the
table above:

- **Alloys and layered targets, full cascades** (replacements +6 % for
  Fe → Cr on Fe, +2 % for 316 SS): TRIM replaces a zero lattice or surface
  binding energy of every target atom after the first by 2 eV and 4 eV
  (0x44efc0, docs/SRIM_MODE.md §2.5), which changes the replacement collisions of the later
  elements. Now within statistics.
- **H → Fe 3 MeV, quick calculation** (vacancies −1.6 %): TRIM adds its
  vacancy and range tables in single precision with a weight of
  1/(bin width), which shifts the totals of long runs by up to ±1.5 %
  (predicted +0.64 % for Fe → Fe 2 MeV, observed +0.63 %). AMCSET
  reproduces the rounding in the TRIM text tables (docs/SRIM_MODE.md §5).
- **Quick calculation, energy to recoils** (−3 %): TRIM books the lattice
  binding energy of each displaced atom in E2RECOIL.txt (0x44bae0).
- **Displacements in the VACANCY.txt header of low-energy runs** (−2 % to
  −16 %): a sputtered atom that was set in motion below Ed counts as a
  displacement as well as a vacancy in TRIM's totals.
- **TRIM.DAT recoils**: code 6 is a full cascade (it had been read as the
  quick calculation: half of TRIM's vacancies) and code 7 moves each atom
  as an ion with monolayer steps.
- **Gas targets with elements that are solids** (He → CO2: range +1.9 %):
  TRIM uses gas stopping for every element of a gas layer; the stopping
  file now has gas tables for all elements (16928 tables).

Other notes:

- **Built-in stopping** differs from SRIM by up to ~25 % near the stopping
  maximum (2x for slow heavy ions). For TRIM-equivalent results use the SRIM
  tables.
- `physics.transport: refined` does not follow TRIM's algorithm: with it,
  sputtering yields are 3-10x below TRIM's, H vacancies 2-10 % low and slow
  light ions (H → W 5 keV) far off.
- TRIM.DAT lists (cascade codes 4-7) follow TRIM.exe's handling as read
  from its code (0x449330, 0x449c10): ions take the TRIM.IN species and
  start at the given point (a zero lateral coordinate is randomised within
  one monolayer, as for the beam); code 5 uses the quick mode's flight rule
  and, like 7, does not restore the random-number state after a cascade.
  TRIM.exe renumbers codes 6 and 7 internally (to 7 and 8, 0x43acb0):
  code 6 goes straight to the primary-recoil step with the given direction
  (a full cascade; RANGE.txt counts the event as an ion that stops where
  the atom starts, and the range table also gets the depth where the atom
  stops), while code 7 runs each atom through the ion loop like a TRIM.DAT
  ion with monolayer steps, so its own losses are booked as the ion's.
  Checked against TRIM: Fe PKAs 0.1-50 keV in Fe, code 6: vacancies/event
  120.8 (TRIM 121.05), code 7: 119.5 (TRIM 119.4), ion ionization 5.72 eV/A
  summed over the table (TRIM 5.61), ion phonons 1.99 (TRIM 1.98; ion steps
  of the free flight path would give 2.09). AMCSET's quick calculation for
  recoils has no TRIM code (TRIM.IN gets code 6). Lines above the TRIM.IN
  energy or outside the target are skipped, as in TRIM.

### 2.2 Second round (November 2026)

Reference runs: the TRIM-2013 cases of `tools/srim/gen_cases.py` (pure
metals, alloys and compounds, layered targets in both orders, gas and solid
phase of the same gas, Ed/Eb/Es variants, oblique incidence, thin targets,
TRIM.DAT ions and recoils, KP, full cascades and the surface-sputtering
mode, 1 keV - 50 MeV), with up to seven TRIM runs of the same input with
different seeds. AMCSET: the Windows build with its default stopping (the
bundled SRIM-2013 tables, section 4), 8 seeds per case (12 for cases with
TRIM replicates), run by `tools/windows/compare_with_trim.bat` and analysed
by `tools/trim_compare/analyze.py`. TRIM runs of one input form a unit
(75 units, 161 TRIM runs, 704 AMCSET runs; the final run added two gas
cascades, Ar -> Ar at 100 keV and Xe -> Xe at 1 MeV). For each quantity, z
is the normal deviate of the two-sample t statistic with the pooled
run-to-run scatter (TRIM's own scatter where it has several runs). Each
depth profile gets a χ² over its bins with a per-bin variance k · (mean), k
taken from AMCSET's seeds; its z is the upper tail of χ² as a normal
deviate.

With z normally distributed, 95 % of the units would be within |z| < 2 and
99.7 % within 3; with the run-to-run scatter estimated from a few runs the
tails are wider (t with 7-17 degrees of freedom). 1031 of the
1032 comparisons of integral quantities are within |z| < 3 (RMS of z
1.10); the one beyond is the replacements of Xe -> Xe gas
(below):

| Quantity | Units | \|z\| < 2 | \|z\| < 3 | χ²/n of z | Largest \|z\| |
|---|---|---|---|---|---|
| Mean range | 74 | 65 | 74 | 1.37 | 2.5 (`lay_fe_fecr_rev_2000keV_fc`) |
| Range straggling | 74 | 72 | 74 | 0.85 | 2.0 (`opt_h_fe_10000keV_kp`) |
| Lateral range | 75 | 70 | 75 | 1.32 | 2.6 (`gas_h_n_1000keV_kp_gas`) |
| Backscattered fraction | 75 | 73 | 75 | 0.77 | 2.4 (`opt_h_w_3000keV_thin10000_kp`) |
| Transmitted fraction | 74 | 74 | 74 | 0.04 | 1.6 (`opt_h_w_3000keV_thin10000_kp`) |
| Vacancies (tables) | 75 | 68 | 75 | 1.55 | 2.6 (`opt_fe_fe_1000keV_kp_window`) |
| Vacancies (header) | 31 | 27 | 31 | 2.13 | 2.6 (`fe_fe_1000keV_fc_ed25eb3`) |
| Displacements (header) | 31 | 27 | 31 | 1.95 | 2.7 (`fe_fe_1000keV_fc_ed25eb3`) |
| Replacements | 38 | 31 | 37 | 2.62 | 5.0 (`gas_xe_xe_1000keV_fc_gas`) |
| Sputtering yield | 36 | 33 | 36 | 1.36 | 2.8 (`aly_fe_fecr50_1000keV_fc`) |
| Ionization by ion | 75 | 72 | 75 | 1.04 | 2.9 (`ni_h_ni_1000keV_fc`) |
| Ionization by recoils | 75 | 69 | 75 | 1.13 | 2.4 (`fe_fe_10keV_fc`) |
| Phonons by ion | 75 | 68 | 75 | 1.18 | 2.8 (`h_fe_3000keV_kp`) |
| Phonons by recoils | 75 | 68 | 75 | 1.37 | 2.6 (`opt_fe_fe_1000keV_kp_window`) |
| Energy to recoils | 75 | 71 | 75 | 1.09 | 2.3 (`lay_fe_fecr_rev_2000keV_fc`) |
| Mean depth of vacancies | 74 | 71 | 74 | 1.21 | 2.5 (`fe_fe_1000keV_fc_ed25eb3`) |
| Vacancy profile (χ² per bin) | 74 | 56 | 69 | median χ²/ndf 1.17 | χ²/ndf 2.21 (`fe_fe_10keV_fc`) |
| Ion range profile (χ² per bin) | 74 | 69 | 72 | median χ²/ndf 0.99 | χ²/ndf 2.00 (`fe_fe_10keV_fc`) |

The largest deviations are single units whose quantities move together
(a deeper range with fewer vacancies, as in `fe_fe_1000keV_fc_ed25eb3`
and `lay_fe_fecr_rev_2000keV_fc`): that is how one run-to-run fluctuation
of the range looks, and the other units of the same ions and energies do
not show it. Sputtering yields of MeV ions come from a few near-surface
cascades and are heavy-tailed; the single TRIM runs of the alloy cases
scatter by ±20 % on both sides. The vacancy profiles of five units have
a χ² tail beyond 3σ (χ²/ndf 1.7-2.2). Three are light ions with few
vacancies per ion (H 10 MeV, the thin H → W and He → Fe targets): their
vacancies per bin come from rare hard collisions, for which a variance
proportional to the bin value is too small. The other two (Fe → Fe 10 keV,
Fe → Fe 1 MeV with Ed 25 eV) differ by a few standard errors in single bins
of the deep tail.

**Open: replacements in Xe gas.** `gas_xe_xe_1000keV_fc_gas` (Xe -> Xe gas
at 1 atm, 1 MeV, full cascades; one TRIM run of 500 ions, which took TRIM
8 hours) has replacements +4.6 % (z = +5.0), vacancies +0.7 % (z = +2.1),
displacements +0.7 % (z = +2.0) and recoil ionization -1.6 % (z = -2.0); ranges,
sputtering and the depth profiles agree. The same comparison for Ar -> Ar
gas at 100 keV, with the same density and zero Eb and Es, agrees
(replacements -0.8 %). A gas monolayer is 33 Å, so a fast Xe recoil strikes
an atom every 33 Å over millions of Å and the recoil stack grows to 10^5
pending atoms (135 000 in 30 ions), where a solid needs hundreds. Checked
and ruled out so far: TRIM's limit on the stack (TRIM.cfg, 200 000 in the
SRIM installation used; not reached in 100 AMCSET ions), the R <= 200 clamp
and 14-step Newton iteration of Tmagic's closest approach (negligible for
these impact parameters), and the zero-energy defaults. Not yet excluded:
TRIM keeps its recoil slots in single precision (Variant type 4, `fstp
dword` at 0x435df2), which at depths of 10^6 Å resolves positions only to
0.25 Å. A second TRIM run of this input would also tell whether the single
run is unusual.

AMCSET minus TRIM in % (z) per unit (units with one TRIM run compare
with that run):

| Case | TRIM runs | AMCSET seeds | Range Δ % (z) | Vacancies Δ % (z) | Replacements Δ % (z) | Sputtered Δ % (z) | Ion ionization Δ % (z) | Recoil phonons Δ % (z) |
|---|---|---|---|---|---|---|---|---|
| `aly_fe_fecr50_1000keV_fc` | 1 | 8 | +1.4 (+1.0) | -0.7 (-1.6) | -1.4 (-1.9) | -24.4 (-2.8) | +1.9 (+1.8) | -0.7 (-1.5) |
| `aly_fe_fecr50_1000keV_kp` | 1 | 8 | +0.0 (+0.1) | +0.1 (+0.9) | – | – | +0.0 (+0.1) | +0.1 (+0.9) |
| `aly_fe_feni50_1000keV_fc` | 1 | 8 | -0.0 (-0.0) | -0.3 (-0.6) | -0.2 (-0.3) | +20.3 (+1.8) | +0.2 (+0.3) | -0.3 (-0.6) |
| `ar_cu_1keV_fc` | 4 | 12 | +0.2 (+0.4) | -0.1 (-0.3) | -2.5 (-2.5) | -0.4 (-0.5) | -0.4 (-0.7) | +0.0 (+0.2) |
| `ar_cu_1keV_sputt` | 4 | 12 | -0.9 (-2.1) | +0.1 (+0.6) | – | +1.1 (+1.4) | -0.0 (-0.1) | +0.0 (+0.0) |
| `dat_fe_fe_ions_fc` | 1 | 8 | +1.7 (+1.7) | +0.1 (+0.3) | +0.0 (+0.0) | -14.3 (-1.0) | -0.0 (-0.0) | +0.2 (+0.4) |
| `dat_fe_fe_ions_kp` | 1 | 8 | -1.2 (-1.8) | +0.3 (+1.3) | – | – | +0.2 (+0.2) | +0.3 (+1.3) |
| `dat_fe_fe_pka_fc` | 1 | 8 | -0.0 (-0.1) | +0.3 (+1.1) | +0.1 (+0.2) | – | +1.1 (+0.6) | +0.0 (+0.1) |
| `dat_fe_fe_pka_kp` | 1 | 8 | – | +0.4 (+1.9) | +1.6 (+2.4) | – | – | -0.1 (-1.0) |
| `dat_fecr_pka_fc` | 1 | 8 | +0.0 (+0.5) | -0.2 (-0.8) | -0.4 (-0.4) | – | -1.6 (-0.5) | -0.1 (-0.6) |
| `fe_316ss_2000keV_fc` | 6 | 12 | +0.3 (+1.0) | -0.3 (-1.4) | -0.4 (-1.2) | +7.6 (+1.3) | +0.1 (+0.6) | -0.3 (-1.3) |
| `fe_crfe_2000keV_fc` | 6 | 12 | +0.5 (+1.3) | -0.4 (-1.6) | -0.5 (-1.9) | +5.1 (+0.7) | +0.5 (+1.8) | -0.4 (-1.6) |
| `fe_fe_1000keV_fc` | 3 | 12 | +0.2 (+0.3) | -0.2 (-1.1) | -0.4 (-1.5) | +5.0 (+1.2) | +0.1 (+0.3) | -0.2 (-0.9) |
| `fe_fe_1000keV_fc_ed25eb3` | 3 | 12 | +1.9 (+2.4) | -0.6 (-2.6) | -0.8 (-2.9) | +2.5 (+0.6) | +1.3 (+2.6) | -0.6 (-2.6) |
| `fe_fe_1000keV_kp_ed25eb3` | 5 | 12 | -0.2 (-1.2) | +0.0 (+0.7) | – | – | -0.2 (-1.1) | +0.0 (+0.8) |
| `fe_fe_100keV_fc` | 5 | 12 | -0.5 (-1.1) | +0.0 (+1.3) | -0.0 (-0.1) | +2.2 (+2.2) | -0.3 (-0.7) | +0.0 (+2.1) |
| `fe_fe_100keV_kp` | 3 | 12 | -0.3 (-0.8) | +0.1 (+2.5) | – | – | -0.4 (-1.4) | +0.1 (+2.4) |
| `fe_fe_10keV_fc` | 5 | 12 | +1.1 (+2.4) | -0.1 (-1.5) | -0.0 (-0.1) | -2.1 (-1.8) | +0.7 (+1.6) | +0.0 (+0.5) |
| `fe_fe_1keV_fc` | 5 | 12 | +0.5 (+1.3) | +0.0 (+0.0) | +0.3 (+0.4) | -0.1 (-0.2) | +0.4 (+1.0) | -0.1 (-0.6) |
| `fe_fe_1keV_sputt` | 5 | 12 | -0.0 (-0.0) | -0.1 (-0.6) | – | -0.3 (-0.3) | +0.4 (+0.9) | +0.1 (+1.3) |
| `fe_fe_2000keV_fc` | 7 | 12 | +0.3 (+1.0) | -0.3 (-1.4) | -0.5 (-2.2) | +4.5 (+0.8) | +0.2 (+1.0) | -0.2 (-1.3) |
| `fe_fe_500keV_fc` | 3 | 12 | +1.5 (+2.0) | -0.4 (-1.8) | -0.5 (-1.7) | +1.8 (+0.5) | +1.1 (+1.9) | -0.3 (-1.6) |
| `fe_fe_500keV_kp` | 3 | 12 | -0.2 (-0.7) | +0.1 (+1.2) | – | – | -0.1 (-0.5) | +0.1 (+1.2) |
| `gas_ar_ar_100keV_fc_gas` | 1 | 8 | +1.0 (+0.6) | -0.4 (-1.5) | -0.8 (-0.7) | -14.0 (-1.1) | +0.4 (+0.3) | -0.4 (-1.6) |
| `gas_h_he_100keV_kp_gas` | 1 | 8 | -0.1 (-0.8) | +0.5 (+1.2) | – | – | -0.0 (-1.2) | +0.2 (+0.9) |
| `gas_h_he_100keV_kp_solidph` | 1 | 8 | +0.1 (+1.0) | +0.7 (+1.6) | – | – | -0.0 (-1.2) | +0.3 (+1.3) |
| `gas_h_n_1000keV_kp_gas` | 1 | 8 | -0.0 (-2.4) | +2.6 (+1.7) | – | – | -0.0 (-1.0) | +1.9 (+1.8) |
| `gas_h_n_1000keV_kp_solidph` | 1 | 8 | +0.0 (+1.4) | -1.2 (-0.9) | – | – | +0.0 (+0.1) | -0.8 (-0.9) |
| `gas_he_alaral_5000keV_kp_mixed` | 1 | 8 | -0.0 (-0.0) | +1.0 (+0.6) | – | – | -0.0 (-0.9) | +0.9 (+0.7) |
| `gas_he_ar_1000keV_fc_gas` | 1 | 8 | -0.1 (-0.5) | -1.7 (-0.9) | – | – | +0.0 (+0.4) | -1.4 (-0.8) |
| `gas_he_ar_1000keV_kp_gas` | 1 | 8 | -0.1 (-1.0) | +1.5 (+1.9) | – | – | -0.0 (-1.2) | +1.3 (+1.9) |
| `gas_he_ar_1000keV_kp_solidph` | 1 | 8 | -0.0 (-0.1) | +0.4 (+0.6) | – | – | -0.0 (-1.0) | +0.4 (+0.7) |
| `gas_he_co2_2000keV_kp_gas` | 1 | 8 | +0.0 (+0.3) | -0.7 (-0.9) | – | – | +0.0 (+0.3) | -0.5 (-0.9) |
| `gas_n_n_500keV_kp_gas` | 1 | 8 | +0.1 (+1.3) | -0.3 (-0.8) | – | – | +0.0 (+0.3) | -0.3 (-0.8) |
| `gas_xe_xe_1000keV_fc_gas` | 1 | 8 | -0.7 (-0.3) | +0.7 (+2.1) | +4.6 (+5.0) | -12.0 (-1.0) | +1.7 (+0.7) | +0.3 (+0.9) |
| `h_crfe_1000keV_kp` | 4 | 12 | -0.0 (-1.1) | +0.8 (+1.6) | – | – | -0.0 (-0.5) | +0.5 (+1.5) |
| `h_fe_1000keV_fc` | 5 | 12 | -0.0 (-0.2) | -0.3 (-0.2) | +0.6 (+0.4) | – | +0.0 (+1.5) | -0.3 (-0.3) |
| `h_fe_1000keV_kp` | 7 | 12 | +0.0 (+1.1) | -0.1 (-0.3) | – | – | -0.0 (-0.1) | -0.1 (-0.2) |
| `h_fe_3000keV_kp` | 7 | 12 | +0.0 (+1.9) | -0.2 (-0.3) | – | – | +0.0 (+0.3) | -0.1 (-0.3) |
| `h_nb_1000keV_kp` | 4 | 12 | -0.0 (-0.5) | +1.0 (+1.4) | – | – | +0.0 (+1.0) | +0.8 (+1.3) |
| `h_ni_1000keV_kp` | 5 | 12 | +0.0 (+0.2) | +0.4 (+0.7) | – | – | +0.0 (+1.1) | +0.4 (+0.8) |
| `he_fe_10keV_fc` | 5 | 12 | -0.3 (-0.8) | +0.2 (+0.8) | +0.7 (+1.2) | -8.4 (-1.6) | +0.1 (+0.4) | +0.2 (+0.8) |
| `lay_fe_cr_2000keV_fc` | 2 | 12 | -0.2 (-0.3) | -0.5 (-1.2) | -0.7 (-1.4) | -1.5 (-0.1) | +0.1 (+0.2) | -0.4 (-1.1) |
| `lay_fe_crfe_2000keV_kp` | 1 | 8 | -0.3 (-0.6) | -0.0 (-0.1) | – | – | +0.0 (+0.1) | -0.0 (-0.1) |
| `lay_fe_crfe_order_2000keV_fc` | 2 | 12 | -0.2 (-0.3) | -0.2 (-0.6) | -0.8 (-1.9) | -3.3 (-0.2) | -0.3 (-0.5) | -0.1 (-0.4) |
| `lay_fe_fe2el_2000keV_fc` | 1 | 8 | +1.2 (+1.6) | -0.9 (-1.6) | -1.1 (-1.6) | -10.6 (-2.4) | +1.2 (+1.9) | -0.9 (-1.6) |
| `lay_fe_fecr_rev_2000keV_fc` | 2 | 12 | +1.6 (+2.5) | -1.1 (-2.5) | -1.2 (-2.5) | -0.4 (-0.0) | +1.2 (+2.4) | -1.0 (-2.4) |
| `lay_fe_fefe_2000keV_fc` | 2 | 12 | -0.1 (-0.2) | -0.2 (-0.6) | -0.0 (-0.0) | -4.4 (-0.6) | +0.1 (+0.4) | -0.1 (-0.5) |
| `ni_fe_ni_1000keV_fc` | 1 | 8 | +0.6 (+0.4) | -0.6 (-1.2) | -1.0 (-1.6) | +4.2 (+0.5) | +0.9 (+0.9) | -0.6 (-1.2) |
| `ni_h_ni_1000keV_fc` | 1 | 8 | +0.1 (+0.8) | -4.1 (-1.4) | -6.0 (-1.5) | -75.0 (-1.2) | +0.0 (+2.9) | -3.1 (-1.4) |
| `ni_he_ni_100keV_kp` | 1 | 8 | +0.1 (+0.7) | +0.4 (+2.0) | – | – | -0.0 (-0.7) | +0.4 (+2.1) |
| `ni_ni_ni_100keV_fc` | 1 | 8 | -1.3 (-1.2) | -0.0 (-0.7) | -0.7 (-2.2) | +4.9 (+1.6) | -0.8 (-0.8) | -0.0 (-0.2) |
| `opt_ar_cu_1keV_fc_75deg` | 1 | 8 | -0.0 (-0.0) | +0.3 (+0.4) | +3.1 (+1.5) | -0.4 (-0.4) | +2.0 (+1.4) | +0.5 (+0.6) |
| `opt_ar_cu_5keV_fc_45deg` | 1 | 8 | +1.1 (+1.0) | +0.2 (+0.5) | -0.2 (-0.2) | +0.3 (+0.2) | +0.6 (+0.8) | -0.1 (-0.3) |
| `opt_au_c_1000keV_fc` | 1 | 8 | +0.5 (+0.9) | -0.6 (-1.6) | -0.8 (-0.9) | -14.5 (-1.8) | +1.0 (+1.1) | -0.5 (-1.8) |
| `opt_fe_316ss_2000keV_kp_bragg09` | 1 | 8 | -0.4 (-2.1) | +0.1 (+0.6) | – | – | -0.2 (-1.2) | +0.1 (+0.6) |
| `opt_fe_fe_1000keV_kp_window` | 1 | 8 | +0.7 (+2.5) | -1.4 (-2.6) | – | – | +0.4 (+1.1) | -1.4 (-2.6) |
| `opt_fe_fe_100keV_fc_coll2` | 1 | 8 | -1.7 (-0.1) | -0.9 (-1.2) | -3.3 (-1.4) | -30.6 (-0.8) | +0.9 (+0.1) | -0.3 (-0.3) |
| `opt_fe_fe_100keV_fc_eb3` | 1 | 8 | -0.8 (-0.8) | +0.1 (+2.4) | +0.5 (+1.7) | +2.8 (+0.8) | -0.7 (-0.8) | +0.1 (+2.5) |
| `opt_fe_fe_10keV_fc_es0ed20` | 1 | 8 | -0.3 (-0.5) | -0.3 (-1.9) | +0.2 (+0.4) | -0.9 (-0.6) | -0.7 (-0.9) | -0.1 (-0.9) |
| `opt_fe_fe_1keV_sputt_coll2` | 1 | 8 | +0.3 (+0.1) | -2.2 (-1.1) | – | +0.9 (+0.1) | -2.1 (-0.5) | -0.0 (-0.0) |
| `opt_fe_fe_50000keV_kp` | 1 | 8 | -0.2 (-1.4) | +0.1 (+0.5) | – | – | -0.0 (-1.4) | +0.1 (+0.5) |
| `opt_h_fe_10000keV_kp` | 1 | 8 | +0.1 (+2.2) | -3.0 (-1.3) | – | – | -0.0 (-0.2) | -2.3 (-1.3) |
| `opt_h_si_100keV_kp` | 1 | 8 | +0.1 (+0.9) | -0.0 (-0.0) | – | – | -0.0 (-0.4) | -0.0 (-0.1) |
| `opt_h_w_100keV_kp` | 1 | 8 | -0.3 (-1.2) | -0.1 (-0.3) | – | – | -0.0 (-0.8) | -0.2 (-0.6) |
| `opt_h_w_1keV_kp` | 1 | 8 | -0.1 (-0.2) | – | – | – | +0.1 (+0.3) | -0.2 (-0.4) |
| `opt_h_w_3000keV_thin10000_kp` | 1 | 8 | – | -11.2 (-1.5) | – | – | -0.0 (-0.3) | -8.7 (-1.6) |
| `opt_he_c_1000keV_kp` | 1 | 8 | +0.0 (+0.0) | +0.3 (+0.5) | – | – | +0.0 (+0.5) | +0.3 (+0.5) |
| `opt_he_fe_1000keV_fc_exyz` | 1 | 8 | -1.2 (-1.0) | +6.4 (+0.6) | +8.3 (+0.6) | – | -0.1 (-0.7) | +5.4 (+0.6) |
| `opt_he_fe_2000keV_thin5000_kp` | 1 | 8 | – | -2.2 (-0.3) | – | – | -0.0 (-0.5) | -1.9 (-0.3) |
| `opt_he_fe_2keV_fc` | 1 | 8 | -1.3 (-1.9) | -0.1 (-0.4) | +1.6 (+0.6) | -0.6 (-0.2) | +0.6 (+1.8) | -0.1 (-0.4) |
| `opt_he_w_2000keV_kp` | 1 | 8 | -0.1 (-1.5) | -0.2 (-0.5) | – | – | -0.0 (-0.2) | -0.1 (-0.5) |
| `opt_o_al2o3_500keV_kp` | 1 | 8 | -0.2 (-1.5) | +0.3 (+0.7) | – | – | -0.1 (-0.6) | +0.3 (+0.7) |
| `opt_si_sic_1000keV_fc` | 1 | 8 | +0.5 (+0.9) | -0.1 (-0.2) | -0.3 (-0.4) | +21.6 (+0.9) | +0.4 (+1.0) | -0.2 (-0.3) |
| `opt_xe_uo2_1000keV_fc` | 1 | 8 | -0.5 (-0.3) | +0.2 (+0.4) | -0.1 (-0.1) | -1.9 (-0.2) | +0.3 (+0.2) | +0.2 (+0.4) |


### 2.3 Ni ions (batch 4, November 2026)

Why: SRIM-2013's stopping of Ni ions is out of line (§3), and §9 offers a
correction; these runs show that AMCSET reproduces TRIM.exe for Ni ions
with the correction off, and how much it changes with it on. Cases:
`tools/srim/gen_cases4.py` (`<SRIM>\_amcset\trim_cases4`, results in
`results\trim4`): Ni into Fe at 1, 3, 5 and 11.5 MeV, into Mo at 5 MeV,
into Ni at 2 and 5 MeV, into Si at 3 MeV, and Co and Cu into Fe at 5 MeV as
controls; quick calculations of 3000-4000 ions, and 5 MeV Ni into Fe also
with full cascades (500 ions); TRIM's defaults (Ed 40 eV, Eb 0, Es from
TRIM's table); one TRIM run each. AMCSET: `batch.py` with 8 seeds, once as
TRIM and once with `--correct-nickel`; `analyze.py`; results in
`build/trim_validation/trim4/` and `trim4_nickel/` (`analysis.txt`,
`fix_effect.txt`). TRIM was run while the workstation had a GPU load; its
run times are not a benchmark.

**As TRIM** (`trim4/analysis.txt`): over the 11 runs, 166 of 169 compared
quantities and 21 of 22 depth profiles lie within 3σ (mean projected ranges
10 of 11 within 2σ, the largest 2.8σ; vacancy profiles median χ²/ndf 1.17,
ion-range profiles 1.10). The four beyond 3σ belong to one run, 5 MeV Ni into Fe (quick): ionization
by the ion +0.70 % (z 3.4), by the recoils −2.6 % (z −3.7), energy to
recoils −1.4 % (z −3.4) and its vacancy profile (χ²/ndf 1.83). They are one
effect seen four ways (energy moved from the recoils to the ion's
ionization), and the same run has the largest range excess, +0.9 % (z 2.8).
The full-cascade run of the same ion and energy leans the same way less
significantly (range +1.1 %, z 1.2; ionization by the ion +1.4 %, z 2.4), as
do Co and Cu into Fe at 5 MeV (range +0.35-0.40 %, ionization by the ion
+0.25-0.35 %, z 1.4-2.3); Ni into Fe at 1 and 3 MeV leans the other way
(range −0.7 and −0.4 %). Whether this ~1 % excess at 5 MeV in Fe is real or
a fluctuation of the one TRIM run per case was not settled (more TRIM seeds
would). It is well below the Ni correction's effect.

**With `nickel_ions`** (`trim4_nickel/fix_effect.txt`; AMCSET with the
correction over TRIM.exe):

| Case | Mean projected range | Range straggling | Vacancies / ion | Ionization by the ion |
|---|---|---|---|---|
| 1 MeV Ni → Fe | 0.973 | 0.975 | 0.978 | 1.040 |
| 3 MeV Ni → Fe | 0.845 | 0.828 | 0.897 | 1.082 |
| 5 MeV Ni → Fe (quick) | 0.789 | 0.734 | 0.848 | 1.083 |
| 5 MeV Ni → Fe (full cascades) | 0.790 | 0.716 | 0.838 | 1.091 |
| 11.5 MeV Ni → Fe | 0.788 | 0.720 | 0.837 | 1.037 |
| 5 MeV Ni → Mo | 0.786 | 0.737 | 0.846 | 1.074 |
| 2 MeV Ni → Ni | 0.926 | 0.888 | 0.940 | 1.082 |
| 5 MeV Ni → Ni | 0.785 | 0.762 | 0.855 | 1.082 |
| 3 MeV Ni → Si | 0.842 | 0.831 | 0.898 | 1.064 |
| 5 MeV Co → Fe, Cu → Fe | 1.004, 1.004 (unchanged by the switch) | 1.000, 1.008 | 0.998, 0.995 | 1.003, 1.004 |

So TRIM.exe itself gives the long Ni ranges that Gurbich et al. measured
against (§4.1): with the correction, ranges are 21 % shorter at 5-11.5 MeV
(where SRIM-2013 is 1.3-1.5 × too long), 7-16 % at 2-3 MeV and 3 % at
1 MeV, the straggling shrinks more than the range, and the vacancies per
ion fall by 15 % (more energy goes to ionization). Co and Cu are untouched.

### 2.4 Against published SRIM calculations

Other people's SRIM runs, as a check that AMCSET's SRIM mode gives what SRIM
users report (`literature_cases.py published_srim`; the papers' Ed, lattice
binding 0, TRIM's Es; 1000 ions × 2, 20 000 for protons). Damage energy here
is the incident energy minus the ionization by the ion and the recoils.

| Source | Case | Mode | Vacancies / ion, paper | AMCSET | Ratio | Damage energy, paper (keV) | AMCSET |
|---|---|---|---|---|---|---|---|
| Stoller et al. 2013, Table 4 (SRIM-2008, Fe, Ed 40 eV) | 69 keV Fe → Fe | quick / full | 469 / 953 | 468.8 / 933.4 | 1.000 / 0.979 | 47.171 / 50.097 | 47.154 / 50.050 |
| | 5 MeV Fe → Fe | quick / full | 7893 / 16 333 | 7879 / 16 261 | 0.998 / 0.996 | 793.255 / 875.652 | 791.849 / 871.874 |
| | 2 MeV H → Fe | quick / full | 14 / 27 | 13.8 / 26.2 | 0.985 / 0.971 | 2.512 / 2.602 | 1.904 / 2.153 |
| Weber and Zhang 2019, Table 1 (full cascades) | 10 MeV Fe → Fe-10Cr | full | 19 197 | 18 303 | 0.953 | | |
| | 50 keV Ni → Ni | full | 877 | 817.4 | 0.932 | | |
| | 5 MeV Ni → Ni | full | 28 148 | 26 915 | 0.956 | | |
| | 3 MeV Au → Ni | full | 33 214 | 33 320 | 1.003 | | |
| | 3 MeV Au → SiC | full | 14 577 | 14 766 | 1.013 | | |
| | 0.9 MeV Au → SrTiO3 | full | 4339 | 4349 | 1.002 | | |

Stoller et al.'s heavy-ion cases agree to 0-2 %, vacancies and damage energy
alike. For the proton their damage energy (the erratum's corrected Eq. 5b,
NIM B 459, 196 (2019)) counts more than the energy not lost to ionization,
so it is not the same quantity; the vacancies agree. Weber and Zhang give
neither their SRIM version, nor their ion count, nor Es; their self-ion
cases are 4-7 % above AMCSET's, their Au cases within 1.3 %. These are
calculations by other users, not measurements; §2.2 is the controlled
comparison with TRIM.exe.

## 3. The SRIM-2013 stopping tables

AMCSET's SRIM mode uses SRIM-2013's electronic stopping, extracted from the
SR Module for every ion in every element (`data/srim2013`). For H and He
ions SRIM's stopping rests on thousands of measurements: Titze et al.
(2021, §4.2) quote agreement to better than 5 % for H and He in all targets
investigated. For heavy ions in most targets there are few or no
measurements, and SRIM scales its proton stopping by an effective charge of
the ion, with ion-specific parameters.

### 3.1 Scan of the tables

`tools/validation/srim_table_scan.py`: for every
ion Z1 = 3-90 in C, Si, Fe, Ni, Cu, Mo, Ag, W and Au, the electronic
stopping at a given energy per mass (0.01-10 MeV/u) is compared with the
geometric mean of its neighbours' (Z1 ± 1, and separately Z1 ± 2, at the
same velocity). Stopping changes smoothly with Z1 apart from the Z1
oscillations at low velocity, which a ratio of both neighbours keeps small.
An ion is listed when both ratios leave 1 ± 0.15 on the same side at one
velocity of 0.1 MeV/u or more:

| Ion | Ratio to Z1 ± 1 | Ratio to Z1 ± 2 | Velocity | Targets |
|---|---|---|---|---|
| Ni (28) | 0.688 | 0.658 | 0.1 MeV/u | all nine |
| Au (79) | 0.753 | 0.750 | 0.1 MeV/u | all nine |
| Li (3) | 1.15-1.23 | 1.95-2.04 | 0.25-1 MeV/u | C, Si, Ni, Cu, Ag, W, Au |

The ratios are the same in every target to three digits: the dip lies in
SRIM's ion-specific scaling, not in a target's data. The Ni dip spans about
0.02-0.3 MeV/u (1-17 MeV for Ni) and is deepest near 0.06 MeV/u (ratio 0.61
to Co/Cu, 0.57 to Mn/Zn, in Fe); the Au dip spans a similar band and is
deepest near 0.08 MeV/u (0.74). Li sits between fully stripped light ions
and heavier ones, where neighbours are not a fair reference (Li stopping is
measured). In iron, at the same energy (bundled tables, `amcset stopping`):

| Ion in Fe | S_e at 2.15 MeV (eV/Å) | S_e at 4.64 MeV (eV/Å) | CSDA path at 10 MeV (µm) |
|---|---|---|---|
| Mn | 277.3 | 494.8 | 2.228 |
| Fe | 248.6 | 450.8 | 2.336 |
| Co | 230.1 | 443.0 | 2.344 |
| **Ni** | **162.6** | **281.7** | **2.971** |
| Cu | 234.3 | 411.7 | 2.364 |
| Zn | 274.6 | 471.6 | 2.082 |

For Ni the measured ranges (§4.1) are too long by about what the dip
predicts, and §9 corrects it. For Au the measured stopping (§3.2) lies above
SRIM's in the dip's band, but the measured ranges of Au in light targets
(§4.3) are already longer than SRIM's: raising Au's electronic stopping
would shorten them further. No correction is offered for Au.

### 3.2 Against the IAEA stopping power database

`tools/validation/stopping_reference.py`, results in
`build/validation/stopping/results.txt` (every point in `iaea_ratios.csv`).
The IAEA's electronic stopping power database (release 2026-09, full
download) holds 274 092 rows with stopping in MeV/(mg/cm²). Of these,
203 025 are SRIM-2013's own table and 4230 are other calculations (MSTAR,
DPASS, CasP, ESPNN, ICRU 49, TDDFT, ...); only the remaining measurements
are compared, in elemental targets (solids, and gases with `--gas`): 38 268
values of 901 ion-target pairs. AMCSET's electronic stopping is interpolated
(log-log, 50 points per decade of `amcset stopping`) at each measured
energy per mass.

**The bundled tables are SRIM-2013's.** Against the database's copy of
SRIM-2013 (120 304 values of 793 pairs) the bundled tables give a median
ratio of 0.9998, with a median deviation of 0.16 %. A few whole pairs of the
copy differ by a constant factor at every energy (O, Si and Mg ions in gases
by 11-17 %, Li and Ag in gases by 5-9 %, Au ions by 5.9 % in every target,
W in Fe by a factor of 270, C and B in Sn and Pb by 2-3), while the bundled
tables of those pairs run smoothly with the neighbouring targets: they look
like unit or phase conversions in the copy, and were not traced further.

**Measurements** (AMCSET with the SRIM tables over measured):

| Energy (MeV/u) | Points | Median | Median \|r − 1\| | Within 10 % | Within 20 % | builtin model: median, within 20 % |
|---|---|---|---|---|---|---|
| all | 38 268 | 0.996 | 0.036 | 81.3 % | 93.9 % | 0.940, 64.8 % |
| below 0.025 | 4923 | 1.000 | 0.109 | 47.3 % | 70.6 % | 0.859, 44.3 % |
| 0.025-0.25 | 12 256 | 0.993 | 0.041 | 79.5 % | 95.3 % | 0.852, 53.6 % |
| 0.25-2.5 | 16 137 | 0.996 | 0.029 | 90.9 % | 98.8 % | 0.971, 70.5 % |
| above 2.5 | 4952 | 0.998 | 0.022 | 88.5 % | 97.8 % | 1.027, 94.5 % |

This agreement is partly circular: SRIM's tables were fitted to many of the
same measurements, above all for H, He and Li (medians 0.97-1.00 in every
band). It says that the tables reproduce the data they rest on, and that
AMCSET's `builtin` stopping (without the tables) is good to a few per cent
only above 2.5 MeV/u. Per ion (the table in `results.txt`), the medians stay
within 10 % except where few points exist, and for heavy ions at 0.025-0.25
MeV/u, where SRIM tends to lie below the measurements: Kr 0.880, Ag 0.854,
Au 0.855, U 0.900, Ni 0.898, Cu 0.903. In the band of SRIM's Au dip (§3.1),
Au measured by Jokinen (ToF-ERDA, 1997) gives 0.82 in Al, 0.93 in C and 0.74
in Cu at 57-187 keV/u, and 0.88 in Ni (2008); Pb, three elements on, lies
within 8 % in eight targets (Geissel 1982).

**Ni ions** (the measurements that test §9's correction), by data set:

| Data set | Target, energy | Points | SRIM tables | With `nickel_ions` | Neighbours in the same data set |
|---|---|---|---|---|---|
| Zhang 2002 (ToF-ERDA; NIM B 196, 1) | C, 80-246 keV/u | 14 | 0.815 | 0.966 | Fe 0.966, Co 1.008, Cu 0.927 |
| | Al, 70-242 keV/u | 14 | 0.848 | 1.029 | Fe 1.055, Co 1.088, Cu 0.968 |
| Harikumar 1997 | C, 140-210 keV/u | 2 | 0.673 | 0.792 | Fe 0.893, Cu 0.865 |
| Anthony and Lanford 1981 | Cu, 33-228 keV/u | 14 | 1.138 | 1.298 | — |
| | Ag, 39-215 keV/u | 12 | 0.970 | 1.161 | — |
| all Ni measurements | 0.025-0.25 MeV/u | 56 | 0.898 | 1.042 | |
| | 0.25-2.5 MeV/u | 225 | 0.963 | 0.990 | |
| | above 2.5 MeV/u | 32 | 1.014 | 1.041 | |

The one modern data set that measured Ni with its neighbours in the same
apparatus (Zhang 2002) finds SRIM's Ni 15-18 % too low where Fe, Co and Cu
are within 9 %, and the correction brings Ni to its neighbours' level. The
older data of Anthony and Lanford lie at or below SRIM's Ni: they are
consistent with the dip, and SRIM's Ni may rest on them. The measured
stopping cannot settle this alone; the measured ranges (§4.1) can, and side
with the correction. No measurement of Ni-ion stopping in Fe or Ni, where
the ranges were measured, is in the database.

### 3.3 Against NIST PSTAR and ASTAR

Protons and alpha particles (ICRU Report 49) in Si, Fe, Cu, W, Au and water,
at NIST's default energies; median ratio of AMCSET's electronic stopping to
NIST's per band, in brackets the ratio farthest from 1 in the band:

| Target | H 0.01-0.1 MeV | H 0.1-1 MeV | H 1-10 MeV | H 10-1000 MeV | He 0.01-0.1 MeV/u | He 0.1-1 MeV/u | He 1-1000 MeV/u |
|---|---|---|---|---|---|---|---|
| Si | 0.980 (0.943) | 0.996 (0.983) | 1.007 | 1.003-1.004 | 0.993 (1.020) | 1.014 (0.970) | 1.004-1.011 |
| Fe | 1.120 (1.149) | 1.027 (1.085) | 1.008 | 1.002-1.004 | 1.021 (1.062) | 0.995 (1.034) | 1.004-1.008 |
| Cu | 1.016 (1.037) | 1.003 (1.041) | 1.000 | 0.998 | 1.005 (0.967) | 0.993 (0.984) | 0.996-1.002 |
| W | 1.079 (1.110) | 1.019 (0.958) | 1.030 | 1.000-1.007 | 1.025 (1.050) | 1.002 (1.034) | 1.002-1.030 |
| Au | 0.981 (0.975) | 1.013 (1.026) | 1.004 | 1.001-1.002 | 0.950 (0.941) | 0.971 (0.952) | 1.003-1.007 |
| water (H2O, 1 g/cm³) | 1.100 (1.113) | 0.960 (1.078) | 1.026 | 1.010-1.014 | 1.013 (1.048) | 1.019 (1.046) | 1.012-1.033 |

Above 0.1 MeV (per nucleon for He) the medians of the SRIM tables lie
within 3 % of NIST in the elements (single energies up to 8.5 %, H in Fe at
0.1-1 MeV) and within 4 % in water; below, where ICRU 49 rests on fits of
its own, up to 12 % (H in Fe). Water is a compound by Bragg's rule here:
SRIM's tables carry no phase or bond correction. CSDA ranges (srim) lie
within 4 % of NIST's above 1 MeV (per nucleon), up to 18 % off below (He in
Fe at 0.01-0.1 MeV/u). Without the tables, `builtin` is 0.5-6 % high at
10-100 MeV (up to 11 % at single energies, W) and up to 60 % off below 1 MeV
(W, Au); `bethe-bloch` (only above 1 MeV/u) is 0.3-6 % high at 10-100 MeV
and up to 14 % at 1-10 MeV (W).

## 4. Ranges against measurements

### 4.1 Ni ions in Ni, Fe, steels and Mo

**Measured peak depths** (`literature_cases.py ni_peaks`). The papers:

- Chauhan, Pena, Jiang, Hu, Garner, Mejia Diaz, Ullah, Correa, Zhu, Kao and
  Shao, J. Nucl. Mater. 633, 157044 (2026), doi 10.1016/j.jnucmat.2026.157044
  (open access): ⁶²Ni at 2.0, 3.5 and 5.0 MeV into pure Ni and pure Fe,
  normal incidence, room temperature, 1.5-4.65·10¹⁷ cm⁻²; ToF-SIMS, peak of
  a bi-Gaussian fit; their SRIM-2013 quick calculation with Ed 40 eV. At
  5 MeV their SRIM peak lies 625 nm (about 36 %) beyond the measured one in
  Ni and 501 nm (26 %) in Fe.
- Pena, Wang, Jiang, Chauhan, Garner, Nastasi and Shao, Nucl. Instrum.
  Methods B 565, 165738 (2025): ⁶²Ni into Ni at 2.0 and 3.5 MeV,
  4.6·10¹⁷ cm⁻², ToF-SIMS peaks 364 and 1008 nm (SRIM 720 and 1290 nm).
- Song, Was and Garner, J. Nucl. Mater. 633, 157002 (2026), Table 2 (read in
  its SSRN preprint): the pure-Fe peaks of Chauhan et al.; 8 MeV ⁵⁸Ni in Fe
  and Fe-10Cr (Lin et al., EDX); 11.5 MeV in Fe (Gurbich et al. 2025,
  resonance profiling, 1.92 ± 0.16 µm); 8 MeV into 316 and 9 MeV into 304
  stainless steel (STEM-EDX), with their SRIM-2013 and TRIM-98 peaks.

AMCSET: quick calculation, Ed 40 eV (Mo 60 eV), 20 000 ions × 4 seeds,
targets 2.2 × the papers' SRIM peak deep (bins of 2.2 % of it); peak of the
stopped-ion profile by a parabola through the logarithm of the smoothed
histogram above 60 % of its maximum. Steels with Song et al.'s main
elements (at.%) and handbook densities (316: 7.99, 304: 7.93 g/cm³). The
measured peaks are not corrected for sputtering during the implantation
(Chauhan et al. say so); the last columns estimate it: Y from full-cascade
runs (1000 ions × 2, with the correction), d = Y·fluence/N the depth
sputtered away, and the peak of the corrected profile seen from the final
surface (the profile averaged over shifts from 0 to d).

| Case (source) | Measured (µm) | Paper's SRIM | AMCSET, SRIM tables | / paper's SRIM | With `nickel_ions` | / measured | Y | d (nm) | Sputtered, / measured |
|---|---|---|---|---|---|---|---|---|---|
| 2.0 MeV → Ni (P) | 0.364 ± 0.05 | 0.720 | 0.694 | 0.964 | 0.646 | 1.776 | 2.14 | 108 | 1.622 |
| 3.5 MeV → Ni (P) | 1.008 ± 0.05 | 1.290 | 1.256 | 0.973 | 1.044 | 1.036 | 1.42 | 72 | 0.998 |
| 5.0 MeV → Ni (C) | 1.111 ± 0.05 | 1.736 | 1.731 | 0.997 | 1.356 | 1.220 | 1.37 | 23 | 1.212 |
| 2.0 MeV → Fe (C) | 0.584 ± 0.010 | 0.785 | 0.753 | 0.959 | 0.698 | 1.195 | 1.49 | 41 | 1.159 |
| 3.5 MeV → Fe (C) | 1.137 ± 0.010 | 1.380 | 1.361 | 0.986 | 1.123 | 0.988 | 1.12 | 62 | 0.960 |
| 5.0 MeV → Fe (C) | 1.389 ± 0.010 | 1.890 | 1.870 | 0.989 | 1.456 | 1.048 | 0.90 | 49 | 1.030 |
| 8 MeV → Fe, Fe-10Cr (Lin, in S) | 1.750 | 2.50 | 2.492 | 0.997 | 1.904 | 1.088 | 0.55 | — | — |
| 11.5 MeV → Fe (Gurbich, in S) | 1.92 ± 0.16 | 2.80 | 2.988 | 1.067 | 2.326 | 1.211 | 0.52 | — | — |
| 8 MeV → 316 SS (S) | 1.815 ± 0.003 | 2.44 | 2.395 | 0.982 | 1.827 | 1.006 | 0.64 | 9 | 1.004 |
| 9 MeV → 304 SS (S) | 1.925 ± 0.030 | 2.62 | 2.543 | 0.971 | 1.952 | 1.014 | 0.59 | 9 | 1.012 |

(C: Chauhan et al.; P: Pena et al.; S: Song et al. The 5 MeV measurement in
Ni is the paper's SRIM peak, read off its figure, minus the 625 nm of its
text; it had the lowest fluence and the noisiest profile. Standard errors of
AMCSET's peaks are 1-2 nm.)

What it shows:

- AMCSET with the SRIM-2013 tables reproduces the papers' own SRIM peaks to
  -4 to +7 % (their bins were 40-50 nm), and lies 1.2-1.9 × deeper than the
  measured peaks, as SRIM-2013 does.
- With `nickel_ions`, six of the ten measured peaks (3.5-9 MeV, in Ni, Fe and
  both steels) are reproduced within 9 %, three within 2 % after the
  sputtering estimate. 5 MeV in Ni and 11.5 MeV in Fe stay 21-22 % too deep
  (the latter 2.5 σ of its measurement), and 2 MeV 1.2 × (Fe) and 1.8 × (Ni)
  too deep.
- The sputtering estimate moves the peaks by 0.2-9 %; TRIM-type sputtering
  yields can be off by a factor up to 2 (§4.5), so the 2 MeV peaks,
  implanted to the highest dose per depth, are the least certain. At 2 MeV
  the correction itself is small: a 2 MeV Ni ion (0.034 MeV/u) starts near
  the low edge of the dip of §3.1 and spends most of its path below it (the
  correction shortens 2 MeV ranges by 7 %, §2.3). Song et al.'s TRIM-98
  peaks in Fe (0.615, 0.962, 1.255 µm at 2, 3.5, 5 MeV) lie 5 % beyond the
  measured peak at 2 MeV and 10-15 % short of it at 3.5-5 MeV; the
  correction is the other way round. Neither fits every energy.
- Chauhan et al. also find that RT-TDDFT stopping (with the projectile's
  semicore electrons) with MD ion transport reproduces their Ni-in-Ni
  peaks; Song et al. and Chauhan et al. recommend TRIM-98's stopping, or an
  energy-dependent correction rather than a single factor (Lin et al.'s
  1.5 at 8 MeV).

**Mean projected ranges** (`literature_cases.py ni_ranges`), from Gurbich,
Prusachenko, Bobrovskiy and Bokhovko, "On the SRIM prediction of ranges for
Ni ions", Appl. Radiat. Isot. 215, 111590 (2025), whose full text was not
available; its abstract: SRIM-2013 overestimates the projected range of
5 MeV Ni in Mo by 1.3 times and of 11.5 MeV Ni in Fe by 1.5 times; TRIM-98,
with radiation-enhanced diffusion taken into account, agrees with the data;
a SRIM-2013 calculation for Co was close to TRIM-98's for Ni, which points
at SRIM-2013's Ni stopping. AMCSET (quick calculation, 2000 ions × 4 seeds,
8 µm targets, TRIM's default Ed/Eb/Es), mean projected range:

| Case | SRIM-2013 tables | `nickel_ions` (§9) | Ratio | Paper: SRIM-2013 / measured |
|---|---|---|---|---|
| 5 MeV Ni → Mo | 1.493 ± 0.005 µm | 1.171 ± 0.004 µm | 1.275 | 1.3 |
| 11.5 MeV Ni → Fe | 2.852 ± 0.003 µm | 2.240 ± 0.002 µm | 1.273 | 1.5 |

The correction accounts for the Mo discrepancy (1.275 of 1.3) and for 1.27
of the Fe case's 1.5. The rest of the Fe case is not explained here: the
paper's agreement with TRIM-98 needed radiation-enhanced diffusion of the
implanted Ni, which no transport code includes.

### 4.2 Sb in Si (SIMS)

Titze, Pacheco, Byers, Van Deusen, Perry, Weathers and Bielejec, J. Vac.
Sci. Technol. A 39, 063222 (2021), doi 10.1116/6.0001406 (OSTI 1830531):
Sb implanted into Si with a 2.9 ± 1.2 nm native oxide at 7°, depth by SIMS
and RBS, fluences below saturation. At 50 keV the SIMS depth is 27.5 ± 1.3
nm with a straggle of 9.7 nm; their SRIM-2013 run (full cascades, 99 999
ions; SiO2 2.32 g/cm³, Si 2.329 g/cm³; Ed/Eb/Es O 28/3/2 eV, Si 15/2/4.7 eV)
gives 32.8 nm. They find SRIM's range too long by 2-6 nm from 10 to 150 keV
(up to 25 % at 10 keV) and the straggle right, and attribute it to SRIM's
low-energy stopping of heavy ions (only 5 Sb data points in SRIM).

AMCSET with the same target and energies (`literature_cases.py sb_si`, full
cascades, 2000 ions × 4 seeds): mean depth 33.53 ± 0.16 nm, straggling
9.76 nm. That is SRIM-2013's result (+2.2 % over their 32.8 nm, which is the
centre of a Gaussian fit rather than the mean) and 22 % deeper than measured.
The other energies are only in the paper's figures and are not used.

### 4.3 Ranges in Si, SiO2, C, B and SiC

`literature_cases.py ranges`, data and sources in
`tools/validation/range_data.py`, every point in
`build/validation/literature/ranges.csv`. Quick calculation, bundled
SRIM-2013 tables, 4000 ions × 2 seeds per point (peaks: 20 000 × 2), the
papers' nominal densities, amorphous targets at normal incidence unless
stated. Ratio AMCSET / measured of the mean projected range Rp (Zhang et
al.: of the peak):

| Data | Points | Rp: median | mean ± sd | Within 10 % | Straggling: median | Median per ion |
|---|---|---|---|---|---|---|
| Behar et al. 1995, Tables 3-4: Cu-Bi in amorphous Si, 10-390 keV (RBS) | 93 | 1.025 | 1.019 ± 0.124 | 66 % | 0.992 | Cu 1.03, Ga 1.04, Br 1.08, Rb 1.14, Pd 1.02, Sn 1.14, Cs 1.22, Eu 0.95, Yb 0.96, Au 0.86, Bi 1.08 |
| Behar, Table 5: As-Bi in SiO2, 10-400 keV (RBS) | 51 | 1.063 | 1.065 ± 0.100 | 63 % | 0.842 | Bi 1.11, Au 1.01, Yb 1.00, Eu 0.94, Cs 1.13, Xe 1.13, As 1.05 |
| Behar, Tables 6-9: H, Li, B in amorphous Si, 4.5 keV-2 MeV (nuclear reactions) | 28 | 0.990 | 0.996 ± 0.052 | 96 % | 0.912 | H 1.01, Li 1.02, B 0.96 |
| Oosterhoff 1988, Table 1: B, P in (100) Si at 7°, 0.1-1.5 MeV (SIMS; amorphous here) | 21 | 1.055 | 1.055 ± 0.043 | 95 % | 0.924 | B 1.07, P 1.06 |
| Behar, Tables 11-13: K-Bi in C films, 10-300 keV (RBS) | 88 | 0.860 | 0.876 ± 0.095 | 33 % | 0.735 | Bi 0.88, Pb 0.76, Au 0.77, Yb 0.84, Er 0.79, Eu 0.90, Sn 0.93, Ga 0.81, Cu 0.86, Cs 0.96, Xe 0.94, Rb 0.94, Kr 0.95, Ar 0.87, K 0.90 |
| Behar, Table 15: Rb-Bi in B films, 20-300 keV (RBS) | 24 | 0.851 | 0.863 ± 0.090 | 33 % | 0.603 | Pb 0.85, Bi 0.86, Au 0.75, Yb 0.82, Cs 0.99, Rb 0.98 |
| Behar, Table 16: F-Au in SiC films, 20-300 keV (RBS; F by a resonance) | 29 | 0.953 | 0.962 ± 0.126 | 45 % | 0.839 | Au 0.83, Er 0.87, Cs 1.03, Ga 0.98, F 1.11 |
| Zhang et al. 2009, Tables I-III: peak depth of Pb and Au in SiC, 20 keV-10 MeV (RBS, SIMS) | 16 | 0.802 | 0.779 ± 0.069 | 0 % | — | Pb 0.84, Au 0.78 |

The sources: Behar, Fichtner, Grande and Zawislak, Mater. Sci. Eng. R 15,
1-83 (1995) (Si amorphized by Ar or Si beams; C, B and SiC films 200-300 nm;
values as printed, with two of its misprints left out: Yb 20 keV in SiO2 and
Cs 250 keV in SiC); Oosterhoff, Nucl. Instrum. Methods B 30, 1 (1988)
(as-implanted SIMS columns, 2·10¹⁴-10¹⁵ cm⁻², 7° tilt and about 27° twist
to avoid channeling); Zhang, Bae, Sun et al., J. Appl. Phys. 105, 104901
(2009) (20-250 keV from Fichtner et al. 1992; 0.7-10 MeV their SIMS, at 5°,
2 MeV at 60°).

What it shows:

- **Si and SiO2: agreement.** Medians 0.99-1.06, as Behar et al. found for
  ZBL; the light ions (H, Li, B, P) within 1-7 % per ion. Single heavy ions
  deviate by up to 22 % (Cs in Si +22 %, Au −14 %), the "Z1 oscillation" of
  low-energy ranges, which Behar et al. remove for Au, Yb and Eu below
  70 keV with a cut-off interatomic potential and Wittmaack (§4.4)
  attributes to transport of the implanted atoms. The straggling agrees in
  Si (0.99) and is 8-16 % small in SiO2 and for the light ions.
- **Heavy ions in light targets: AMCSET, like SRIM, is too short.** For Pb,
  Bi, Au, Yb and Er in C, B and SiC, AMCSET's ranges are 12-25 % short of
  the measured ones, and the straggling 16-40 % smaller than measured
  (group medians); Cs, Xe, Rb and Kr in C, Cs and Rb in B, and F, Ga and Cs
  in SiC agree within 2-11 %. Behar et al. remove the deficit with
  a correlation of electronic and nuclear energy loss in close collisions
  (one parameter, α = 0.5, a modified TRIM); Zhang et al. instead find
  SRIM's electronic stopping of slow heavy ions in light compounds about
  twice too high (their reciprocity estimate is 25-47 % of SRIM's).
  Weber and Zhang (2019) reproduce the measured Au profile in SiC by
  lowering the density to 2.78 g/cm³ in SRIM. Neither remedy is in AMCSET.
- AMCSET reproduces Zhang et al.'s SRIM-2008 peaks to 1 % up to 250 keV
  and 2-4 % at 0.7-10 MeV.

### 4.4 What the literature says about SRIM ranges

- Medium-heavy ions in light targets (C, B, Be, SiC, polymers), 10 keV-7.5
  MeV: measured ranges exceed ZBL's by up to 40 %; in Si, ions with
  29 ≤ Z1 ≤ 83 at 10-390 keV agree to better than 10 %, except Au, Yb and Eu
  below about 70 keV (Behar et al. 1995, §4.3).
- Slow heavy ions in light targets: SRIM's electronic stopping is too high
  (Zhang et al. 2009); Weber and Zhang (2019) recommend full cascades for
  damage in compounds, stating density, Ed and binding energies, and
  checking ranges against measurements.
- Heavy ions at 0.1-50 keV/u: Wittmaack (Nucl. Instrum. Methods B 380
  (2016)) finds SRIM's heavy-ion electronic stopping up to 4-5 × the
  Lindhard-Scharff value and Pt, Au, Pb ranges 15-33 % longer than SRIM's;
  his second paper (NIM B 388, 15 (2016)) re-adjusts four groups' RBS
  ranges in Si below 0.4 keV/u by up to ±10 % and finds the "Z1
  oscillations" of those ranges one-sided: implanted noble metals and
  lanthanides moved deeper by the bombardment, rare gases and alkalis
  towards the surface. That matches §4.3's Si deviations (Au 14 % deeper
  than calculated, Cs 22 % shallower); a binary-collision code without
  transport of the implanted atoms cannot follow either.
- Ni ions: §4.1. Fe ions are reported too deep by about 18 % at 3.5 MeV in
  HT9 and too shallow by 18 % at 5 MeV in Cr (as quoted by Chauhan et al.
  2026), so range errors are not Ni's alone, but only Ni's stopping is out
  of line with its neighbours (§3.1).

AMCSET's SRIM mode inherits all of these with SRIM's tables: a statement
about SRIM's ranges is a statement about AMCSET's.

### 4.5 Sputtering yields

`literature_cases.py sputtering`, data in `tools/validation/sputter_data.py`,
every point in `build/validation/literature/sputtering.csv`. Measured yields
at normal incidence read off the figures of Yamamura and Tawara,
NIFS-DATA-23 (1995; = At. Data Nucl. Data Tables 62, 149 (1996)), for Ar on
Si, Fe, Cu, Ag, W and Au (0.03-500 keV), and of Yamamura, Sakaoka and
Tawara, NIFS-DATA-31 (1996), for H, D and He on Fe and W (0.06-100 keV):
only the experimental points, read to about 20-30 % (more near the
threshold); where data sets disagree the main cluster was read; two points
far below the others are left out. AMCSET: full cascades, TRIM's default
Ed, Eb (3 eV) and Es, a target sized by `--auto-depth`, enough ions for about
2000 sputtered atoms (at most 200 000).

| Ion on target | Points | AMCSET / measured: median | Above 0.5 keV: median | Above 0.5 keV: range | No yield at |
|---|---|---|---|---|---|
| Ar on Si | 18 | 0.64 | 0.80 | 0.54-0.89 | — |
| Ar on Fe | 15 | 1.50 | 1.45 | 1.20-1.64 | — |
| Ar on Cu | 17 | 0.97 | 0.94 | 0.85-1.14 | — |
| Ar on Ag | 15 | 0.85 | 0.85 | 0.74-0.88 | 0.03 keV |
| Ar on W | 15 | 1.26 | 1.32 | 1.13-1.48 | — |
| Ar on Au | 15 | 0.79 | 0.68 | 0.53-0.80 | — |
| H on Fe | 10 | 0.67 | 0.92 | 0.66-1.20 | ≤ 0.3 keV |
| D on Fe | 9 | 1.03 | 1.05 | 0.84-1.15 | ≤ 0.1 keV |
| He on Fe | 10 | 1.33 | 1.43 | 1.31-1.82 | ≤ 0.1 keV |
| H on W | 10 | 1.72 | 1.72 | 0-2.40 | ≤ 1 keV |
| D on W | 12 | 1.60 | 1.72 | 0-2.34 | ≤ 0.5 keV |
| He on W | 12 | 1.81 | 1.88 | 1.76-3.27 | ≤ 0.25 keV |
| all | 158 | 1.04 | 1.05 (98 points, 88 % within a factor 2) | | |

The median over all pairs is 1.05, but each pair has its own factor, from
0.68 (Ar on Au) to 1.88 (He on W), steady over one to three decades of
energy: the yields follow the surface binding energy and the cascade near
the surface that TRIM's model gives, not the measurements' scatter. Near the
threshold AMCSET (as TRIM) sputters nothing where small yields are measured
(H on W below 1 keV, He on W below 0.3 keV, H on Fe below 0.5 keV): TRIM's
cascade ends at Ed and Eb, and its surface model needs more energy. Ar on
Si near 50-200 eV is 3-8 × low. This is SRIM's sputtering model; the
comparison with TRIM.exe (§2.2, 36 sputtering cases) shows that AMCSET
reproduces it. No correction is offered: Es per element can be set
(`layers[].elements[].Es`), which is how users usually adjust TRIM's
yields.

## 5. Damage models against MD, experiments and SRIM post-processing

Not a TRIM comparison: these runs test the `dpa` table (docs/SRIM_MODE.md §2.8) and ask which
damage model describes the damage best. Scripts: `tools/damage_models/`;
results on the workstation in `build/damage_models/` (`eb0/`, `fe_eb3/`,
`mitsi2023/`).

**1. The post-processing is right.** Mitsi, Koutsomitis and Apostolopoulos
(Nucl. Instrum. Methods B, 2023; arXiv 2307.12867) applied arc-dpa to every
primary in COLLISON.txt of SRIM-2013 quick calculations (their method M1,
10 000 ions) for 1 MeV H, 1 MeV He, 3 MeV Al, 5 MeV Fe and 10 MeV Au into
Fe, Ni, Cu, Pd, W and Pt (their Ed: 40, 40, 29, 41, 90, 44 eV). AMCSET's
quick calculation with the same Ed (2 × 2000 ions per case,
`arc_dpa_vs_srim_postprocessing.py`) gives their arc/NRT ratio in all 30
cases within 0.43 % rms (largest −0.92 %, for 1 MeV H into Pd, whose
reference ratio has three digits and AMCSET's a standard error of 2 %).
Per ion, Fe into Fe: NRT 529.5 and arc-dpa 216.6 at 78.7 keV (theirs: 530,
217), 7867 and 2891 at 5 MeV (7870, 2900). In every case `nrt` equals TRIM's
quick-calculation vacancies.

**2. Damage functions.** PKAs of one energy started in the middle of 20 µm
pure targets, isotropic directions, full cascades, Eb = 0, the Ed of the
arc-dpa fits, 400 PKAs (fewer above 10 keV, at least 50) × 4 seeds per point
(`pka_damage_functions.py`). Per PKA, as ratios to NRT:

| Target (Ed) | PKA energy | TRIM full-cascade vacancies | arc-dpa | TRIM replacements | rpa |
| --- | --- | --- | --- | --- | --- |
| Fe (40 eV) | 1 / 10 / 100 keV | 2.00 / 2.01 / 2.10 | 0.502 / 0.348 / 0.304 | 0.107 / 0.108 / 0.113 | 4.5 / 8.7 / 9.9 |
| Cu (33 eV) | 1 / 10 / 100 keV | 2.24 / 2.31 / 2.54 | 0.336 / 0.199 / 0.169 | 0.113 / 0.116 / 0.128 | 7.6 / 25.5 / 35.0 |
| Ni (39 eV) | 1 / 10 / 100 keV | 2.22 / 2.30 / 2.53 | 0.320 / 0.240 / 0.231 | 0.097 / 0.099 / 0.108 | 5.8 / 18.1 / 25.1 |
| Pd (41 eV) | 1 / 10 / 100 keV | 2.11 / 2.12 / 2.22 | 0.283 / 0.169 / 0.153 | 0.130 / 0.131 / 0.138 | 7.3 / 21.5 / 26.0 |
| Pt (42 eV) | 1 / 10 / 100 keV | 2.16 / 2.11 / 2.12 | 0.194 / 0.117 / 0.111 | 0.250 / 0.241 / 0.242 | 5.4 / 18.9 / 29.4 |
| W (70 eV) | 1 / 10 / 100 keV | 1.97 / 1.89 / 1.91 | 0.481 / 0.223 / 0.150 | 0.142 / 0.136 / 0.139 | 2.9 / 9.9 / 18.3 |

Kinchin-Pease is 1.25 × NRT above 2.5 Ed by definition. The standard
errors of the full-cascade vacancies are at most about 0.5 %. With
Eb = 3 eV (Fe), the full-cascade vacancies fall by 4-5 % (1.92 / 1.92 /
2.01 × NRT) and the replacements by 11-12 %. The whole table, 0.1-200 keV:
`build/damage_models/eb0/damage_functions.txt`.

**3. Against MD.** arc-dpa and rpa are fits to MD cascades (Wigner-Seitz
defect analysis, ZBL-joined potentials, ZBL96 electronic friction); their
authors report that both agree "very well with MD simulation results over
more than three orders of magnitude in energy". So, against MD:

- TRIM's full-cascade vacancies overestimate the surviving Frenkel pairs by
  4× (Fe, W at 1 keV) to 19× (Pt at 100 keV). They count every atom moved
  above Ed in a binary-collision cascade without recombination.
- NRT overestimates them by 2-9× (1/arc-dpa in the table), Kinchin-Pease
  by another 25 %.
- TRIM's replacements are 21-275× fewer than the replaced atoms of MD
  (rpa/replacements); a binary-collision cascade has no heat spike, which
  is where most mixing happens.
- An independent, newer MD study of W (Byggmästar et al., arXiv 2503.02710,
  2025, machine-learned potential, PKAs up to 2 MeV) finds arc-dpa (with a
  low-energy correction and Ed = 95 eV) to describe the sublinear regime up
  to about 20-30 keV, then a superlinear regime (overlapping subcascades,
  large clusters) and a linear one above about 300 keV that "deviates from
  all previous models". arc-dpa may therefore underestimate the primary
  damage of the highest-energy PKAs in W.

**4. Against experiments.** Averback, Benedek and Merkle, Phys. Rev. B 18,
4156 (1978), measured the resistivity increase of Cu and Ag foils
irradiated at 6 K with several ions, relative to a modified Kinchin-Pease
damage function: efficiency about 1 for H, falling with the projectile
mass, and about 0.35 (Cu) and 0.4 (Ag) for ions heavier than Ne, constant
above recoil energies of a few keV. Nordlund et al. (2018) summarise the
experimental and MD data for Cu as about 1/3 of NRT for energetic cascades.
For Cu, AMCSET's arc/NRT is 0.63 for 1 MeV H, 0.53 for 1 MeV He and 0.19-0.26
for 3-10 MeV Al, Fe and Au (Ed = 29 eV, §5 item 1); per PKA 0.20 at
10 keV and 0.17 at 100 keV (Ed = 33 eV). TRIM's full cascades give 2.2-2.5×
NRT. (The paper's ions, energies and Ed are behind a paywall and were not
reproduced here; the comparison is of efficiencies, not of the same runs.)

**Conclusion.**

- **Surviving defects in metals (primary damage)**: arc-dpa is the most
  accurate of the four, by construction against MD, and the closest to the
  6 K experiments in Cu for cascade-producing ions (0.19-0.26 against about
  0.35, where NRT gives 1 and TRIM's full cascades 2.2-2.5). Its Cu
  constants put it about 25-45 % below those experiments, for light and
  heavy ions alike. No damage model is accurate for light ions at the
  ~10 % level.
- **Atomic mixing**: rpa; TRIM's replacements are not a measure of it.
- **Exposure / dose for comparing irradiations**: NRT-dpa remains the
  standard (ASTM E521). Its authors recommend it for comparing irradiations;
  Stoller et al., Nucl. Instrum. Methods B 310, 75 (2013), recommend TRIM's
  quick calculation (= AMCSET's `nrt`) when ion results are to be compared
  with neutron dpa.
- **TRIM's full-cascade vacancies**: neither a dose standard nor a defect
  count; use them for where atoms are displaced, not for how many survive.
- None of them includes thermally activated recombination, clustering or
  dose overlap.

## 6. The crystal transport against measurements

The crystal transport (`physics.transport: crystal`) is AMCSET's own, not
TRIM's. Its checks without measurements (random direction = amorphous,
symmetry, the normalisation of the stopping division, the cut-off at Ed
against NRT) are in docs/CLI.md "The crystal transport"; why ranges in a
crystal exceed the amorphous ones, with seven ion-crystal pairs, is
docs/CLI.md "Why ranges in a crystal are longer" (`range_study.py`).

**Against measurements** (`tools/crystal/calibration.py`; SRIM-2013
stopping; a perfect crystal for every ion; the numbers are in
`build/crystal_calibration/results.txt` on the workstation).

- *Thermal amplitudes.* The Debye model reproduces the tabulated u₁(293 K)
  of all 19 elements of Gemmell's Table IV from their Debye temperatures
  (test `CrystalMaterials`).
- *Channeling dip: 1 MeV He along ⟨110⟩ of Ge at 293 K* (atoms displaced by
  the ion at depths of 200-1200 Å, relative to amorphous Ge). Minimum yield
  0.022 and half-angle 0.99° = 0.80 ψ₁, against the measured 0.03 and 0.95°
  (Picraux et al., as quoted by Gemmell) and Barrett's formulas, 0.025 and
  0.93°. The shoulder reaches 1.56-1.61 at 1.7-2.0°.
- *Boron in (100) Si, as-implanted SIMS profiles*, for five divisions of
  the stopping:

| Measurement | Measured | Crystal transport, fraction 0 / 0.25 / 0.5 / 0.75 / 1 | `trim2013` (amorphous) |
|---|---|---|---|
| 15 keV, 7° tilt, 3·10¹⁵ cm⁻² (Jung 2010, Table 2): Rp (nm) | 64.2 | 62.4 / 62.7 / 62.5 / 63.0 / 62.8 | 60.5 |
| ... ΔRp (nm) | 31.8 | 32.3 / 33.1 / 33.5 / 34.3 / 34.7 | 25.7 |
| ... skewness | 0.98 | 0.59 / 0.69 / 0.77 / 0.89 / 1.00 | 0.00 |
| ... kurtosis | 5.56 | 3.45 / 3.78 / 4.05 / 4.56 / 5.15 | 2.52 |
| 30 keV, 1.5·10¹⁴ cm⁻² (Shao 2003, Fig. 33): depth of 10¹⁸ cm⁻³ (nm) | about 210 | 225 / 223 / 223 / 223 / 225 | 201 |
| ... depth of 10¹⁷ cm⁻³ (nm) | about 290 | 285 / 289 / 306 / 315 / 319 | 226 |
| 0.5 keV, 10¹⁴ cm⁻² (Schroer 2000, Fig. 1): depth of 10¹⁹ cm⁻³ (nm) | 14.2 | 14.6 / 14.8 / 15.0 / 15.2 / 15.4 | 9.7 |
| ... depth of 10¹⁸ cm⁻³ (nm) | 25.2 | 19.1 / 19.9 / 20.6 / 21.4 / 22.1 | 12.0 |
| ... depth of 10¹⁷ cm⁻³ (nm) | 35 | 21.2 / 22.5 / 24.2 / 26.4 / 28.7 | 13.8 |
| ... depth of 10¹⁶ cm⁻³ (nm) | 57 | 22.7 / 24.4 / 27.4 / 31.1 / 37.8 | 15.1 |

  The papers: Jung, Trans. Electr. Electron. Mater. 11, 120 (2010); Shao et
  al., Mater. Sci. Eng. R 42, 65 (2003); Schroer et al., Mater. Sci. Eng. B
  71, 219 (2000). The 15 keV moments are those of the profile to 450 nm, as
  in the paper, which gives 59.7 nm, 30.1 nm, 0.72 and 4.62 for UT-Marlowe;
  it does not state the twist (22° is assumed; 45° gives Rp 65.5 nm, and 0°,
  a planar channel, 91.6 nm). The 30 keV and 0.5 keV values were read off
  the published figures (a few per cent of the depth), which state neither
  tilt nor twist: 7°/22° is assumed at 30 keV and 0° at 0.5 keV (7°/22°
  gives 10.7 nm for the first 0.5 keV level; 10.1 nm was measured at
  10¹⁵ cm⁻²).

What this shows. The peak and the width are reproduced whatever the
division: Rp to 2-3 %, ΔRp to 2-9 %, the first decade of the 0.5 keV profile
to 3-8 %, where the amorphous calculation misses the width by 19 % and the
0.5 keV depth by 32 %. The division itself is not fixed by these data: the
30 keV tail prefers about 0.25; the skewness and kurtosis at 15 keV and the
deep 0.5 keV tail prefer 1. The unknowns of the measurements are as large as
its effect: tilt and twist; the damage of the earlier ions at 10¹⁵ cm⁻² and
above, which dechannels and which the perfect crystal of the simulation
lacks; and the deep SIMS tail at 0.5 keV, measured to 57 nm where no
division gives more than 38 nm. The default therefore stays 0.5. Fitting it
needs low-dose implants with stated tilt and twist.

**Channeling tails of As in Si against tilt** (`literature_cases.py
as_si_channeling`). Tasch, Yang, Morris and Lim, J. Vac. Sci. Technol. B 12,
166 (1994): 15 keV As into (100) Si, 10¹⁴ cm⁻², rotation 0, SIMS, tilt 0-10°,
and into amorphous Si; wafers with about 1 nm of native oxide, beam
divergence below 0.5°. Values read off their Fig. 3 (log scale, about a
factor 2): peak 6·10¹⁹ cm⁻³ at about 12 nm; at 100 and 200 nm 3·10¹⁷ and
3·10¹⁶ cm⁻³ at 0°, 10¹⁶ and 2·10¹⁵ at 10°, below 10¹⁵ in amorphous Si.
AMCSET: 200 000 ions, the crystal transport with its defaults (perfect
crystal, 300 K, stopping division 0.5; the reference direction ⟨110⟩ in
the surface, twist 0, assumed to be their rotation 0), a 300 nm layer;
concentration over the peak concentration:

| Case | Peak (nm) | c(100 nm)/peak | Measured | c(200 nm)/peak | Measured |
|---|---|---|---|---|---|
| amorphous Si (trim2013) | 15.7 | 0 | < 1.7·10⁻⁵ | 0 | < 1.7·10⁻⁵ |
| crystal, 0° | 15.0 | 3.5·10⁻¹ | 5·10⁻³ | 1.8·10⁻² | 5·10⁻⁴ |
| crystal, 10° | 14.2 | 6.0·10⁻³ | 1.7·10⁻⁴ | 0 | 3.3·10⁻⁵ |
| 1 nm SiO2 on the crystal, 0° | 13.7 | 4.5·10⁻² | 5·10⁻³ | 7.2·10⁻⁴ | 5·10⁻⁴ |
| 1 nm SiO2 on the crystal, 10° | 13.6 | 3.4·10⁻³ | 1.7·10⁻⁴ | 4.8·10⁻⁶ | 3.3·10⁻⁵ |

The peak (12-16 nm) and the fall of the tails with tilt are reproduced,
and the native oxide matters: one nanometre of amorphous oxide scatters
the slow As ions out of the ⟨100⟩ channel and lowers the 0° tail at 100 nm
eightfold. Still, with the oxide, the tails at 100 nm are 9 × (0°) and
20 × (10°) the measured ones, and at 10° the calculated tail ends before
200 nm, where 3·10⁻⁵ of the peak was measured: AMCSET channels too many ions
to about 100 nm and dechannels the deeper ones too abruptly. The causes in
reach of the model are the perfect crystal (10¹⁴ cm⁻² As is near the
amorphization dose of Si; Tasch et al. find the tails shrink with dose),
the beam divergence and wafer miscut (below 0.5° each, not modelled) and the
stopping division. Like the boron profiles, these data point to the dose
dependence of channeling (docs/CRYSTAL_DESIGN.md, phase 5) as the next step
for slow heavy ions; until then, channeling tails of heavy ions at doses
near amorphization are an upper estimate.

## 7. Ionization by slow recoils

The share of a recoil's energy that ends as ionization (the "quenching
factor" of particle detectors) is measured down to sub-keV recoils in Ge
and Si. It tests the electronic stopping of slow atoms and the energy
partition of the whole cascade (TRIM's ionization by the ion and by the
recoils).

AMCSET (`literature_cases.py quenching`): 500 PKAs of one energy, isotropic,
started 500 nm deep in a 1 µm pure target (a recoil source), full cascades,
TRIM's default Ed/Eb/Es, 4 seeds; the yield is (ionization by the ion +
by the recoils) / recoil energy. Once with the SRIM-2013 tables, once with
Lindhard-Scharff stopping (`--stopping lss`). Lindhard: k g(ε)/(1 + k g(ε)),
g = 3ε^0.15 + 0.7ε^0.6 + ε, ε = 11.5 E/keV Z^(-7/3) (Lindhard et al., Mat.
Fys. Medd. Dan. Vid. Selsk. 33, 10 (1963), in the form of Lewin and Smith,
Astropart. Phys. 6, 87 (1996)).

**Ge in Ge**: Bonhomme et al., Eur. Phys. J. C 82, 815 (2022),
arXiv:2202.03754, Table 5 (quenching factors between 0.16 and 0.23 from 0.4
to 6.3 keV, compatible with Lindhard k = 0.162 ± 0.004):

| E (keV) | Measured | AMCSET, SRIM-2013 | AMCSET, LSS | Lindhard k = 0.162 |
|---|---|---|---|---|
| 0.97 | 0.181 ± 0.013 | 0.128 | 0.264 | 0.175 |
| 1.17 | 0.192 ± 0.019 | 0.131 | 0.270 | 0.179 |
| 1.54 | 0.190 ± 0.011 | 0.137 | 0.279 | 0.186 |
| 1.58 | 0.193 ± 0.010 | 0.136 | 0.280 | 0.186 |
| 1.92 | 0.190 ± 0.020 | 0.140 | 0.286 | 0.191 |
| 1.99 | 0.187 ± 0.006 | 0.141 | 0.287 | 0.192 |

**Si in Si**: Chavarria et al., Phys. Rev. D 94, 082007 (2016),
arXiv:1608.00957, Table 1 (C; E_e/E_r, the uncertainty from that of E_r);
Izraelevitch et al., JINST 12, P06014 (2017), arXiv:1702.00873, Table 1 (I;
systematic uncertainties). Lindhard with k = 0.15, the value Chavarria et
al. compare with:

| E (keV) | Measured | AMCSET, SRIM-2013 | AMCSET, LSS | Lindhard k = 0.15 |
|---|---|---|---|---|
| 0.68 (C) | 0.088 ± 0.013 | 0.261 | 0.223 | 0.203 |
| 1.05 (C) | 0.114 ± 0.011 | 0.281 | 0.241 | 0.216 |
| 1.36 (C) | 0.132 ± 0.008 | 0.294 | 0.252 | 0.224 |
| 1.66 (C) | 0.145 ± 0.005 | 0.303 | 0.261 | 0.230 |
| 1.94 (C) | 0.155 ± 0.005 | 0.311 | 0.269 | 0.236 |
| 2.28 (C) | 0.171 ± 0.002 | 0.319 | 0.275 | 0.242 |
| 1.79 (I) | 0.168 ± 0.034 | 0.307 | 0.264 | 0.233 |
| 2.97 (I) | 0.215 ± 0.020 | 0.333 | 0.289 | 0.252 |
| 4.87 (I) | 0.258 ± 0.011 | 0.361 | 0.313 | 0.272 |
| 7.67 (I) | 0.280 ± 0.006 | 0.388 | 0.339 | 0.294 |
| 10.77 (I) | 0.308 ± 0.005 | 0.411 | 0.360 | 0.312 |
| 14.48 (I) | 0.329 ± 0.004 | 0.434 | 0.381 | 0.331 |
| 20.67 (I) | 0.360 ± 0.004 | 0.462 | 0.408 | 0.355 |

The standard errors of AMCSET's values are 0.001 or less.

What this shows:

- With SRIM-2013's stopping AMCSET gives too little ionization in Ge
  (0.13-0.14 against 0.18-0.19, 25-32 % less) and too much in Si (+28 % at
  20 keV, 3× at 0.7 keV).
- Lindhard's integral theory with the measured k describes Ge and, above
  about 4 keV, Si (Izraelevitch et al. say the same); below that the Si
  data fall below every model.
- With Lindhard-Scharff stopping in the cascade AMCSET overshoots both,
  by 10-55 % (Si at 20 keV to Ge at 2 keV): the binary-collision cascade
  books more ionization than Lindhard's integral theory although its
  stopping is nearly the same (AMCSET's Lindhard-Scharff k for these pairs
  is close to the k values above). The cause was not traced; candidates
  are the treatment of the slowest atoms (TRIM's stop energy, Ed and Eb)
  and of the electronic loss in close collisions.

So no stopping switch fixes this: the ionization of keV recoils from
AMCSET (like TRIM's) is uncertain by tens of per cent, and up to a factor
of three for sub-keV Si. Use the measured yields for detector response.
Above some tens of keV per recoil, where most of the ionization comes from
fast atoms, this matters less (the TRIM comparison, §2, checks only that
AMCSET books ionization as TRIM does).

## 8. The native engine

The native engine (AMCSET's original code, ions and electrons) has not been
validated against measurements or another code. Its physics bugs P1-P6 were
fixed in 0.1.0 (CHANGELOG.md); docs/KNOWN_ISSUES.md §1 lists what is known.
Use SRIM mode for ions.

## 9. Deviations and the corrections that can be switched on

| Deviation | Cause | What AMCSET offers | Default |
|---|---|---|---|
| Ni ions go 1.2-1.9 × too deep (§4.1) | SRIM-2013's Ni stopping 30-40 % below Co and Cu at 0.02-0.3 MeV/u (§3.1) | `physics.stopping_corrections.nickel_ions` (CLI `--correct-nickel`, GUI Calculation tab): Ni projectiles take the geometric mean of Co's and Cu's electronic stopping at the same velocity, in every target, also as recoils. Peaks at 3.5-9 MeV then within 9 %; 2 MeV still 1.2-1.8 × deep | off (SRIM mode reproduces TRIM) |
| Au's stopping 26 % below its neighbours (§3.1), and 7-26 % below Au measurements in the same band (§3.2) | SRIM's Au scaling | none: Au ranges in light targets are already too short (§4.3) | — |
| Heavy ions in Si: single ions ±10-22 % (§4.3), Sb +22 % at 50 keV (§4.2) | the interatomic potential at low energy; transport of the implanted atoms (Wittmaack, §4.4) | none: no correction fits all ions; the deviation is documented per case | — |
| Heavy ions in light targets (C, B, SiC) 12-25 % too short, straggling 16-40 % too small (§4.3) | SRIM's stopping of slow heavy ions in light targets (Behar et al.: correlated inelastic losses; Zhang et al.: electronic stopping too high) | none in AMCSET | — |
| Sputtering yields off by a pair-specific factor 0.68-1.9; no yield near light-ion thresholds (§4.5) | TRIM's surface and cascade model | the surface binding energy Es per element | TRIM's Es |
| Channeling tails of slow heavy ions at doses near amorphization 9-20 × too high at 100 nm (§6) | perfect crystal; beam divergence and miscut not modelled | an amorphous surface layer (the native oxide) brings 0° within 10 × | — |
| Ionization of keV recoils (§7) | binary-collision energy partition | none: use measured yields | — |
| TRIM's full-cascade vacancies are not surviving defects (§5) | no recombination in a BCA | arc-dpa and rpa in the `dpa` table (docs/OUTPUT.md "dpa") | reported alongside |
| Replacements in Xe gas +4.6 % against TRIM (§2.2) | open; suspect: TRIM's single-precision recoil slots | — | — |
| TRIM.exe's own bugs (docs/SRIM_MODE.md §7) | TRIM.exe | 11 switches (`physics.trim_bugs`) | reproduced |

**How `nickel_ions` was derived.** The evidence is (1) the table scan (§3):
Ni is the one ion whose stopping departs from both pairs of neighbours by
far more than the Z1 oscillations, by the same factor in every target;
(2) the measured ranges (§4.1): SRIM-2013's Ni ranges are too long by about
what that dip gives; (3) Gurbich et al.'s observation that SRIM-2013's Co
behaves as TRIM-98's Ni, which agreed with their data. The correction
replaces only what is out of line (Ni as the projectile; Ni as a target is
untouched) with the smooth Z1 trend of SRIM-2013's own tables, the
geometric mean of Co and Cu at the same velocity. It shortens 5 MeV Ni in
Mo by 1.275 × (measured: 1.3) and 11.5 MeV Ni in Fe by 1.273 × (measured:
1.5). It also smooths whatever Z1 oscillation Ni really has below 0.02
MeV/u, where the ratio is 0.94-0.95, a few per cent of the stopping.

Tested since against what it was not built from: ten measured peak depths
of Ni in Ni, Fe and two steels at 2-11.5 MeV (§4.1; six within 9 % with the
correction, all 1.2-1.9 × too deep without it), and the measured Ni-ion
stopping of the IAEA database (§3.2; the one data set that measured Ni with
its neighbours agrees with the corrected values, an older one does not).
What it does not fix: 2 MeV (1.2-1.8 × too deep with or without it), and
5 MeV in Ni and 11.5 MeV in Fe (21-22 % too deep). A correction fitted to
these peaks, or TRIM-98's Ni stopping (which Song et al. and Chauhan et al.
recommend; its tables are in Song et al.'s supplementary material, not
obtained), would be the next step; it is not made here because the same
ten peaks would then be both the fit and the test.

## 10. Not covered yet

- **Sputtering** beyond yields at normal incidence (§4.5): angular
  dependence, energy distributions of sputtered atoms, compounds. The yields
  are read off figures (20-30 %); tabulated values (Behrisch and Eckstein,
  Sputtering by Particle Bombardment, Springer 2007) were not available.
- **Backscattering and reflection** coefficients against measurements, and
  **energy-loss and angular straggling** of transmitted ions: only against
  TRIM.exe (§2). Measured reflection coefficients of H and He were not found
  in open sources (candidates: Thomas et al., INDC(NDS)-287 (1993); Eckstein
  et al., IPP 9/82 (1993)).
- **Stopping in compounds** against the IAEA database's compound
  measurements (§3.2 covers elemental targets only).
- **Ni ions**: the full Gurbich et al. (2025) paper (behind a paywall) and
  Song et al.'s TRIM-98 stopping tables (their supplementary material) were
  not obtained; the latter would allow a TRIM-98 option to be compared with
  `nickel_ions` (§9).
- **More TRIM.exe seeds** for the Ni-ion cases (one run each, §2.3), to
  settle the ~1 % range excess of the 5 MeV cases in Fe.
- **Damage against experiments**: Averback, Benedek and Merkle (1978) was
  not read (§5 item 4); defect ranges by TEM (Lu et al., Sci. Rep. 6, 19994
  (2016): 3 MeV Au in Ni, damage seen to 1.1 µm where SRIM's ends near
  0.5 µm) are a matter of defect migration, outside a binary-collision code.
- **Electrons** (native engine): §8.
