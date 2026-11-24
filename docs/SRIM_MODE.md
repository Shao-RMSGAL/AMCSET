# SRIM mode (TRIM-compatible engine)

AMCSET's **SRIM mode** is a binary-collision Monte Carlo engine written to
reproduce SRIM-2013 TRIM: same physics models, same input parameters, same
tallies and the same output files. ("SRIM" and "TRIM" are used
interchangeably in the field; strictly, TRIM is the transport program in the
SRIM package.) By default it follows TRIM.exe's own algorithm step by step
(`physics.transport: trim2013`), as established by disassembling TRIM.exe;
`docs/TRIM2013_ALGORITHM.md` describes that algorithm with the addresses of
the routines it comes from. Random numbers differ from TRIM's, so individual
ions differ, but ranges, damage, sputtering and the energy partition agree
with TRIM within the statistical uncertainty of the comparison (§6; one
open difference, replacements in Xe gas, §6.2), at a small fraction of
TRIM's run time (§8).

SRIM mode is the default engine (`"mode": "srim"` in a scenario, the
"Engine" box in the GUI). The original AMCSET engine is still available as
the **native** mode (it also transports electrons); see `AGENTS.md`.

Code: `common/include/trim/` and `common/src/trim/` (engine, TRIM.IN
reader/writer, TRIM-format output), `common/include/physics/` and
`common/src/physics/` (stopping and scattering, shared with the native
engine and the CLI).

## 1. Using it

- GUI: choose *SRIM mode* (Calculation tab), set Ed/Eb/Es per element in the
  layer tables, pick the damage calculation and the stopping source, run.
  The plots of TRIM's calculation window are live while it runs; the files
  go to the output folder (docs/GUI.md, docs/OUTPUT.md).
- CLI: `amcset run scenario.json`, or `amcset trim TRIM.IN` to run an
  existing TRIM input file directly (`docs/CLI.md`).
- Library: `trim::TrimSimulation(TrimConfig)`; `run()`, `tallies()`,
  `write_outputs(dir)`, `summary_json()`.

To match TRIM, use SRIM's stopping tables:

```
amcset trim TRIM.IN --stopping srim --stopping-table srim2013_stopping.amcset -o out
```

Without them AMCSET uses its built-in stopping (§2.2); results are then
"TRIM-like" rather than TRIM-identical. (TRIM.exe uses the ZBL universal
potential for every collision and never reads SRIM's pair potentials, so
`--srim-data` is not needed.)

## 2. Physics, as implemented

### 2.1 Two transports

| | `trim2013` (default) | `refined` |
|---|---|---|
| Ion flight path | TRIM-85 rule with the layer's mean Z and M (E_min = 5 eV), no limit on the energy lost per flight; first flight one monolayer; every flight a random fraction R L (at least a monolayer) in full-cascade mode | FFP from the largest p transferring 2 eV, limited to a 5 % electronic loss |
| Collision energy | after the flight's electronic loss | before it |
| Straggling | uniform deviate with Bohr's variance (TRIM) | Gaussian or gamma (SRIM book Eq. 7-47) |
| Recoils | monolayer steps, TRIM's slot stack and bookkeeping | free flights like the ion |
| Surface | sub-threshold atoms near the surface set in motion; atoms below the barrier reflected | optional surface-atom branch |
| Kinchin-Pease | layer-averaged Robinson partition (TRIM) | per-element partition |

`trim2013` reproduces TRIM's results, including TRIM's own quirks (energy
booked twice in a few places, interstitials that exceed vacancies by the
replacements, the direction of atoms struck by recoils; see
`docs/TRIM2013_ALGORITHM.md`). `refined` is AMCSET's earlier SRIM-mode
transport: its flights converge better and it conserves energy exactly, but
it departs from TRIM where TRIM's steps are coarse (slow light ions,
near-surface cascades and sputtering).

### 2.2 Electronic stopping

Models (`stopping.model`):

| Model | Description |
|---|---|
| `srim` | SRIM-2013 tables (SR Module output for every projectile/target-element pair), interpolated log-log as TRIM does; velocity scaling for other projectile masses; Bragg's rule times the layer's compound correction. Above the table (2 GeV of projectile energy: 2 GeV/u for H, 37 MeV/u for Fe) the curve continues with Bethe-Bloch, scaled to be continuous. |
| `builtin` | No data files: Lindhard-Scharff and a Bethe-Bloch form with a low-velocity term joined as S = (S_LS^-4 + S_BB^-4)^(-1/4). Within ~10 % of SRIM for H, ~25 % for He, up to 2x high for slow heavy ions. |
| `bethe-bloch` | Relativistic Bethe-Bloch with W_max, Barkas-type effective charge z* = Z1(1 - exp(-125 β Z1^(-2/3))) and the Sternheimer-Peierls density effect (high-energy particles, up to 100 GeV). |
| `lss` | Lindhard-Scharff only (low-energy heavy ions). |

The electronic loss of a flight is S(E) L at the energy before the flight.
**Straggling**: `trim2013` adds (R − 1/2) Ω √L, limited to ±S L, with
Ω = min(√12 × relativistic Bohr, k_low E) (TRIM's SELOSS); `refined` uses
Ω²/L = min(relativistic Bohr, (Z1+Z2)^(8/3) E[keV] N[Å⁻³] / (3.2 M1)) eV²/Å
(SRIM book Eq. 7-47) as a truncated Gaussian or a gamma deviate
(`straggling_distribution`).

### 2.3 Nuclear scattering

- ZBL universal potential for every collision, screening length
  a_U = 0.8854 a0 / (Z1^0.23 + Z2^0.23), as TRIM.exe (it never reads
  POTCOEF03.DAT). SRIM's pair-specific potentials can be used for recoils
  (`physics.pair_potentials: recoils`) or for all collisions (`all`) with
  `srim_data`; neither is TRIM's behaviour.
- Scattering angle: Biersack-Haggmark "magic formula"
  (`physics.scattering: magic`; `integral` evaluates the classical
  scattering integral by Gauss-Mehler quadrature instead). For ε > 10 TRIM
  switches to screened Rutherford scattering,
  sin²(Θ/2) = 1 / (1 + (1 + b(1+b))(2εb)²), and so does AMCSET.
  Note: the book prints the magic formula's G with numerator and
  denominator exchanged; the form used is G = 1/(γ(√(1+A²) − A)), which
  agrees with the scattering integral to 1-2 %.
- Nuclear stopping for the stopping tables and the CLI: ZBL universal fit,
  or SRIM's SNUC03 pair fits with `srim_data`.

### 2.4 Flight paths

`trim2013`: L = 1/(π p_max² N) with p_max = a / (EEG + √EEG + 0.125 EEG^0.1),
EEG = √(ε ε_min) and ε_min from E_min = 5 eV, using the layer's mean mass and
integer mean Z; flights that would cross a layer boundary, leave the depth
window or the surface are shortened by random factors. Impact parameter
p = √(R/(π L N)). Recoils move one monolayer L0 = N^(-1/3) per step with
p = √R L0/√π (smaller for light atoms at the surface). Details:
`docs/TRIM2013_ALGORITHM.md`.

`refined`: the largest impact parameter that transfers more than T_min =
`physics.min_energy_transfer` (2 eV), p_max = max_j a_j b_max(ε_j),
b_max = 1/(ξ + ξ^0.5 + 0.125 ξ^0.1), ξ = √(ε ε_min/γ) (SRIM book Eqs. 7-28,
7-40); L = 1/(π p_max² N), limited so that the electronic loss over L is
≤ 5 % of E; `physics.max_impact_parameter_ws` can cap p_max; near a surface
or interface the flight is truncated to a random fraction of the distance;
`physics.free_flight_path: false` uses monolayer steps.

`damage: monolayer` (TRIM's "surface sputtering / monolayer collision
steps", cascade code 3) uses monolayer flights for the ion as well.

### 2.5 Damage

- **Full cascades** (`damage: full-cascade`, TRIM's "Detailed calculation
  with full damage cascades"): a struck atom with T ≥ Ed is displaced (a
  vacancy) and leaves with T − Eb (TRIM subtracts Eb only deeper than two
  monolayers). Atoms below Ed are not followed, except within a random
  0-9 monolayers of the surface, where they are set in motion without being
  displaced (and can be sputtered). A displaced mover that stops right after
  striking a displaced atom of its own element replaces it (the vacancy is
  removed). Ions and recoils stop at E_f = `physics.final_energy` (1 eV, or
  `ion_stop_energy` / `recoil_stop_energy` if set) and become interstitials;
  as in TRIM, a replacing atom is still counted as an interstitial.
- **Quick calculation** (`damage: kinchin-pease`): primary recoils are not
  followed; each gets ν = 1 if Ed ≤ E_d(T) < 2Ed/0.8 and ν = 0.8 E_d(T)/(2Ed)
  above, with Robinson's damage energy E_d(T) = T/(1 + k_d g(ε_d)),
  g = ε + 0.4 ε^(3/4) + 3.4 ε^(1/6). `trim2013` uses TRIM's constants with
  the layer's integer mean Z and mean M (ε_d = 0.01 Z^(-7/3) T[eV],
  k_d = 0.1334 Z^(2/3)/√M); `refined` uses the struck element's
  (ε_d = 0.01014 Z^(-7/3) T). TRIM applies the displacement threshold to
  the **damage energy**, not to T (for Fe with Ed = 40 eV the smallest
  displacing recoil in COLLISON.txt is 45.27 eV), and AMCSET does the same.
  T − E_d goes to recoil ionization, E_d to recoil phonons.
- Per-element and per-layer Ed, Eb and Es, as in TRIM. TRIM.exe replaces a
  displacement energy of 0 (or less) by 20 eV, a surface binding energy by
  4 eV and a lattice binding energy by 2 eV, but only for its target atoms
  2..n (its loop starts at the second atom; TDATA.txt shows the values it
  used). A TRIM.IN with Eb = 0 for an alloy therefore runs with Eb = 2 eV
  for every element but the first. `trim2013` does the same, in TRIM's atom
  order (kept from a TRIM.IN as `trim_atom`), and prints a warning.
- In the quick calculation TRIM subtracts Eb once per displaced atom from
  E2RECOIL's "from IONS" column (and ν Eb from the recoil phonons); its
  second column ("Absorbed by Recoils") keeps the whole transfer.

### 2.6 Sputtering

An atom reaching the front surface escapes if its normal energy E cos²θ is
at least Es (planar barrier); otherwise it is reflected back into the
target. `summary.json` and the yields give the sputtered energy after the
barrier (E − Es). SPUTTER.txt lists the atoms as TRIM does in `trim2013`:
the energy inside the surface ("SBE is NOT subtracted") and the direction
before the barrier; `refined` lists the energy and the refracted direction
after it. The file header says which.

### 2.7 Random numbers and parallelism

L'Ecuyer's combined multiplicative generator (1988), with two
independent streams per ion (ion and recoils) seeded from `seed` and the ion
number via SplitMix64. Ions are processed in fixed blocks of 16 that are
merged in order, so a result depends only on the seed and the number of
ions, never on the thread count.

### 2.9 The crystal transport (not TRIM)

`physics.transport: "crystal"` is a third transport, for layers that have a
crystal (`layers.N.crystal`): collision partners on lattice sites with
thermal vibrations, simultaneous collisions, electronic stopping divided
between the path and the collisions, after Robinson and Torrens (1974). It
shares the stopping powers, the potential, the tallies and the outputs of
SRIM mode but reproduces nothing of TRIM.exe. It is compared with a
measured channeling dip and measured boron profiles in silicon, which
bracket its division of the electronic stopping without fixing it. What it
does, its settings, its checks and those comparisons:
docs/CLI.md "Crystals"; the design and its differences from the paper:
docs/CRYSTAL_DESIGN.md §11. The comparisons with TRIM in §6 concern
`trim2013` only.

### 2.8 Damage models and dpa (not TRIM)

Besides TRIM's vacancies, `trim2013` books four damage functions of every
primary knock-on atom's damage energy (NRT, Kinchin-Pease, arc-dpa, rpa)
by depth and converts them to dpa at `output.dpa_fluence` (the `dpa` table
and `summary.json`'s `damage_models`; definitions and limits in
docs/OUTPUT.md "dpa"). The bookkeeping only adds to its own tallies: no
random number is drawn and every TRIM output stays bit-identical. Which
model agrees with MD and experiments: §6.3.

## 3. Parameters (`physics` in a scenario)

| Key | Default | Meaning |
|---|---|---|
| `transport` | `trim2013` | `trim2013` (or `trim`): TRIM.exe's algorithm; `refined`: §2.1 |
| `pair_potentials` | `none` | `none` (TRIM), `recoils`, `all` (needs `srim_data`) |
| `scattering` | `magic` | `magic` or `integral` |
| `electronic_straggling` | `true` | Electronic energy-loss straggling |
| `final_energy` | `1 eV` | E_f (stop energy) |
| `ion_stop_energy`, `recoil_stop_energy` | `0 eV` | Extra stop energies (0 = E_f) |
| `max_recoils_per_ion` | `0` | Cap on followed recoils per ion (0 = none) |
| `fast_math` | `false` | Fast approximations (§8.2): statistically the same results, about 2x faster, not number for number |
| `bethe_bloch.effective_charge`, `.density_effect` | `true` | Bethe-Bloch corrections |
| `straggling_distribution` | `truncated-gaussian` | `refined` only: or `gamma` |
| `free_flight_path` | `true` | `refined` only: FFP algorithm; `false` = monolayer steps |
| `min_energy_transfer` | `2 eV` | `refined` only: T_min of the FFP algorithm |
| `max_impact_parameter_ws` | `0` | `refined` only: cap of p_max in Wigner-Seitz radii (0 = none) |
| `kp_subtract_binding_energy` | `false` | `refined` only: Kinchin-Pease on T − Eb |
| `surface_atom_depth` | `0 A` | `refined` only: depth of the surface-atom sputtering branch |

## 4. SRIM stopping tables and data

AMCSET comes with the SRIM-2013 stopping tables of every ion (Z1 = 1-92)
in every element (Z2 = 1-92), solid and gas phase: 16928 tables in
`data/srim2013/srim2013_stopping.amcset` (42 MB of text). SRIM mode uses
them by default (stopping model `srim` with no table given): TRIM.exe
computes its stopping from the same SRIM-2013 coefficients
(`Data/SCOEF03.DAT`) as SR Module, which made the tables, and with them
AMCSET's ranges agree with TRIM's within the statistics (§6.2). If the file
is not found, a run warns and uses the built-in stopping, which differs
from TRIM's by up to 25 %. It is looked for in `$AMCSET_DATA_DIR/srim2013/`,
in `data/srim2013/` next to the executable or one or two folders above it,
in `share/amcset/srim2013/` (installed), and in the source tree. SRIM is
free for non-commercial use; the tables carry the same terms
(`data/srim2013/README.md`).

`tools/srim/` contains the scripts that made them from a SRIM-2013
installation (Windows):

1. `gen_cases.py` writes `sr_jobs.csv` (projectile/target pairs) and the
   TRIM validation cases.
2. `extract_stopping_tables.ps1` runs SR Module for every pair (in parallel
   copies of the "SR Module" folder).
3. `pack_stopping.py` packs the outputs into one text file
   (`srim2013_stopping.amcset`). A folder of raw SR Module outputs named
   `sr_ZZZ_ZZZ_s.txt` also works as `stopping.table`.

The optional pair potentials and nuclear-stopping fits are read directly
from SRIM's `Data` folder (`POTCOEF03.DAT`, `SNUC03.dat`), given as
`srim_data`.

See `tools/srim/README.md` for details.

## 5. Tallies

Everything TRIM tallies, per depth bin (100 bins over the target or the
`output.depth_window`) and in total: ion and recoil final positions; range
moments (mean, straggling, skewness, kurtosis), lateral and radial ranges;
vacancies (by element; by ion and by recoils in KP mode), replacements,
displacements, interstitials, knock-ons; ionization and phonons by ions and
by recoils; energy to recoils and energy absorbed per element;
backscattered, transmitted and sputtered atoms with their energies,
positions and directions; collision details of the first ions. AMCSET adds
recoil energy spectra (PKA and all recoils) and `summary.json`.

TRIM.exe sums its count tables (vacancies and knock-ons in VACANCY.txt,
replacements in NOVAC.txt, stopped ions and atoms in RANGE.txt) in single
precision, adding 1/(bin width) per event to a VB `Single`. Once a bin's sum
is large, every addition is rounded to the sum's precision, so the tables
are off from the true counts by up to about ±1.5 %, depending on the bin
width and the number of events (+0.6 % for 2000 Fe ions at 2 MeV into Fe,
−0.2 % at 1 MeV, +1.3 % for W → W at 1 MeV); the totals in the file headers
are exact. AMCSET's TRIM-format files reproduce this
(`output.trim_single_precision`, default `true`); its CSV, JSON and Parquet
tables and `summary.json` hold the exact values.

TRIM's definitions are kept where they define a tally, so that numbers can
be compared directly (e.g. E2RECOIL counts only transfers that displace an
atom, and TRIM's IONIZ and PHONON tables add up to slightly more than the
incident energy because TRIM books some energy twice; `trim2013` does the
same).

## 6. Validation against SRIM-2013 TRIM

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

### 6.1 Differences found after the first round, and their causes

The second round (§6.2) traced the differences that were left in the
table above:

- **Alloys and layered targets, full cascades** (replacements +6 % for
  Fe → Cr on Fe, +2 % for 316 SS): TRIM replaces a zero lattice or surface
  binding energy of every target atom after the first by 2 eV and 4 eV
  (0x44efc0, §2.5), which changes the replacement collisions of the later
  elements. Now within statistics.
- **H → Fe 3 MeV, quick calculation** (vacancies −1.6 %): TRIM adds its
  vacancy and range tables in single precision with a weight of
  1/(bin width), which shifts the totals of long runs by up to ±1.5 %
  (predicted +0.64 % for Fe → Fe 2 MeV, observed +0.63 %). AMCSET
  reproduces the rounding in the TRIM text tables (§5).
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

### 6.2 Second round (November 2026)

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


### 6.3 Damage models against MD, experiments and SRIM post-processing

Not a TRIM comparison: these runs test the `dpa` table (§2.8) and ask which
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
for 3-10 MeV Al, Fe and Au (Ed = 29 eV, §6.3 item 1); per PKA 0.20 at
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

## 7. Bugs of TRIM.exe

Reverse engineering TRIM.exe (`docs/TRIM2013_ALGORITHM.md`) turned up
behaviour that is wrong rather than a modelling choice. The trim2013
transport reproduces all of it by default, so that AMCSET's results agree
with TRIM's; each bug can be corrected separately (`physics.trim_bugs` in a
scenario, `--fix-trim-bug NAME` / `--fix-trim-bugs` on the command line, the
"TRIM.exe bugs" group on the GUI's Calculation tab, whose tooltips give the
details below). **A run with a corrected bug no longer reproduces TRIM: do
not compare it with TRIM runs or published TRIM numbers.** The refined
transport has none of these bugs. Nor has the crystal transport (§2.9),
which is not TRIM: it runs with every one corrected whatever the switches
say (`TrimConfig::effective_bugs`), with exact sums in its TRIM-format
tables, and records that in results.json; the scenario keeps its switches
for the trim2013 transport.

| Name | TRIM.exe (address) | Corrected | Effect of the correction |
| --- | --- | --- | --- |
| `stale_primary_energy` | An ion collision with T <= Ed displaces nothing, yet the cascade runs with the energy recoil slot 1 had at the end of the previous cascade (STARTRECOIL) | No cascade; T stays in the ion's phonons | Recoil phonons -8 to -11 %, recoil ionization -6 to -9 %, sputtering -3 to -5 % at 1 keV (Fe -> Fe, Ar -> Cu); about -1 % at 1 MeV |
| `recoil_energy_created` | A recoil's electronic loss DEE over a step and the transfer T are both computed at its energy before the step, can together exceed it, and are booked in full; the mover keeps max(0, E - T - DEE) (0x436de5). Fe -> Fe at 1 keV deposits 105 % of the incident energy | T <= E - DEE | Replacements -3 to -7 %, sputtering -3 to -7 % at 1-10 keV, recoil ionization up to -4 %, vacancies up to -1 % |
| `struck_atom_direction` | The atom struck by a recoil is turned by DIRCOS with sin x = 0 (0x4377d5): its direction depends only on the mover's x cosine | Turned from the mover's direction, like the primary of an ion collision | Sputtering -17 % (Ar -> Cu 1 keV), -4 % (Fe -> Fe 10 keV); vacancies up to -3 % |
| `replacement_interstitial` | An atom that makes a replacement collision is still counted as an interstitial (RECSTOPDATA, 0x457970) | Not an interstitial | Interstitials -2 to -5 % |
| `subthreshold_sputter_counted` | A sputtered atom that was set in motion below Ed adds a displacement and a vacancy to the totals (not to the depth tables) (SPUTTER, 0x441ff0) | Not counted | Total vacancies and displacements -9 % (Fe -> Fe 1 keV), -17 % (Ar -> Cu 1 keV), -2.5 % at 10 keV, none at MeV energies |
| `ion_phonons_counted_twice` | An ion transfer below the stop energy is booked as ion phonons and again as recoil phonons (0x44d4e0, 0x4575e0) | Once | Recoil phonons -0.02 % |
| `cascade_random_replay` | The random-number state is saved before and restored after every cascade (0x434e89, 0x438140): cascades replay the ion's next numbers and are correlated | One stream per cascade | Totals within statistics (< 1 %); TRIM's run-to-run scatter is larger than independent cascades would give |
| `zero_energy_defaults_skip_first_atom` | Zero Ed, Es, Eb become 20, 4, 2 eV for target atoms 2..n; the loop (0x44efc0) skips atom 1 | Every atom | Only inputs with zero energies |
| `integer_layer_z` | The layer's mean Z is rounded at every element, z = CInt(c Z + z) (0x44f8c0) | Stoichiometric mean | Compounds and alloys, below 0.5 % in the cases tried |
| `exyz_skips_first_step` | EXYZ.txt: the next record after the start is set a step too low (0x434745) | Every crossing written | EXYZ.txt only |
| `recoil_event_counted_twice` | TRIM.DAT recoils (code 6) count each event twice in RANGE.txt's ion column | Once | That table only |
| `single_precision_tables` (`output.trim_single_precision`) | VACANCY/NOVAC/RANGE tables summed in single precision | Exact sums | Table totals of long runs up to +-1.5 %; TRIM-format tables only |

The effects are AMCSET runs with that one correction against runs with
none (three seeds of 200-5000 ions each, the inputs of the TRIM validation
cases; MeV sputtering yields are too noisy to quote). With every
correction, Fe -> Fe at 1 keV has 12 % fewer vacancies and 15 % less recoil
phonon energy, and Ar -> Cu at 1 keV 19 % fewer vacancies and 23 % less
sputtering; at 1 MeV the vacancies change by -3 % and the replacements by
-13 %.

Not reproduced: TRIM keeps its recoil slots in single precision; AMCSET
keeps double precision (the difference is far below the statistics in
solids; see the open Xe gas case in section 6.2). TRIM also discards the
recoils of a cascade beyond the stack size set in Data/TRIM.cfg ("Minor
TRIM Error: Cascade > N Atoms", RECOILLOOP before 0x437de4); AMCSET's stack has no limit.

## 8. Performance

### 8.0 Against TRIM.exe

Measured on the validation runs of section 6.2 (AMCSET 0.1.0, before the
optimisations of §8.1, which make AMCSET's transport 1.35 times as fast
and cut its fixed cost per run from about 0.9 s to about 0.2 s), on one PC (AMD Ryzen 7
7800X3D, 8 cores / 16 threads, Windows). TRIM.exe ran in batch mode
(TRIMAUTO = 1), one instance per case and about 14 instances at a time
(TRIM is single-threaded); the runner's start and end times give the number
of instances alongside each run. AMCSET ran one case at a time on all 16
threads, through `tools/trim_compare/batch.py`; its times include start-up
(loading the 42 MB stopping file and writing TRIM's output files, about
0.9 s per run). Speed-up per instance = TRIM time / AMCSET time; speed-up
per machine = TRIM time / (instances running) / AMCSET time, which charges
TRIM only for its share of the PC.

| Case | Input | Ions | TRIM (s) | TRIM instances | AMCSET (s) | Per instance | Per machine |
|---|---|---|---|---|---|---|---|
| Fe→Fe 1 keV, FC | `434_rep_fe_fe_1keV_fc_s1` | 5000 | 42 | 15.6 | 0.91 | 46× | 2.9× |
| Fe→Fe 10 keV, FC | `422_rep_fe_fe_10keV_fc_s1` | 5000 | 132 | 14.7 | 1.11 | 119× | 8.1× |
| Fe→Fe 100 keV, FC | `377_rep_fe_fe_100keV_fc_s1` | 5000 | 1028 | 14.0 | 2.46 | 418× | 29.9× |
| Fe→Fe 1 MeV, FC | `328_rep_fe_fe_1000keV_fc_s1` | 2000 | 1784 | 14.0 | 4.83 | 370× | 26.4× |
| Fe→Fe 2 MeV, FC | `307_rep_fe_fe_2000keV_fc_s1` | 2000 | 3008 | 14.0 | 6.39 | 471× | 33.7× |
| Fe→316SS 2 MeV, FC | `330_rep_fe_316ss_2000keV_fc_s1` | 1000 | 1186 | 14.0 | 3.72 | 319× | 22.8× |
| Fe→Fe 1 MeV, KP | `324_rep_fe_fe_1000keV_kp_ed25eb3_s1` | 10000 | 2407 | 14.0 | 1.12 | 2155× | 154.0× |
| H→Fe 1 MeV, KP | `341_rep_h_fe_1000keV_kp_s1` | 20000 | 686 | 14.0 | 0.93 | 742× | 53.0× |
| He→W 2 MeV, KP | `321_opt_he_w_2000keV_kp` | 20000 | 2106 | 14.0 | 0.90 | 2339× | 167.1× |
| Ar→Ar gas 100 keV, FC | `385_gas_ar_ar_100keV_fc_gas` | 1000 | 8212 | 6.4 | 28.99 | 283× | 44.2× |
| Xe→Xe gas 1 MeV, FC | `367_gas_xe_xe_1000keV_fc_gas` | 500 | 28587 | 2.9 | 158.10 | 181× | 61.7× |
| Fe→Fe TRIM.DAT ions, FC | `386_dat_fe_fe_ions_fc` | 1000 | 6993 | 7.0 | 1.55 | 4512× | 642.7× |

Over all 75 cases TRIM used 4.7 hours of the PC (30.9 instance-hours)
and AMCSET 5.3 minutes: 53× less machine time. The geometric mean
speed-up is 30× per machine (0.5-643×) and 406× per
instance. Kinchin-Pease runs finish within AMCSET's start-up time, so their
speed-ups are lower bounds; short keV cascades that take TRIM under a
minute gain least (3-8× per machine), and the smallest runs (30-100 ions)
are dominated by start-up on both sides. AMCSET gives identical results on
any number of threads (ions run in blocks of 16, merged in order).

For a controlled measurement - TRIM alone on the machine, AMCSET on one
thread and on all threads, same inputs - run
`tools/windows/benchmark_vs_trim.bat` on an otherwise idle PC; it writes
`build/benchmark/benchmark.csv`.

### 8.1 Optimisations that keep the results (the default)

The default calculation gives bit-identical results to 0.1.0: the 13
reference cases of `tools/perf/bench.py` (Fe → Fe 1 keV to 1 MeV, 316 SS,
layers, Ar gas, TRIM.DAT ions and recoils, KP and full cascades; 377 output
files) are identical file by file (`bench.py compare`). Measured on the
idle workstation with `bench.py ab` (one thread, the builds interleaved,
fastest of 3 runs; `build/perf/final`), the transport (`summary.json`'s
time) is 1.35 times as fast as 0.1.0 (geometric mean; 1.29-1.40 times per
case) and the whole process 1.67 times (1.46-2.30): start-up, output and
the CLI's wait at the end fell from 0.93 s to 0.23 s per run (and to about
0.14 s once the CLI stopped polling every 200 ms).

Where the time went (sampling profile of Fe → Fe 100 keV, full cascades,
0.1.0): 28 % in `exp`, 12 % in the closest-approach iteration, 5 % in the
screening function's derivative and 4 % in the magic formula itself: about
half of the run computes scattering angles. A short run also spent 21 % of
its time parsing the stopping file. The changes:

- The screening function and its derivative from the same four
  exponentials (`ScreeningFunction::phi_and_derivative`, same expressions
  in the same order), inlined into the Newton iteration of the closest
  approach and the magic formula: half the `exp` calls.
- The live tallies of the GUI are no longer accumulated a second time per
  event; snapshots are derived from the per-block `TrimTallies`.
- Stopping tables: the 42 MB file is read in one piece and only the tables
  of the run's ion/element pairs are parsed, on first use (thread-safe), so
  start-up fell from about 0.8 s to about 0.05 s. This matters most for
  short runs (KP, low energies, the GUI's test of a scenario).
- Output files are written in parallel (each TRIM file and each table
  format on its own thread), and their numbers are printed with
  `std::to_chars` into large buffers (`output/text.h`, `format_number_to`)
  instead of `snprintf` and a `std::string` per number: the same characters
  (`OutputText.*` compare them with `snprintf` over random values, ties,
  denormals and extremes). For 2 million H ions on 16 threads (RANGE_3D.txt
  92 MB, CSV and JSON about 100 MB each) the time after the transport fell
  from 2.4 s to 1.0 s (single runs on a busy machine).
- Per-flight details: the depth-bin edges that spread a flight's energy
  over the bins are computed once per block (the same expression, so the
  same doubles), and the random-number generator reduces modulo
  m = 2^31 - c with a multiplication and one conditional subtraction
  instead of a 64-bit division (the same numbers,
  `TrimRandom.ModularReductionMatchesTheRemainder`).
- Release builds keep PDBs (`/Zi`, `/DEBUG /OPT:REF /OPT:ICF`) for crash
  reports and profiling, at no cost in speed.

Tried and not kept: `/arch:AVX2` (same speed: the hot loops are scalar
chains of `exp` and divisions), profile-guided optimisation (`AMCSET_PGO`,
+2 %, identical results; the option stays for experiments), polynomial
`exp`/`log` (slower than the C runtime's). GCC 15 (`-O3 -march=native
-flto`, `docs/DEVELOPMENT.md` §2; WSL Ubuntu on the same PC) runs the same code
6 % faster than MSVC 14.5 in the exact mode and 14 % faster in the fast
mode. (Its C library's `exp` and `log` differ from Microsoft's in the last
bits, so Linux and Windows runs of the same seed are not identical
number for number; each platform is deterministic.)

Parallelism: ions run in blocks of 16 on every core with results
independent of the thread count (§2.7); the table of §8.2 is built on all
cores; output files are written concurrently.

### 8.2 Fast approximations (`physics.fast_math`, off by default)

`physics.fast_math: true` (CLI `--fast`, GUI *Fast approximations*)
replaces the most expensive exact computations by close approximations.
The random-number draws and the algorithm are TRIM's, but scattering
angles differ in the last digits, so a history can eventually take another
course: results agree with the exact calculation statistically, not number
for number, and a fast run is not a step-by-step reproduction of TRIM.exe.
Use the exact mode for comparisons with TRIM and for reference results.

- **Scattering angle from a table** (`physics::MagicTable`): ln θ of the
  magic formula for the universal potential on 800 × 800 nodes, log-spaced
  in ε from 10⁻⁹ to 12 and spaced in √b up to b = 256 (gases reach about
  250), built once per process on all cores (5 MB), interpolated bicubically.
  The relative error is below 4 × 10⁻⁵ in sin²(θ/2) and 10⁻⁶ in cos²(θ/2)
  (measured 3.5 × 10⁻⁵ over 20 000 random (ε, b),
  `Scattering.MagicTableReproducesTheMagicFormula`). Outside the table,
  for SRIM pair potentials and for the scattering integral the exact
  formula is used; above ε = 10 TRIM's screened Rutherford form, as before.
- **Tables indexed by the bits of the argument** (`fast::OctaveTable`:
  128 bins per factor of two, read from the exponent and top mantissa bits,
  so no logarithm; interpolated linearly): the electronic stopping of every
  species in every layer (1 meV to 137 GeV, relative error below 4 × 10⁻⁵,
  `FastMath.TabulatedStoppingPowerMatchesTheCurve`), TRIM's Kinchin-Pease
  partition ε + 0.4 ε^3/4 + 3.4 ε^1/6 (2 × 10⁻⁶), the flight path's
  e + √e + 0.125 e^0.1 (3 × 10⁻⁶) and the log10 of the recoil spectra's
  bins: no `exp`, `log`, `pow` or `cbrt` per lookup.
- The collision azimuth's sine and cosine from a polynomial of the quarter
  turn (no branches; error below 10⁻¹¹), the scattering rotation with one
  reciprocal instead of five divisions and renormalised by one Newton step
  instead of a square root and a division, and one division per flight
  instead of one per depth bin when a flight's energy is spread over the
  bins.

Speed (13 reference cases, one thread, `bench.py ab`, idle machine),
measured before the tables of the stopping power, the Kinchin-Pease
partition and the flight path were added: 2.08 times the exact mode with
MSVC (geometric mean; 1.63-2.31 per case: about 2.25 for full cascades,
1.6-1.7 for the Kinchin-Pease runs, 1.85 for the Ar gas case), 2.81 times
0.1.0 (3.13 times in wall time). GCC 15 is a further 14 % faster
(an earlier interleaved measurement: 2.39 times MSVC's exact mode). The
tables removed what the profiles then showed in the fast mode: the
stopping power's `exp` and `log` (about 9.5 % of a full-cascade run), and
`cbrt` (11 %), the flight path's `exp` and `log` (8 %) and the stopping's
(4.5 %) in a Kinchin-Pease run of H ions; spreading flights over the depth
bins fell from 12 % to 5 % of that run. Their speed-up is to be measured on
an idle machine (`docs/ROADMAP.md`).

Accuracy (`tools/perf/fast_vs_exact.py`: the 13 cases at half their
benchmark ion counts, 6 seeds each in both modes, 421 quantities of
`summary.json`; with every approximation above): the fast and exact means
agree as independent samples of the same distribution should (median |z|
0.65, 0.67 expected; |z| > 2 for 20, about 19 expected). The four
quantities with |z| > 3 are one case's backscattering and ion ionization
(Fe → Fe 100 keV): the two seed sets differ there, because with the same
seeds fast and exact give the same backscatter count and ionization within
4 × 10⁻⁶. With the same seeds the two modes share most histories: the
relative difference of a quantity is 6 × 10⁻⁸ (median) and below
4.3 × 10⁻⁴ for 90 % of them (below 6 × 10⁻⁴ for every quantity of the ten
larger cases); it exceeds 10⁻³ only in the three smallest runs (6 to 200
ions), where one history that takes another course shifts the totals.
(Before the tables: 2 × 10⁻⁹ median, 8 × 10⁻⁵ for 90 %.) A comparison of the fast mode with
the TRIM reference runs (§6, about an hour) has not been made yet.

The quick test runs of the automatic target depth (docs/CLI.md "Checks and
warnings") always use the fast approximations: only the depth of the
deepest ion matters there.
