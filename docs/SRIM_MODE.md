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
open difference, replacements in Xe gas, docs/VALIDATION.md §2.2), at a small fraction of
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

SRIM-2013's tables have one known error that AMCSET can correct: the
electronic stopping of **Ni ions** is up to 40 % below that of Co and Cu at
0.02-0.3 MeV/u, in every target, and measured Ni peak depths are 1.2-1.9
times shallower than SRIM-2013's. `stopping_corrections.nickel_ions` (off by
default, since TRIM.exe has the error) replaces it with the geometric mean
of Co's and Cu's at the same velocity, which brings the measured peaks at
3.5-9 MeV within 9 %: docs/VALIDATION.md §3, §4.1 and §9.

The electronic loss of a flight is S(E) L at the energy before the flight.
**Straggling**: `trim2013` adds (R − 1/2) Ω √L, limited to ±S L, with
Ω = min(√12 × relativistic Bohr, k_low E) (TRIM's SELOSS); `refined` uses
Ω²/L = min(relativistic Bohr, (Z1+Z2)^(8/3) E[keV] N[Å⁻³] / (3.2 M1)) eV²/Å
(SRIM book Eq. 7-47) as a truncated Gaussian or a gamma deviate
(`straggling_distribution`). The crystal transport uses the same variance,
always as a gamma deviate: its steps of about 1 Å have a straggling larger
than the mean loss, where a Gaussian cut at zero raises the mean.

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
model agrees with MD and experiments: docs/VALIDATION.md §5.

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
| `stopping_corrections.nickel_ions` | `false` | Not TRIM: Ni projectiles (and recoils) take the geometric mean of Co's and Cu's electronic stopping at the same velocity instead of SRIM-2013's Ni stopping, which is up to 40 % too low (ranges 1.2-1.9 × too long; docs/VALIDATION.md §9). SRIM tables only. CLI `--correct-nickel` |
| `straggling_distribution` | `truncated-gaussian` | `refined` only: or `gamma` (the crystal transport always uses `gamma`) |
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
AMCSET's ranges agree with TRIM's within the statistics (docs/VALIDATION.md §2.2). If the file
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

## 6. Validation

Moved to **docs/VALIDATION.md**, which holds every comparison of AMCSET
with TRIM.exe (§2 there), with measurements (stopping, ranges, ionization,
crystals) and with molecular dynamics (damage models, §5 there), with the
parameters and scripts of each.

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
