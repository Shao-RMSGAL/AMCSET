# common/src/trim — SRIM mode (TRIM-compatible engine)

Read `docs/SRIM_MODE.md` (what the engine does and how it is validated) and
`docs/TRIM2013_ALGORITHM.md` (TRIM.exe's algorithm with addresses) before
changing anything here.

## Files

| File | Contents |
|---|---|
| `trim_engine.cpp` | Configuration helpers (`to_string`, defaults, `TrimConfig::validate`, `TrimBugs::members`); `TargetModel` construction in the `TrimSimulation` constructor (TRIM's zero-energy defaults, per layer/element models, TRIM's layer averages `Trim2013Layer` = TRIMSTART, stopping curves, straggling tables, collision pairs); the `Transporter` class (both transports); the driver (`run`, `run_block`, snapshots, pause, cancel, `continue_from`) |
| `trim_output.cpp` | TRIM-format text files (RANGE, VACANCY, NOVAC, IONIZ, PHONON, E2RECOIL, LATERAL, TDATA, BACKSCAT, TRANSMIT, SPUTTER, TRIMOUT, RANGE_3D, COLLISON, EXYZ, PKA_SPECTRUM) and `summary.json`; `write_outputs` also writes CSV/JSON/Parquet through `common/src/output/` |
| `trim_tables.cpp` | `make_dataset`: TrimTallies → `output::Dataset` (every table and scalar, with dimensions and units). The GUI's plots and tables read this, so a new table appears in the GUI and in every output format |
| `trim_in.cpp`, `trim_dat.cpp` | TRIM.IN and TRIM.DAT (and BACKSCAT/TRANSMIT/SPUTTER as particle lists) readers and writers |
| `single_precision.cpp` | Reproduces TRIM's single-precision accumulation of its count tables (`output.trim_single_precision`) |
| `trim_export.h`, `trim_export.cpp` | Not TRIM: the trajectory file and the MD event file (`TrimConfig::trajectory_export`, `md_events`): their columns (keep `docs/OUTPUT.md` in step), the per-block rows and the writers. The events themselves are decided in the Transporter (`md_birth`, `md_crossing`, `md_ion_entry`) |
| `target_depth.cpp` | Not TRIM: the ions' path-length estimate (`estimate_penetration`, `estimate_deepest`), `TrimConfig::plausibility_warnings` and the automatic target depth (`size_target`: quick fast-math test runs that size the last layer) |

## Transporter map (trim2013, the default transport)

| TRIM.exe | Function | Notes |
|---|---|---|
| IONCYCLE, 0x433fd0 | `run_ion_2013` | start, per-flight loop, tallies (0x44d4e0), KP or cascade start |
| 0x44ab60 | `ion_flight_path_2013` | free flight path L = 1/(π p_max² N), randomisation counter |
| 0x44b330 / 0x44b7a0 | `interface_2013` / `edge_2013` | shortening at layer boundaries / the depth window |
| SELOSS 0x45d4e0 | `electronic_loss_2013` | S(E)·L + uniform straggling; draws a number every call |
| 0x45f7c0 | `pick_element` | draws a number even for pure elements |
| Tmagic 0x467600 | `physics::trim_sin2_half_theta` (`common/src/physics/scattering.cpp`) | magic formula; ε > 10 screened Rutherford |
| DIRCOS 0x45ce40 | `rotate`, `rotate_about_x_2013` | the second is TRIM's struck-atom bug |
| 0x44bae0 | `kp_2013` | modified Kinchin–Pease with layer-averaged Robinson partition |
| STARTRECOIL | `start_cascade_2013` | RNG save/restore (0x434e89, 0x438140) |
| RECOILLOOP / NEXTRECOIL | `cascade_2013` | recoil slot stack, monolayer steps, SPUTTER 0x441ff0 |
| RECSTOPDATA 0x457970 | `recoil_stop_2013` | replacements, interstitials |
| TRIM.DAT recoils | `run_recoil_event_2013` | cascade codes 6 (7 runs through the ion loop) |

`follow` / `cascade` / `run_ion` (without `_2013`) are the **refined**
transport (`physics.transport: refined`), AMCSET's earlier algorithm. It is
kept for comparison; it is not TRIM.

`run_ion_crystal` / `crystal_follow` / `crystal_strike` /
`crystal_path_loss` are the **crystal** transport (`physics.transport:
crystal`): the binary-collision scheme of Robinson and Torrens, Phys. Rev.
B 9, 5008 (1974), in layers that have a crystal (`LayerModel::lattice`, a
`crystal::LatticeSearch`), random partners in the others. It is AMCSET's
own physics: rule 1 (fidelity to TRIM.exe) does not apply to it, rules 3
(determinism: the lattice's configuration is a hash of the seed, the ion
and the site, never shared state) and 4 (units) do. Cite the paper's
section or equation for what follows it and say where it departs
(`docs/CRYSTAL_DESIGN.md` §11 lists the departures). A change to it must
keep `tools/crystal/channeling_checks.py`'s trends (random direction =
amorphous; the normalisation of the stopping division) and must not touch
the other two transports (`bench.py compare`).

## Rules

1. **Fidelity first.** The trim2013 functions reproduce TRIM.exe, including
   the order and number of random draws (the RNG sequence decides every
   later collision). Do not reorder draws, skip "unused" draws, merge
   branches, or change a constant TRIM uses (14.4, 0.5292, 0.8853, 931494,
   177, …) without evidence from TRIM.exe.
2. **TRIM bugs are switches.** To make TRIM's behaviour correct, never edit
   it in place: add a `TrimBugs` member (default `true` = as TRIM) and branch
   on `cfg_.bugs.<name>`. Recipe: `docs/HOWTO.md` "Add a TRIM bug switch".
   The 11 existing switches and their measured effects:
   `docs/SRIM_MODE.md` §7. They belong to the trim2013 transport. The
   crystal transport takes none of TRIM's bugs: the engine runs it with
   `TrimConfig::effective_bugs()` (all corrected) and
   `effective_single_precision_tables()` (false), so code shared between
   the transports must read `config_.bugs` / `cfg_.bugs` of the engine,
   never the scenario's.
3. **Determinism.** Results must not depend on the thread count. Each ion
   uses `make_stream(seed, ion_index, k)` (k = 0 ion, 1 legacy recoils,
   2+ cascades when `cascade_random_replay` is fixed). Blocks of 16 ions are
   merged in order. Do not share mutable state between blocks except through
   `TrimTallies::add` in `run()`.
4. **Units.** `Mover` and the model use Boost.Units quantities; some hot
   loops work in Å and eV as doubles (`in_units(...)` at the boundary, the
   unit in the variable name or comment). Keep that pattern; do not
   introduce unlabelled doubles.
5. **Tallies.** A new tally goes into `TrimTallies` (`trim.h`), its
   `resize`, `without_lists` and `add`, the engine, `make_dataset` (a table
   or scalar), `trim_output.cpp` if TRIM writes it, and `docs/OUTPUT.md`.
   The damage-model tallies of the `dpa` table (not TRIM) are booked only by
   `Transporter::book_pka`, at the three places a primary gets its energy
   (`kp_2013`, the ion loop's full-cascade branch, `run_recoil_event_2013`);
   it reads, never changes, TRIM's state.
6. **Evidence.** Physics changes need a TRIM.exe routine or a significant
   comparison (`tools/trim_compare`); see AGENTS.md rule 4. After any change
   to the transport, run at least `common-test --gtest_filter='Trim*'` and,
   for physics, the full comparison (`docs/DEVELOPMENT.md` §6).
7. **Speed without changing results.** An optimisation of the default
   (exact) path must give bit-identical output: same expressions in the
   same order (no reassociation, no reciprocal instead of a division, no
   other math function), checked with `tools/perf/bench.py run` on the old
   and new executable and `bench.py compare` (all files identical). Anything
   that changes the last bits — tables, polynomial approximations,
   reciprocals — goes behind `cfg_.fast_math` (`docs/SRIM_MODE.md` §8.2)
   through the Transporter helpers `scatter`, `azimuth` and `turn`. The
   trajectory/MD bookkeeping (`Mover::time`, `parent_*`, `md_mask`) must
   never draw a random number or touch a tally, except the handoff mode's
   `Mover::handed` (record mode is tested to leave every tally identical). The
   transport loops never log, allocate per step or take locks; workers
   update their crash breadcrumb (`diag::Activity`) once per ion.

## Gotchas

- `TrimConfig::straggling_distribution` applies to the refined transport;
  trim2013 uses TRIM's uniform deviate (R − ½)·Ω·√L (with √12 folded into
  the constants).
- TRIM's stop energy is always 1 eV (`final_energy`), not min(Ed, Es, Eb).
- Recoil electronic stopping uses the recoil species' stopping; the
  straggling constants are the ion's (as TRIM).
- Depth tallies are per bin of the depth window; totals (TDATA header,
  summary) include events outside the window.
- Replacements need the struck atom of the slot below to be displaced and
  of the same species (`recoil_stop_2013`).
- Not reproduced on purpose: TRIM's single-precision recoil slots and its
  recoil-stack limit (`docs/TRIM2013_ALGORITHM.md`). The open Xe-gas
  difference may come from the first.
