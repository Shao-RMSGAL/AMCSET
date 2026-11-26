# Architecture

How AMCSET is put together: modules, data flow, threading and the
boundaries that keep it maintainable. For the physics see
`docs/SRIM_MODE.md`; for TRIM.exe's algorithm `docs/TRIM2013_ALGORITHM.md`.

## 1. Modules and dependencies

```
            ┌──────────────┐      ┌──────────────┐
            │  amcset-gui  │      │    amcset    │   front ends
            │  gui/ (Qt 6) │      │ cli/main.cpp │
            └──────┬───────┘      └──────┬───────┘
                   │              ┌──────┴───────┐
                   │              │amcset-cli-lib│   option parsing, stopping command
                   │              └──────┬───────┘
                   └───────────┬─────────┘
                        ┌──────┴────────┐
                        │ amcset-common │            static library, no Qt
                        └──────┬────────┘
     ┌──────────┬──────────────┼─────────────┬───────────────┬──────────────┐
  config::   trim::        physics::      output::        common::        data::
  run_config SRIM mode     stopping,      tables, units,  native engine,  bundled
  (scenario) engine, I/O   scattering     writers         units, tallies  data files
```

| Namespace | Directory | Responsibility | May use |
|---|---|---|---|
| `amcset::common` | `common/include/*.h`, `common/src/*.cpp` | units (`amcset_units.h`), engine interface, live `Tallies`, element/isotope data, native engine | Boost, abseil |
| `amcset::physics` | `common/*/physics/` | electronic stopping (models, SRIM tables, `StoppingCurve`), straggling, screening, magic formula, kinematics | `common` |
| `amcset::trim` | `common/*/trim/` | SRIM-mode engine, TRIM.IN/TRIM.DAT, TRIM files, tallies → dataset | `common`, `physics`, `output` |
| `amcset::crystal` | `common/*/crystal/` | unit cells (lattice, sites with occupancies, built-in structures), orientation (Miller indices, Euler angles, tilt/twist), Debye amplitudes, the common crystals of Gemmell 1974 (`materials.h`), cell files (CIF, POSCAR, LAMMPS data, extended XYZ), and the collision-partner search of a lattice (`LatticeSearch`: the atoms ahead of a mover, thermal displacements and occupancy from a hash of the site). Geometry in Angstrom doubles. Used by the crystal transport of `trim::` (`docs/CRYSTAL_DESIGN.md`) | `common` |
| `amcset::output` | `common/*/output/` | typed tables (`Dataset`), output units, `OutputSpec` shaping, CSV/JSON/Parquet writers | `common` |
| `amcset::config` | `run_config.*` | `RunConfig` (scenario), JSON, unit parsing, overrides, `make_engine` | all of the above |
| `amcset::data` | `data_files.*` | locate `data/` (bundled stopping tables) | std only |
| `amcset::diag` | `diagnostics.*` | diagnostic log, crash reports, per-thread breadcrumbs | abseil (sink, stack traces), DbgHelp on Windows |
| `amcset::cli` | `cli/` | command line → `RunConfig`; `stopping` and `cell` commands | `config`, `physics`, `crystal` |
| `amcset::gui` | `gui/src/` | Qt application | everything in `common` (never the reverse) |

Rules: `common` never includes Qt or CLI headers. Physics used by both
engines lives in `physics`. Front ends never compute physics or format
result files themselves.

Third-party: Boost (Units, JSON, Asio thread pool in the native engine), abseil (logging,
stack traces and symbols on Linux),
GoogleTest, Qt 6.11 (Widgets, Gui, QRhi; GUI only). Parquet is written by
AMCSET's own minimal writer (`common/src/output/parquet.cpp`, no Arrow).

## 2. The scenario pipeline

```
scenario.json ──config::load_config──┐
TRIM.IN (+TRIM.DAT) ─trim::read_trim_in─┤
CLI options ──cli::apply_options─────┤──► config::RunConfig ──config::make_engine──► SimulationEngine
GUI forms ──gui::Scenario::modify────┘          │
                                                └─ config::to_json / save_config / trim::write_trim_in (+ write_trim_dat)
```

- `config::RunConfig` = `mode` (srim | native) + `trim::TrimConfig trim` (ion,
  layers with per-element Ed/Eb/Es, counts, seed, threads, stopping,
  physics options, bug switches, output options — used by **both** engines)
  + `config::NativeOptions native` (native-only models and files).
- JSON keys are documented in `docs/CLI.md`. Values with dimensions are
  strings with units (`"2 MeV"`, `"1.2 um"`); `parse_config` rejects unknown
  keys and names the offending key in every error. `to_json` writes the
  canonical form; `apply_override("a.b.0.c=value")` edits the JSON tree and
  re-parses, so every key is overridable without extra code.
- `make_engine` validates and builds a `trim::TrimSimulation` (SRIM mode) or
  converts to the native `Simulation::Settings` + `Volume`.
- TRIM.IN parsing (`trim_in.cpp`) produces a `TrimConfig` plus
  `TrimInExtras` (fields without a `TrimConfig` equivalent: cascade code,
  plot type, disk-file flags, …) so that a TRIM.IN can be written back;
  cascade codes 4–7 load the TRIM.DAT next to it.

## 3. SRIM mode engine (`trim::TrimSimulation`)

### Construction

`TrimSimulation(TrimConfig)` validates the configuration and builds a
`TargetModel` (internal, `trim_engine.cpp`):

- TRIM's zero-energy defaults (Ed 20, Es 4, Eb 2 eV; skipping atom 1 unless
  that bug is fixed);
- per layer: elements with fractions, Ed/Eb/Es, Robinson KP constants,
  density, monolayer N^(-1/3), and TRIM's layer averages (`Trim2013Layer`:
  integer-accumulated Z, screening length, F_I, EPSDG, p_max0, KP and
  straggling constants) — TRIMSTART;
- species (the ion and each distinct target atom) × layer:
  `physics::StoppingCurve` (from the bundled SRIM tables by default),
  straggling table, collision pairs (universal screening; SRIM pair
  potentials only as an AMCSET option);
- the particle source (beam or TRIM.DAT list), the depth window and bins;
- with `fast_math`, `physics::MagicTable::universal()` (the tabulated magic
  formula, 5 MB, built once per process on all cores in ~0.1 s).

The stopping tables come from `SrimStoppingDatabase::load`: the 42 MB
bundled file is read in one call and only indexed (one entry per "table"
line); a table is parsed the first time `find()` asks for it, under the
database's mutex. A calculation needs a few of the 16 928 tables, so start-up
went from ~0.8 s to ~0.05 s. Loaded databases are kept per process (path,
size, time).

### Run

```
run()
  threads = config.threads or hardware_concurrency (fewer if threads cannot be created)
  diag::set_context(engine, seed, ions, ...)
  blocks of 16 ions, handed out by an atomic counter
  worker (diag::Activity breadcrumb): for each block
      Transporter tr(config, model, block TrimTallies)
      for each ion: breadcrumb = ion; tr.run_ion(i)   ── trim2013: run_ion_2013
      pending[block] = tallies; merge every finished block in block order
  write_outputs(output_directory) if set (files written concurrently)
```

- **Determinism.** Ion i uses `make_stream(seed, i, k)` (SplitMix64-seeded
  L'Ecuyer). Blocks are merged in order, so every floating-point sum is the
  same for any number of threads; `continue_from` + `first_ion` make
  "N ions, then M more" identical to one N+M run.
- **Live data.** `snapshot()` copies the merged `TrimTallies` (optionally
  without the long lists) for the GUI; `snapshot_tallies()` derives the
  coarse `common::Tallies` shared with the native engine from that snapshot
  (on TRIM's depth bins; nothing is filled per collision for it);
  `try_pop_trajectory()` returns collision points of the first
  `trajectory_ions` ions (converted to +z = depth).
- **Failures.** An exception in a worker names the ion, the seed and the
  block start to rerun it (`first_ion`), cancels the other workers and is
  rethrown by `run()` after every thread has joined. Crashes (no exception)
  are reported by `diag` (§9).
- **Pause/cancel.** Workers wait on a condition variable while paused
  (checked between blocks) and stop taking ions once `cancel_` is set
  (checked between ions); a cancelled run still merges the finished blocks
  and writes its outputs.
- **Transporter** (`trim_engine.cpp`): `run_ion_2013` (ion loop: flight
  path, electronic loss, collision, tallies, KP or cascade),
  `cascade_2013` (TRIM's recoil slot stack), `recoil_stop_2013`,
  `kp_2013`, `electronic_loss_2013`, `ion_flight_path_2013`,
  `interface_2013`, `edge_2013`, `pick_element`; the refined transport is
  `follow`/`cascade`. Map to TRIM.exe routines:
  `common/src/trim/CLAUDE.md`.
- **Tallies** (`TrimTallies`, `trim.h`): depth histograms (range, vacancies
  per element, replacements, interstitials, ionization, phonons, energy to
  recoils, …), moments, counts, surface-energy histograms, recoil spectra,
  exit lists, ion positions, EXYZ records, cascade records. `add()` merges.

### Outputs

```
TrimTallies ──trim::make_dataset──► output::Dataset { metadata, scalars, tables[] }
     │                                   │
     │ trim_output.cpp                   │ output::apply(table, OutputSpec)  (columns, units,
     ▼                                   ▼                                     rows, number format)
TRIM text files (+summary.json)     write_csv / write_json (results.json) / write_parquet
```

`write_outputs` writes every file as its own task (`common::run_tasks`,
`parallel_tasks.h`): the TRIM files, `summary.json`, and the dataset files
(one task per table, one for `results.json`). The files are independent and
byte-identical to sequential writing; the long particle lists (SPUTTER.txt,
TRIMOUT.txt, RANGE_3D.txt) no longer add up.

The `Dataset` is the single source for every non-TRIM format and for the
GUI's plots, tables and statistics. Columns carry a `Dimension` and a unit
(`output/units.h`), so unit conversion and labels are generic. TRIM-format
files reproduce TRIM's layout (and, by default, its single-precision table
sums).

## 4. Native engine (`common::Simulation`)

The original AMCSET code, kept for electron transport and as a second
engine; not validated (`docs/KNOWN_ISSUES.md` §1).

```
make_engine (native) → Simulation::Settings (3 constructor overloads) + Volume (Layers)
Simulation::run()
  boost::asio::thread_pool(threads); one task per Bombardment
    Bombardment::run_bombardment(): Ion/Electron::fire()  (recursive recoils)
      records into the Bombardment's Tallies and collision points
    merge Tallies (mutex), queue trajectories (mutex)
  write histogram CSVs + tallies.json to the output directory ("output" by default)
```

- All native types are in `amcset_common.h`; the model functions of
  `Simulation` (screening, stopping) are implemented in `ion.cpp`; model
  selection maps the Settings enums to function pointers once.
- Each Bombardment owns an `mt19937` seeded by `EntropySeeder`: runs are not
  reproducible and all Bombardments are created up front (memory grows with
  the ion count).
- Coordinates: +z = depth; directions are unit pure quaternions.
- The GUI shows native results through `gui/src/run/datasets.cpp`
  (`native_dataset` from `snapshot_tallies()`).

## 5. GUI (`gui/src/`)

```
MainWindow
 ├─ Scenario          RunConfig + QUndoStack; modify() = validated undoable edit; load/save/import/export
 ├─ setup panels      IonPanel, TargetPanel (+LayerTableModel, ElementTableModel), CalculationPanel,
 │                    OutputPanel, OutputSpecDialog — widgets bound through FormBinder
 ├─ RunController     start/extend/pause/stop; owns the engine and a std::thread
 │    (auto depth: trim::size_target first on that thread → target_sized → engine)
 │    QTimer (100 ms, GUI thread): completed_count → progress; try_pop_trajectory → TrackStore
 │    (≤ 30 ms per tick); every 400 ms a tally snapshot → make_dataset → dataset_changed;
 │    transmitted_so_far → ions_escaped (MainWindow: log, non-modal retry question)
 ├─ TrackStore        append-only float segments (cap 4 M) for the GPU views
 └─ views             TrackPanel/TrackView (QRhiWidget; Vulkan/D3D/Metal/OpenGL), PlotsPanel/PlotWidget,
                      TablesPanel, StatsPanel — all read the shared output::Dataset
```

- Only the GUI thread touches widgets; the engine thread only runs
  `engine->run()`. Communication is by polling the thread-safe engine
  interface.
- "Continue" (`RunController::extend`) runs another `TrimSimulation` with
  `first_ion` = ions so far and `continue_from(previous tallies)`, possibly
  with another ion/energy/angle on the same target (checked by a target
  signature).
- `results.json` files reopen as a `Dataset` (`read_results_json`), so every
  view works on saved results too.
- Automatic target depth: `RunController::start` validates, then runs
  `trim::size_target` on the worker thread (`sizing()`; `stop()` sets its
  cancel flag); the GUI-thread poll joins it, emits `target_sized` (the
  MainWindow applies the sized thickness to the scenario as an undoable
  edit) and builds the engine on the sized configuration.
- Plausibility warnings (`TrimConfig::plausibility_warnings`) are evaluated
  300 ms after the last scenario change (they compute stopping powers) and
  shown in the status bar and, once each, in the Messages tab.
- Unattended mode (`--run --screenshots DIR --quit --timeout S`) drives the
  smoke test.

## 6. Data files at run time

`data::bundled_srim_stopping_table()` searches `$AMCSET_DATA_DIR/srim2013/`,
`data/srim2013/` next to the executable or up to two folders above,
`share/amcset/srim2013/` next to the executable's folder, then the source
tree (compiled in as `AMCSET_SOURCE_DATA_DIR`). `SrimStoppingDatabase::load`
caches each file for the process (keyed by canonical path, size and mtime);
reading and indexing the 42 MB bundled file takes about 0.05 s, and each
table is parsed when first used (§3). The native engine's
electron tables (`common/assets/*.csv`, installed as `share/amcset/electrons/`)
are found like the stopping tables (`data::find_electron_table`); a table
path given in a scenario is relative to the working directory.

## 7. Performance notes

- SRIM mode is CPU-bound in `Transporter` (flight paths, magic formula,
  cascades); stopping and straggling are table lookups. About half of a
  full-cascade run computes scattering angles (the closest-approach
  iteration's exponentials); `physics.fast_math` replaces that by a table
  (`docs/SRIM_MODE.md` §8.2). Start-up takes about 0.05 s (the stopping
  file is indexed, tables parsed on use); the output files are written
  concurrently, their numbers through `std::to_chars` into large buffers
  (`output/text.h`: printf's exact characters without printf).
- Trajectory and MD event files (`trim_export.h/.cpp`,
  `output/stream.h`): every block collects its points (filtered) and events
  in an `ExportBlock`; `run()` writes the blocks in block order under the
  merge lock, through `StreamTableWriter` (CSV, or Parquet with a row group
  per ~256 k rows), so the files are in ion order for any number of threads
  and are never held in memory. The events are decided in the transport
  (`Transporter::md_birth`, `md_crossing`, `md_ion_entry`): a bit per E_c
  in `Mover::md_mask` marks branches that are MD's; in the handoff mode
  `Mover::handed` makes the cascade skip them. `Mover::time`,
  `parent_track` and `parent_point` are bookkeeping only (no random
  draws, no tallies).
- The hot loops do not allocate, lock or log per collision; the exact path
  must stay bit-identical under optimisation (`common/src/trim/CLAUDE.md`
  rule 7; `tools/perf/bench.py compare`).
- Scaling: blocks of 16 ions keep all cores busy for runs of more than a few
  hundred ions; very small runs are dominated by start-up.
- Memory: tallies are small (100 bins × elements); the exit lists, ion
  positions, EXYZ and COLLISON records grow with the ion count (capped by
  `max_exit_records`); heavy-ion gas cascades can hold 10⁵ pending recoils
  in the slot stack.
- Measured against TRIM.exe: `docs/SRIM_MODE.md` §8.

## 8. Extension points

| To add | Touch | Recipe |
|---|---|---|
| a scenario setting | `trim.h`/`run_config.*` → CLI/GUI → docs | `docs/HOWTO.md` §1 |
| a TRIM bug switch | `TrimBugs` + engine + GUI tooltip + test + docs | `docs/HOWTO.md` §2 |
| a result table/column | `TrimTallies` → engine → `make_dataset` → (TRIM file) → docs | `docs/HOWTO.md` §3 |
| a plot | `PlotsPanel` definitions | `docs/HOWTO.md` §4 |
| a stopping model | `physics::ElectronicStoppingModel` + `layer_stopping` + parser + tests | `docs/HOWTO.md` §5 |
| a CLI option | `cli.h`/`cli.cpp`/`main.cpp` + docs + test | `docs/HOWTO.md` §6 |
| an output format | `output/writers.*`, `output::Format`, `write_outputs` | `docs/HOWTO.md` §7 |

## 9. Errors, diagnostics and crash reports

**Exceptions.** Library code throws `std::invalid_argument` for
configuration errors (the message names the key, the value, what is
accepted and what to use instead) and `std::runtime_error` for I/O, data
and transport failures; it never prints. Every value read from a file
(scenario JSON, TRIM.IN, TRIM.DAT, stopping tables, results.json) is
checked where it is read: finite, whole where a count is expected, in
range before any conversion to an integer type, so nothing wraps around or
indexes past an array. `TrimConfig::validate` checks the physical limits
(docs/CLI.md "Checks and warnings"); `validate(false)` (used while a
scenario is being edited) skips the checks that need a complete scenario
(a particle list, existing files). The transport checks the ion's and
recoils' state once per flight: a NaN would never satisfy the stop
conditions, so it throws with the last flight's values instead of looping.
Worker exceptions carry the ion, seed and rerun options (§3, Failures);
`common::run_tasks` (parallel output files) and the native engine's thread
pool rethrow the first error after every task has finished.

**Diagnostic log** (`amcset::diag`, `diagnostics.h`): `diag::initialize`
(front ends) sets the level, the log file and a ring buffer of the last
lines; abseil `LOG` lines go to the same sink. `AMCSET_LOG(level)` costs one
comparison when the level is off. Engines log run-level events only.

**Crash reports.** `diag::install_crash_handlers()` (first thing in `main`)
installs an unhandled-exception filter on Windows (stack walk through
DbgHelp, minidump) and `sigaction` handlers on an alternate stack on POSIX
(abseil stack traces), plus `std::set_terminate`. The report holds the
stack, the run context (`diag::set_context`), every worker's breadcrumb
(`diag::Activity`: the ion and block being transported) and the last log
lines; the process then exits with code 3.

**GUI.** Qt does not propagate exceptions: none may leave a slot, an event
handler or a virtual override. Every connection goes through
`guarded(context, what, f)` / `guarded_quiet(what, f)` (`app/guard.h`),
which catch everything and report it (a message box, at most one at a
time, or the log); `gui::Application::notify` is the last net for events;
overrides (`data()`, `paintEvent`, `closeEvent`, ...) catch and return a
safe value. `Scenario::modify` validates every edit and returns the message
of a rejected one; `FormBinder::refresh` is `noexcept`. Qt's own messages
go to the diagnostic log through a message handler.
