# AGENTS.md — working on AMCSET

This is the entry point for coding agents (Claude Code reads it through
`CLAUDE.md`) and for human contributors. Read it fully before changing code;
it is short on purpose and links to the details. Directory-specific rules
live in `CLAUDE.md` files next to the code (`common/`, `common/src/trim/`,
`gui/`, `cli/`, `test/`, `tools/`); they are loaded when you work there.

State of this file: AMCSET 1.0.0, 2026-11-26.

---

## 0. The ten rules

1. **Units are a type.** Physical values are Boost.Units quantities
   (`common/include/amcset_units.h`). Build them as `40.0 * electron_volt`,
   read them with `in_units(q, unit)`. No bare doubles for dimensional
   values, no hidden unit conversions. Only reduced BCA quantities (ε,
   b = p/a, sin²(θ/2)) and formulas stated for specific units at a clearly
   commented boundary use plain doubles.
2. **SRIM mode reproduces TRIM.exe, bugs included.** The `trim2013`
   transport (`common/src/trim/trim_engine.cpp`) follows SRIM-2013's TRIM.exe
   step by step, as recovered by disassembly (`docs/TRIM2013_ALGORITHM.md`).
   Do not "fix" or "improve" its behaviour in place. A TRIM bug gets a
   switch in `trim::TrimBugs` (default: reproduce) — see
   `common/src/trim/CLAUDE.md`.
3. **Numerical behaviour is a contract.** A refactor must not change
   results. SRIM-mode runs are deterministic for a given seed and
   independent of the thread count, and `trim_test.cpp` checks exact
   tallies. Native-engine refactors are checked with the determinism harness
   (`docs/TESTING.md`). An intended physics change is made alone, with
   evidence, a test and a `CHANGELOG.md` entry.
4. **Evidence for physics.** A change to SRIM-mode physics needs either the
   TRIM.exe routine it reproduces (address + reasoning) or a statistically
   significant comparison with TRIM reference runs (`tools/trim_compare`,
   `docs/VALIDATION.md` §2). Fe → Fe is the priority case; never tune to a
   single case.
5. **Every change builds and passes the tests** (247 CTest tests on Windows, 246 on Linux: the
   GoogleTest suites and the CLI end-to-end runs). Add a focused test for
   every bug fix or feature.
6. **One scenario model.** `config::RunConfig` (`run_config.h`) describes a
   calculation for both engines and both front ends. A new setting is added
   there once (struct field, parse, serialise, validate) and then exposed in
   the CLI and GUI.
7. **License header** (exact text) on every `.h`, `.cpp`, `.py` outside
   `build/` and `generated/`; `python3 scripts/check_headers.py` must pass
   (CI runs it).
8. **No TRIM code in the repository.** Reverse engineering TRIM.exe was done
   to understand its behaviour; AMCSET re-implements that behaviour. Never
   add decompiled code, disassembly listings or SRIM binaries. The only
   SRIM-derived data in the repository is
   `data/srim2013/srim2013_stopping.amcset` (SRIM's non-commercial terms,
   see its README).
9. **Keep docs true.** When behaviour, a key, a file format, a command or a
   number quoted in the docs changes, update the owning doc (§7) in the same
   change.
10. **Do not commit or push unless asked.** Work in progress usually sits
    uncommitted on `develop` (§8).

---

## 1. What AMCSET is

AMCSET (Texas A&M University, GPL-3.0-or-later) is a C++20 Monte Carlo
simulator of ion (and electron) transport and radiation damage in layered
amorphous targets.

| Engine | Code | Status |
|---|---|---|
| **SRIM mode** (default) | `trim::TrimSimulation`, `common/src/trim/` | Reproduces SRIM-2013 TRIM statistically: 1031 of 1032 compared quantities within 3σ over 161 TRIM runs; one open difference (replacements in Xe gas, `docs/VALIDATION.md` §2.2). Bundled SRIM-2013 stopping tables, TRIM.IN/TRIM.DAT input, every TRIM output file, CSV/JSON/Parquet output, 11 switchable TRIM bugs, about 50× less machine time than TRIM.exe. |
| **Native engine** | `common::Simulation`, `common/src/{ion,electron,simulation,...}.cpp` | The original AMCSET binary-collision code (ions and electrons). Physics bugs P1–P6 fixed in 0.1.0; transport results **not validated**. |

Both implement `common::SimulationEngine` (`simulation_engine.h`): `run()`
on a worker thread; progress, live tallies and trajectories polled from
another thread; cancel and (SRIM mode) pause.

Front ends: **`amcset`** (CLI, `cli/`, `docs/CLI.md`) and **`amcset-gui`**
(Qt 6 Widgets with GPU path views through QRhi, `gui/`, `docs/GUI.md`).
Both read and write JSON scenarios and TRIM.IN.

---

## 2. Repository map

```
AGENTS.md, CLAUDE.md        This guide; CLAUDE.md imports it and adds Claude Code notes
README.md, CHANGELOG.md     User overview; changes per version
CMakeLists.txt              Top level: options AMCSET_BUILD_GUI/CLI/TESTS, ENABLE_ASAN,
                            compiler_flags target, find_package, config.h generation
CMakePresets.json           "release" preset (inherits Conan's conan-release)
conanfile.py                Pinned deps: boost 1.89, gtest 1.15, abseil 20260526.0,
                            qt 6.11.1 (shared, Vulkan), gflags, libiconv
qbuild.sh / qbuild.ps1      Conan build helpers (-c install, -b build, -t test, -r run)
Dockerfile, run-dev.sh      Ubuntu 24.04 development container (GCC 14)
common/                     amcset-common (static library): everything except the front ends
  include/amcset_units.h      Quantity types, unit constants, physical constants, in_units
  include/simulation_engine.h Engine interface (both engines, both front ends)
  include/run_config.h        RunConfig (scenario JSON) <-> engines; unit parsing; make_engine
  include/physics/            stopping.h (models, SRIM tables, StoppingCurve),
                              scattering.h (screening, magic formula, kinematics)
  include/trim/               trim.h (TrimConfig, TrimBugs, TrimTallies, TrimSimulation),
                              random.h (L'Ecuyer, per-ion streams), trim_in.h, trim_dat.h,
                              trim_tables.h (tallies -> output::Dataset), single_precision.h
  include/output/             table.h (Dataset/Table/Column, OutputSpec), units.h, writers.h
  include/crystal/            Crystals (not TRIM; followed by the crystal transport only): geometry.h, unit_cell.h
                              (lattice, sites, built-in structures), orientation.h (Miller, Euler,
                              tilt/twist), crystal.h (Debye amplitudes), cell_io.h (CIF, POSCAR,
                              LAMMPS data, extended XYZ), lattice_search.h (collision partners
                              on lattice sites), materials.h (common crystals: Gemmell 1974),
                              grain_map.h (polycrystals: grains, generated Voronoi), grain_io.h
                              (grain files, Neper, EBSD)
  include/data_files.h        Locates data/ (the bundled SRIM tables) at run time
  include/element_data.h      Element table (TRIM's ATOMDATA values + NIST I-values)
  include/isotope_data.h      NIST isotopes (generated database)
  include/amcset_common.h     Native engine: Particle/Ion/Electron/Layer/Volume/Simulation
  include/tallies.h           Welford, Histogram, Tallies (live statistics of both engines)
  src/trim/                   SRIM mode: trim_engine.cpp (model + transport + driver),
                              trim_output.cpp (TRIM files, summary.json), trim_tables.cpp,
                              trim_in.cpp, trim_dat.cpp, single_precision.cpp
  src/physics/                stopping.cpp, scattering.cpp (shared physics)
  src/crystal/                unit_cell.cpp, orientation.cpp, crystal.cpp, cell_io.cpp, lattice_search.cpp,
                              materials.cpp, grain_map.cpp, grain_io.cpp
  src/output/                 table.cpp (shaping), units.cpp, writers.cpp (CSV/JSON), parquet.cpp
  src/run_config.cpp          JSON <-> RunConfig, overrides, make_engine (incl. native adapter)
  src/*.cpp                   Native engine (ion.cpp holds all ion models)
  assets/, generated/, scripts/  Element/isotope sources, generated headers, generators
data/srim2013/              Bundled SRIM-2013 stopping tables (42 MB) + README (terms)
cli/                        amcset-cli-lib (cli.cpp: options -> RunConfig, stopping tables)
                            and amcset (main.cpp: commands, help text)
gui/                        amcset-gui: src/{app,model,run,setup,views}, shaders/ (+ baked .qsb)
test/common/, test/cli/     GoogleTest suites (common-test, cli-test); CLI e2e tests are
                            declared in test/CMakeLists.txt
tools/trim_compare/         TRIM output readers, AMCSET-vs-TRIM statistics (batch, analyze)
tools/damage_models/        Damage functions per PKA and model; arc-dpa vs published SRIM post-processing
tools/validation/           SRIM-2013 table scan over every ion; stopping against the IAEA database and
                            NIST (stopping_reference.py); AMCSET against measurements from the literature
                            (literature_cases.py with range_data.py, sputter_data.py: Ni peaks, 350 ranges,
                            sputtering, published SRIM damage, As channeling, recoil ionization); docs/VALIDATION.md;
                            make_figures.py draws the figures of docs/validation/ (needs matplotlib)
tools/crystal/              Checks of the crystal transport (channeling trends, normalisation) and its
                            comparison with measurements (calibration.py); polycrystal_checks.py;
                            range_study.py (why ranges in crystals are longer)
tools/srim/                 SRIM-2013 table extraction, TRIM reference cases and runs (Windows)
tools/windows/              build_and_test.bat, compare_with_trim.bat, benchmark_vs_trim.bat,
                            gui_smoke_test.bat; make_msix.ps1, sign.ps1 (releases), setup_runner.ps1 (CI)
tools/linux/                make_appimage.sh
packaging/                  Install rules and CPack (CMakeLists.txt); windows/ (manifest, version
                            resource, MSIX manifest, WiX patch); linux/ (desktop entry, AppStream
                            metadata, man pages, Debian copyright); flatpak/ (unbuilt manifest)
resources/icons/            The logo: SVG, PNGs, ICO, MSIX tiles (drawn by scripts/make_icons.py)
tools/determinism_harness/  Bit-for-bit regression harness (native engine)
tools/perf/                 bench.py (speed and bit-identity of a SRIM-mode change), fast_vs_exact.py
                            (physics.fast_math against exact), fuzz_inputs.py, winsample.cpp (profiler)
scripts/                    check_headers.py (CI), add_headers.sh, make_icons.py, cascade animation (py/jl)
docs/                       See docs/README.md (index and reading order)
.github/workflows/          check-headers.yml (Ubuntu), tests.yml (self-hosted Windows),
                            release.yml (packages on GitHub's runners; draft release for a v* tag)
```

Generated or local (all in `.gitignore`): `build*/`, `output/`,
`common/include/config.h` (from `config.h.in`), `compile_commands.json`,
`ConanPresets.json`. `common/generated/*.h` are generated but committed;
never edit them by hand; nor `resources/icons/` (redraw with
`scripts/make_icons.py`).

---

## 3. Build and test

Details and troubleshooting: `docs/DEVELOPMENT.md`. Test catalogue:
`docs/TESTING.md`.

### Linux without Qt (fastest loop: core, CLI and tests)

Needs a C++20 compiler (GCC ≥ 13), CMake ≥ 3.23, Ninja, Boost headers
≥ 1.83, abseil with `absl::log` (≥ 20230125; Ubuntu 24.04's is too old —
build `20250127.1` from source with `-DCMAKE_CXX_STANDARD=20`) and
GoogleTest.

```bash
cmake -S . -B build-noqt -G Ninja -DCMAKE_BUILD_TYPE=Release -DAMCSET_BUILD_GUI=OFF \
      [-DCMAKE_PREFIX_PATH=/path/to/abseil-install]
cmake --build build-noqt -j
ctest --test-dir build-noqt --output-on-failure        # 246 tests, 1-2 min
build-noqt/common-test --gtest_filter='Trim*'         # one suite
```

CTest runs the tests with the repository root as working directory.

### Windows (the main workstation; MSVC + Conan)

```bat
set AMCSET_NO_PAUSE=1
tools\windows\build_and_test.bat            :: configure, build (ninja -k 0), ctest
tools\windows\build_and_test.bat --conan    :: reinstall the Conan dependencies first
```

Log: `build\build_and_test.log`; the last line of
`build\build_and_test.status` is `SUCCESS`, `BUILD FAILED` or `TESTS
FAILED`. `AMCSET_NO_PAUSE=1` skips the 2-minute pause at the end (needed
for non-interactive use). Outputs: `build\Release\amcset.exe`,
`amcset-gui.exe`, `common-test.exe`, `cli-test.exe`. Run them after `call
build\Release\generators\conanrun.bat` (Qt and abseil DLLs on PATH).
Linking `amcset-gui.exe` fails with LNK1104 while the GUI is running.

### Linux with Qt (Conan)

```bash
conan profile detect --force
./qbuild.sh -c -b -t     # conan install (builds Qt the first time: hours), build, test
./qbuild.sh -r           # run the GUI
```

Release flags: `-O3 -march=native -flto` (GCC) or `/O2 /GL` (MSVC), with
`ABSL_MIN_LOG_LEVEL=2`. Clang is not supported (the isotope table is empty
under clang).

### Checks after a change

| Changed | Run at least |
|---|---|
| anything | build + `ctest`; `python3 scripts/check_headers.py` |
| `common/src/trim/`, `common/src/physics/` | `common-test --gtest_filter='Trim*:Stopping*:Scattering*'`; for physics changes the TRIM comparison (§6) |
| native engine (`common/src/*.cpp`) | determinism harness before and after (`docs/TESTING.md`) |
| `run_config.*`, `cli/` | `cli-test`; `ctest -R cli` |
| output tables or writers | `common-test --gtest_filter='Output*:TrimTables*'` |
| `gui/` | build with the GUI; `tools\windows\gui_smoke_test.bat` (screenshots of every view) |

---

## 4. Architecture in one page

Full description: `docs/ARCHITECTURE.md`.

```
scenario.json / TRIM.IN (+TRIM.DAT) / CLI options / GUI forms
        │  config::parse_config, trim::read_trim_in, cli::apply_options, gui::Scenario
        ▼
config::RunConfig { mode, trim::TrimConfig trim, NativeOptions native }
        │  config::make_engine
        ▼
common::SimulationEngine ─┬ trim::TrimSimulation  (SRIM mode)
                          └ common::Simulation    (native engine)
        │  run() on a worker thread; front ends poll completed_count(),
        │  snapshot_tallies(), try_pop_trajectory(); request_cancel(), set_paused()
        ▼
SRIM mode: TrimTallies ── trim::make_dataset ──► output::Dataset (typed tables + scalars)
        │                                             │ output::apply(table, OutputSpec)
        ▼                                             ▼
TRIM text files (trim_output.cpp)         CSV / results.json / Parquet (common/src/output/)
```

SRIM mode: `TrimSimulation` builds a `TargetModel` (per layer and element
constants, per species × layer `StoppingCurve`, collision pairs, TRIM's
layer averages). `run()` processes ions in fixed blocks of 16 on
`std::thread`s; each block has its own `Transporter` and `TrimTallies`;
blocks are merged in block order, so results do not depend on the thread
count. Every ion draws from its own L'Ecuyer stream `make_stream(seed, ion,
k)`. `Transporter::run_ion_2013` is TRIM's ion loop, `cascade_2013` TRIM's
recoil stack.

GUI: `MainWindow` owns `Scenario` (the edited RunConfig with a QUndoStack),
`RunController` (engine on a `std::thread`, a GUI-thread timer polls it and
builds a live `Dataset` from tally snapshots; pause, stop, continue) and
`TrackStore` (path segments for the GPU views). Setup panels bind widgets
to the scenario through `FormBinder`.

---

## 5. Conventions

**Units and constants.** See rule 1. Number densities are atoms per volume
(m⁻³ internally): N = ρ / mean atom mass; there is no Avogadro factor
anywhere. Physical constants come from `amcset_units.h` (Boost.Units, CODATA
2006). Where SRIM mode reproduces TRIM it uses TRIM's own constants
(e² = 14.4 eV·Å, a₀ = 0.5292 Å, 0.8853, 931494 keV per amu, …) and says so
in a comment.

**Coordinates.** SRIM mode follows TRIM: **x is depth**, y and z lateral; a
tilted beam tilts in the x–y plane. The native engine and the
`SimulationEngine` trajectories use **+z = depth** (SRIM mode converts when
it publishes trajectories). Native `Layer::depth` is the cumulative depth of
the layer's bottom; `trim::TargetLayer::width` is a thickness.

**Citing TRIM.exe.** Each block of the trim2013 transport names the TRIM.exe
routine it reproduces: an address (`0x45d4e0`) or TRIM's own location label
(`RECOILLOOP`). Only cite addresses listed in `docs/TRIM2013_ALGORITHM.md`
or verified in TRIM.exe; never invent one.

**Errors.** Configuration errors throw `std::invalid_argument` with a
message naming the offending key (`physics.trim_bugs.foo: ...`); I/O and
data problems throw `std::runtime_error`. The native engine wraps
exceptions with `rethrow(EXCEPTION_MESSAGE(""))` inside function-try-blocks
(only from inside a `catch`). Never let an exception escape a Qt slot.

**Messages to users.** Library code does not print. Engines collect warnings
(missing bundled tables, replaced zero energies, …) and the front ends show
them.

**Logging.** abseil `LOG(INFO|WARNING|ERROR)` / `VLOG(n)`; INFO and WARNING
are compiled out in Release. No side effects (e.g. random draws) inside log
statements.

**Style.** C++20; 2-space indentation; `snake_case` for functions,
variables and files (also in `gui/`); classes in `PascalCase`; data members
end with `_`. Doxygen comments (`/*! \brief ... */`, `//!`). Comments state
facts and reasons. Warning flags `-Wall -Wextra -Wshadow -Wformat=2
-Wunused` (GCC) and `/W3` (MSVC); keep the build warning-free.

**Files.** LF line endings except `.bat` and `.ps1` (CRLF); keep a file's
existing endings when editing. UTF-8. The license header is the 16-line
`// Copyright 2025, Texas A&M University ...` block (`#` comments in
Python; a shebang goes on line 1 above it); `scripts/add_headers.sh` adds
it.

**Tests.** GoogleTest; suites named by area (`Trim*`, `Stopping`,
`Scattering`, `Output*`, `RunConfig*`, `Cli`, `CommonTest` for the native
engine). `common_test.cpp` uses `#define private public` to read private
*data members* of the native engine only — never to call private functions
(MSVC encodes access in mangled names, so such calls do not link).
SRIM-mode tests may compare exact tallies because its runs are
deterministic.

---

## 6. How to make a change

1. Read the directory's `CLAUDE.md` and the doc that owns the topic (§7).
2. For physics, find the evidence first (TRIM routine, published data, or
   a comparison) and state it in the code comment and the change
   description.
3. Make the smallest change; keep refactors separate from behaviour
   changes.
4. Add or adjust tests. When an exact-tally expectation in `trim_test.cpp`
   must change, change it deliberately and say why.
5. Build, run `ctest`, run `scripts/check_headers.py`.
6. For SRIM-mode physics, run the TRIM comparison and confirm nothing moved
   outside the statistics (`tools\windows\compare_with_trim.bat` on the
   workstation, about an hour; `docs/DEVELOPMENT.md` §6). Fe → Fe first.
7. Update the docs (§7) and `CHANGELOG.md`.

Recipes for the common changes — new scenario key, new TRIM bug switch, new
output table or column, new GUI field or plot, new stopping model, new CLI
option, regenerating the stopping tables, adding a TRIM reference case — are
in **`docs/HOWTO.md`**.

### Definition of done

- Builds with MSVC (Windows) and GCC (Linux) without new warnings.
- All CTest tests pass; the new behaviour has a test.
- The header check passes; no generated or machine-local files are added.
- Results are unchanged unless intended; intended changes are measured and
  documented (numbers in `docs/SRIM_MODE.md` and `CHANGELOG.md`).
- Docs updated per §7.

---

## 7. Which document owns what

| Topic | Document |
|---|---|
| Doc index, reading order | `docs/README.md` |
| Architecture, data flow, threading, module boundaries | `docs/ARCHITECTURE.md` |
| Building, environments, workstation, validation runs, benchmark | `docs/DEVELOPMENT.md` |
| Recipes for common changes | `docs/HOWTO.md` |
| Tests, suites, determinism harness, optional SRIM-data tests | `docs/TESTING.md` |
| SRIM-mode physics, parameters, TRIM bugs, performance | `docs/SRIM_MODE.md` |
| Validation: TRIM comparison, measurements, MD, deviations and optional corrections | `docs/VALIDATION.md` (numbers, parameters); `docs/validation/` (the review with figures; redraw with `tools/validation/make_figures.py` after a rerun) |
| TRIM.exe's algorithm with routine addresses | `docs/TRIM2013_ALGORITHM.md` |
| CLI commands and the full scenario schema | `docs/CLI.md` |
| Output files, tables, OutputSpec | `docs/OUTPUT.md` |
| GUI user guide and code layout | `docs/GUI.md` |
| Open bugs and limitations | `docs/KNOWN_ISSUES.md` |
| Improvement backlog | `docs/ROADMAP.md` |
| Crystalline targets: design, phases, what is built | `docs/CRYSTAL_DESIGN.md` (user guide: `docs/CLI.md` "Crystals") |
| Terms (Ed, Eb, Es, KP, PKA, monolayer, ε, ...) | `docs/GLOSSARY.md` |
| Version history | `CHANGELOG.md` |
| Installing the packages | `docs/INSTALL.md` |
| Packages, release procedure, code signing | `docs/RELEASING.md` |
| Windows and Linux platform requirements and how they are met | `docs/PLATFORM.md` |
| Bundled stopping tables and their terms | `data/srim2013/README.md` |
| SRIM data extraction, TRIM reference runs | `tools/srim/README.md` |

---

## 8. Current state (2026-11-25)

- **Git:** branch `develop` (remote `git@github.com:Shao-RMSGAL/AMCSET.git`).
  Everything up to the 1.0 release files and the validation review is
  committed and pushed (2026-11-26). Commit only when asked; then group the
  commits and check `git status` for build outputs first. The data file is under
  GitHub's 50 MB warning size; use Git LFS if it ever grows past it.
- **Release 1.0:** version 1.0.0, packages, release workflow, signing script
  and platform documents are in place (`docs/RELEASING.md`,
  `docs/PLATFORM.md`). 1.0 is published unsigned (docs/INSTALL.md says
  what users do about Windows' warnings); signing through the SignPath
  Foundation is the maintainer's choice for later (RELEASING.md §3.2).
  `.github/workflows/release.yml` built every package of `v1.0.0`
  (2026-11-26: MSI, ZIP, AppImage, CLI tarball, .deb, .rpm, SHA256SUMS.txt)
  into a draft release. Not yet done: the MSIX in the workflow (needs
  signing), the App Certification Kit, AppStream screenshots, the Flatpak
  build.
- **Last full validation:** `build/trim_validation/trim3/analysis.txt` on the
  workstation (75 units, 161 TRIM runs, 704 AMCSET runs), summarised in
  `docs/VALIDATION.md` §2.2.
- **Open physics question:** replacements in Xe gas (+4.6 %, z = 5.0),
  `docs/VALIDATION.md` §2.2; leading suspect: TRIM's single-precision recoil
  slots.
- **Crystals:** phases 1 and 2 of `docs/CRYSTAL_DESIGN.md` are in: unit
  cells, orientation, cell files, `layers.N.crystal`, `amcset cell`, and the
  crystal transport (`physics.transport: crystal`, after Robinson and
  Torrens 1974) for single crystals; also the common crystals of Gemmell
  1974, the GUI's Crystal tab (unit-cell editor with a 3D view) and a first
  comparison with measurements (`tools/crystal/calibration.py`: a He
  channeling dip in Ge, B profiles in Si). Those data bracket the division
  of the electronic stopping (0.25 to 1) without fixing it: the default 0.5
  is not a fit. Ranges in crystals exceed amorphous ones by channeling
  (`tools/crystal/range_study.py`, docs/CLI.md "Why ranges in a crystal are
  longer"). Polycrystals (phase 3, **on hold** by the maintainer's choice:
  single crystals and amorphous targets first) run from a scenario
  (`layers.N.crystal.polycrystal`: generated grains, grain files, Neper,
  EBSD) but have no GUI yet, and the Neper and EBSD readers are untested on
  real files: `docs/CRYSTAL_DESIGN.md` §13 lists what is left. Then: dose
  dependence (phase 5), a fit and wider validation (phase 6: low-dose
  profiles with stated tilt and twist, W and Fe, Cu/Au RBS).
- **Not done yet:** a controlled benchmark (`tools\windows\benchmark_vs_trim.bat`
  on an idle machine), validation of the native engine, a Linux CI job.

Backlog and priorities: `docs/ROADMAP.md`. Known problems:
`docs/KNOWN_ISSUES.md`.
