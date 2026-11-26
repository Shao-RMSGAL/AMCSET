# Roadmap and improvement backlog

Ordered by value within each group. Items reference
`docs/KNOWN_ISSUES.md` (KI). Each item says what "done" means; follow the
rules in `AGENTS.md` (evidence for physics, results unchanged by refactors,
tests, docs). Updated 2026-11-25.

## Done in 0.1.0

- Native engine physics fixes P1–P6 (KI §0), with tests.
- SRIM mode: TRIM-2013 transport rebuilt from TRIM.exe (`trim2013`),
  bundled SRIM-2013 stopping tables (default), TRIM.IN and TRIM.DAT
  (cascade codes 1–7), every TRIM output file plus CSV/JSON/Parquet with
  user-chosen tables, units and formats, seeded and thread-count independent
  results, continue/superpose runs, 11 switchable TRIM.exe bugs.
- Validation against 161 TRIM-2013 runs of 75 inputs (1031 of 1032
  comparisons within 3σ) and run-time comparison (~50× less machine time).
- CLI with JSON scenarios shared with the GUI; new Qt 6 GUI (undoable
  scenario editing, live TRIM plots, tables, statistics, GPU path views).
- Test suite (126 CTest tests), Windows build/compare/benchmark/smoke-test
  scripts, agent-oriented documentation.

## Done for 1.0.0 (`CHANGELOG.md`)

- Exact-mode transport 1.35× faster and the fixed cost per run 0.9 s →
  ~0.2 s with
  bit-identical results; optional fast approximations (~2.1× more,
  statistically checked against the exact mode); parallel output files.
- Crash reports (stack, minidump, the ion being transported and how to
  rerun it), a diagnostic log, exceptions caught before Qt, validated
  inputs with physical limits, plausibility warnings, automatic target
  depth and the transmitted-ions warning; input fuzzer.
- dpa by damage model (NRT, Kinchin-Pease, arc-dpa, rpa); trajectory and
  hybrid BCA-MD event files in every transport.
- Crystals (`docs/CRYSTAL_DESIGN.md`): unit cells and cell files, the
  crystal transport, the GUI's Crystal tab, comparison with measurements,
  the range study; polycrystals in the transport (on hold: no GUI).
- Release: packages for Windows and Linux, the release workflow, signing
  script, platform integration (manifest, signals, desktop files), logo
  (`docs/RELEASING.md`, `docs/PLATFORM.md`).

## P0 — Finish the 1.0 release (`docs/RELEASING.md` §6)

1. **Code signing** through the SignPath Foundation (the maintainer's
   choice; 1.0 goes out unsigned with instructions, docs/INSTALL.md):
   RELEASING.md §3.2; done when a tagged build's Windows files verify with
   `signtool verify /pa`.
2. **First run of `.github/workflows/release.yml`** on a test tag; fix what
   fails; install the MSI, the AppImage and a `.deb`/`.rpm` by hand.
3. The Windows App Certification Kit on the signed MSI (PLATFORM.md).
4. AppStream screenshots; then the Flatpak build and the winget and
   Flathub submissions.
5. Offscreen GUI rendering (KI §2 GUI) so that GUI changes can be checked
   without a window again.

## P1 — SRIM-mode fidelity (finish the TRIM reproduction)

0. **Fast approximations against TRIM.** Run the TRIM comparison
   (`compare_with_trim.bat`, about an hour; ask the user) with AMCSET's
   `--fast` and compare the z statistics with the exact mode's
   (`docs/SRIM_MODE.md` §8.2). Done when the result is in §8.2.

1. **Replacements in Xe gas** (KI §3, `docs/VALIDATION.md` §2.2). Add an
   option that rounds the recoil slots (energy, position, direction) to
   `float` as TRIM does; compare `gas_xe_xe_1000keV_fc_gas` with and without
   it (12+ seeds). Done when the cause is identified (or excluded) and the
   docs say so; if it is single precision, decide with the user whether it
   becomes the trim2013 default.
2. **Controlled benchmark.** Run `tools\windows\benchmark_vs_trim.bat` on an
   idle machine (ask the user); put `benchmark.csv` into
   `docs/SRIM_MODE.md` §8 next to the validation-run figures. Also rerun
   `tools/perf/bench.py ab` (0.1.0, exact, `--fast`) for the fast mode's
   newer tables (§8.2 still quotes the speed from before them).
3. **Use the bundled tables in the optional tests.** Make
   `StoppingTable.RealSrimTablesIfAvailable` and
   `TrimRegression.AgainstTrimWithSrimTablesIfAvailable` fall back to
   `data::bundled_srim_stopping_table()` so they always run (KI §2). Done
   when no test is skipped for lack of data.
4. **TRIM's recoil-stack limit** as an option (`physics.recoil_stack_limit`,
   default 0 = unlimited, TRIM.cfg value to reproduce TRIM), behind a
   documented switch (KI §3).
5. Extend the reference set where TRIM has not been compared yet: more
   compounds with light elements (oxides, carbides, hydrides), multilayers
   with gas layers, oblique incidence above 60°, very high energies with
   KP. One new batch generator in `tools/srim/`, then `compare_with_trim.bat`.

## P2 — Physics beyond TRIM

1. **Native engine: validate or retire.** Compare its ranges and damage with
   SRIM mode / TRIM for Fe → Fe and H/He in metals; either fix it to agree
   (KI 1.4–1.10) or make it a thin layer over the shared physics
   (`common/include/physics/`) for electrons only. Needs a decision by the
   user.
2. **Electron transport validation** against published CSDA ranges and
   backscatter coefficients (NIST ESTAR, Tabata); add a density-effect
   correction (KI 1.7).
3. **Built-in stopping** near the stopping maximum (shell corrections,
   Andersen–Ziegler-type fits) so that results without the SRIM tables are
   within a few percent (KI §3).
4. **Damage beyond Kinchin–Pease:** NRT, Kinchin-Pease, arc-dpa and rpa by
   depth and dpa at a fluence are done (docs/OUTPUT.md "dpa",
   docs/VALIDATION.md §5). Open: the implanted-ion concentration at a
   fluence; spreading each primary's damage over its cascade instead of
   booking it at the start; the low-energy arc-dpa correction of Yang and
   Olsson used by recent W MD; arc-dpa constants for more elements and for
   alloys; the refined transport.

## P3 — Code design (behaviour-preserving refactors)

1. **Split `trim_engine.cpp`** (2670 lines) into model construction,
   `Transporter` (trim2013 and refined in separate files) and the driver,
   with no change in results (exact tallies in `trim_test.cpp`, the full TRIM
   comparison bit-for-bit identical summaries).
2. **Native `Settings` builder:** replace the three positional ~30-argument
   constructor overloads with an aggregate + `validate()`.
3. **Split `amcset_common.h`** into one header per type; move the model
   functions out of `Simulation` (they live in `ion.cpp`).
4. Remove `using namespace boost::units` from `amcset_utilities.h`.
5. Replace the `#define private public` test hack with test accessors.
6. Isotope table as a generated `.cpp` with a plain array indexed by Z (also
   fixes clang).
7. Native engine memory: create Bombardments on the workers, thread-local
   RNG and tallies, merge once (memory O(threads)); keep trajectories only
   when needed; explicit-stack cascades instead of recursion (KI §2).

## P4 — Features

0. Unicode file paths on Windows before 10 1903 (KI §2; newer Windows uses
   the UTF-8 code page of the manifest): `std::filesystem::path` from the
   front ends down.
1. GUI tests with Qt Test: `Scenario` undo/redo, `FormBinder` round trips,
   the automatic-depth sizing and the escape question of `RunController`,
   each panel against `RunConfig` JSON.
2. GUI: expose the native engine's stop energies and electron tables.
3. Lateral target boundaries / 2D–3D geometry; dynamic composition changes
   with fluence (would leave TRIM compatibility: separate transport).
4. Per-layer isotope lists and non-natural isotopes in SRIM mode.
5. Batch/parameter sweeps from the CLI (a list of overrides → one output
   folder per point).

## P5 — Tooling and CI

1. A Linux test job on every push (GitHub-hosted; the release workflow
   tests Linux only for tags): core, CLI and tests without Qt; optionally
   the determinism harness.
2. Apply the optimisation flags in one place (they are set both in
   `CMAKE_CXX_FLAGS_RELEASE` and on `compiler_flags`); configurable
   `ABSL_MIN_LOG_LEVEL` in Release. (`AMCSET_MARCH` already replaces
   `-march=native` for portable builds.)
3. Git LFS for the stopping-table file if it grows past 50 MB.
4. `clang-format` / `clang-tidy` configuration (after the clang fix).
5. Doxygen API docs in CI from the existing `Doxyfile`.
