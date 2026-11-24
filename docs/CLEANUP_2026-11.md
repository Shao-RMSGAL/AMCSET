# Cleanup sweep — November 2026

> Historical record of the first cleanup (before SRIM mode, the CLI and the
> new GUI). File names, the GUI classes and test counts below describe the
> code as it was then; for the current state see `AGENTS.md`,
> `docs/ARCHITECTURE.md` and `docs/KNOWN_ISSUES.md`.

Scope: a full read-through of the repository followed by a cleanup that
corrects logic errors, documentation errors, dead code and sub-optimal code
design **without changing simulation results**. Physics issues that would
change results were documented in `docs/KNOWN_ISSUES.md` instead of being
fixed.

## How behaviour preservation was verified

1. **Bit-for-bit simulation output.** `tools/determinism_harness` fires
   particles with fixed seeds and prints every collision point at 17
   significant digits. It covers ions with all five screening models and all
   three electronic stopping models, oblique incidence, cascades, and
   electrons. Output before and after the cleanup: **57,570 lines,
   byte-identical** (GCC 13.3, Boost 1.83, abseil 20250127.1, CMake Release
   `-O3`, without the project's LTO / `-march=native` flags).
2. **Tallies JSON.** `Tallies::to_json()` and `Tallies::merge()` (rewritten)
   produce byte-identical output to the original for a populated set of
   tallies including empty-histogram adoption.
3. **Unit tests.** `common-test` passes before (22 passed, 2 skipped) and
   after (24 passed — 2 new tests — 2 skipped). The core builds with
   `-Wall -Wextra -Wshadow -Wformat=2 -Wunused` and zero warnings.
4. **GUI.** Not buildable in the cleanup environment (needs Qt ≥ 6.7 with
   QtGraphs). Every GUI translation unit was compiled with `-fsyntax-only`
   against Qt 6.4 headers plus a `QLineSeries` stub; the only errors are the
   three Qt ≥ 6.7 APIs the code uses (`QIcon::ThemeIcon`,
   `QCheckBox::checkStateChanged`), identical to the pre-cleanup tree.
   **Please build and smoke-test the GUI with the real toolchain.**

## Bug fixes (behaviour changes limited to error / edge paths)

| Where | Problem | Fix |
|---|---|---|
| `Ion::closest_approach` (ion.cpp) | The bisection fallback never updated its convergence measure or iteration counter, so entering it was an **infinite loop** (worker thread hangs forever). | Re-evaluate the bracket width each step, cap at `max_iterations`, log on failure. Results of non-hanging cases are unchanged. |
| `Volume::Volume` (volume.cpp) | Ordering check compared every layer against the *first* layer only (`prev_depth` never updated), so e.g. depths {2, 1000, 500} were accepted and produced wrong layer lookups. | Update `prev_depth` each iteration. New test. |
| `Volume::get_layer(size_t)` | `index > size()` check allowed `index == size()` → out-of-bounds read. | `>=`. New test. |
| `Layer::Layer` (layer.cpp) | All-zero fractions gave 0/0 = NaN weights. | Throw `std::invalid_argument`. New test. |
| `IsotopeData::get_isotope_mass` | Invalid-mass-number error message printed the atomic number and `MAX_Z`. | Print the mass number and `MAX_A`. |
| `IsotopeData::element_to_atomic_number` | `"Nh"` (Nihonium, 113) was missing; non-IUPAC symbols `"Jb"`, `"Uu"`, `"Up"`, `"Uh"` mapped to elements (and `"Uu"` had an unreachable duplicate branch). | Added Nh, removed the bogus two-letter aliases (the official 3-letter placeholders Uut…Uuo still work). New round-trip test over all 118 elements. |
| `Electron::initialize_tables` | The "empty path → bundled asset" fallback was unreachable (the empty path was turned into `<cwd>/` first). A redundant mutex was taken inside `std::call_once`. | Resolve the fallback before making the path absolute; drop the mutex. |
| `Bombardment::Bombardment` | A `VLOG(1)` statement drew a random number, so turning on verbose logging changed the random stream. | Removed. |
| `Particle::Particle` | VLOG printed `(r, z, y, z)` instead of `(r, x, y, z)`. | Fixed. |
| `Simulation::run_simulation`, `Bombardment::run_bombardment` | Failures writing output files were logged at `VLOG(1)` (i.e. invisible). | `LOG(ERROR)` / `LOG(WARNING)`. |
| `Quaternion` | Free `operator+(Q,Q)` / `operator/(Q,double)` duplicated the members, making `a + b` and `q / s` ambiguous (compile error on use); `normalize_copy()` was declared but never defined (link error on use). | Removed the duplicate free operators; implemented `normalize_copy()`. New test. |
| `MainWindow::~MainWindow` | `quit()` does not stop a QThread without an event loop, so closing the window during a run blocked until the whole simulation finished. | Call `request_cancel()` first (running bombardments still finish). |
| `MainWindow::handleSimulationError` | `deleteLater()` on a `SimulationThread` that may still be inside `run()` (the error signal is emitted from `run()`), which can abort with "QThread: Destroyed while thread is still running". | `wait()` before `deleteLater()`. |
| `MainWindow` layer remove button | Removing a layer left the starting depths of the rows below it stale; `updateLayerDepths()` would dereference null widgets if called with zero rows. | Recompute after removal; early return on zero rows. |
| `LayerWindow` | Layout was installed on a throw-away child and then passed to `setLayout()` (Qt rejects this), leaving the window unmanaged. | Install the layout on the window. (Window is currently unreachable.) |
| `scripts/animate_cascade.py` | Shebang was below the license header, so `./animate_cascade.py` did not work. | Moved to line 1. |
| `conanfile.py` | `exports_sources` listed a non-existent `server/*` and omitted `test/*` (required by the top-level CMakeLists, so `conan create` would fail). | Fixed. |

## Dead code removed

- `DataFrame<Types...>` (unused template, ~140 lines) and its `<shared_mutex>`,
  `<variant>`, `<unordered_map>`, `<iomanip>` includes.
- `Simulation::DrawFunction`, the `draw_function` constructor parameter,
  member and `get_draw_function()` — never invoked (the result queue replaced
  it). The GUI no longer passes an empty lambda.
- `Ion::screening_function` (an unused copy of `thomas_fermi_screening`).
- `Electron::initialize_tables()` (no-arg, loaded `"temp.csv"`),
  `Electron::mott_table_for`, `Electron::screening_table_for` (unused; the
  latter returned a reference to a function-static vector, a data race).
- The unused `incident` flag in `Ion::fire`; `std::this_thread::yield()` in
  `run_bombardment`; `(void)save_per_bombardment_info_`.
- Duplicated speed/time computations around `CollisionPoint` creation in
  `Ion::fire` and `Electron::fire` (same values, computed once now).
- Redundant empty-histogram branches in `Tallies::merge` (already handled by
  `Histogram::merge`) and the six copy-pasted loops in `Tallies::to_json`.
- `DECLSPEC` / `EXPORTING_AMCSET` (unused; the library is static).
- Duplicate `using namespace amcset::common;` lines, duplicate
  `number_density_quantity` alias, duplicate includes, duplicate forward
  declarations, an orphaned Doxygen block, commented-out members/functions.
- `mutable` on `Layer::discrete_distribution_` and the non-const return of
  `get_discrete_distribution()` (Boost's `operator()` is const).
- GUI: the Qt "Menus example" leftovers (Format menu: Bold, Italic, four
  alignment actions, line/paragraph spacing — slots, actions and members);
  the "About" box text from that example now describes AMCSET.
  `MainWindow::calculateLayerDensity` (declared, never defined) and a stray
  free function of the same name; a loop that created invisible, unused
  editor widgets for every new composition row; a never-shown `QTableView`
  allocated in `openLayerWindow`; a duplicated `setRange`/`setValue` and a
  duplicated `saveState` connection on the range box; the duplicate
  isotope-snapping lambda on the mass box (now just calls
  `updateIncidentIonSelection`); `Widget3DView::requestClearScene` and
  `m_clearSceneRequested` (set but never consumed); `Application`'s unused
  `Simulation*`; adding `formBox` to a splitter it was immediately re-parented
  out of; unused includes (`boost/interprocess`, `boost/lockfree`,
  `boost/thread/mutex`, `<iostream>`, …).
- Build: commented-out ASAN flags (superseded by `ENABLE_ASAN`), a commented
  `find_package(OpenSSL)`, an empty `target_include_directories`, a dead
  `PROJECT_VERSION` check, a `WIN32` block that set the C++ standard after all
  targets were created (no effect), `include(GoogleTest)` in `gui/`; unused
  variables in `qbuild.sh`/`qbuild.ps1`.

## Design / code-quality changes (no behaviour change)

- `SimulationThread::run()` moved from `public slots:` to `protected:`
  (calling it as a slot would have run the simulation on the caller's
  thread).
- `Layer::get_discrete_distribution()` returns a const reference.
- Function-pointer selection helpers in `simulation.cpp`: corrected error
  messages ("screening" vs "stopping").
- Settings validation messages consistently use `EXCEPTION_MESSAGE`.
- Header hygiene: every header now includes what it uses (`<stdexcept>`,
  `<sstream>`, `<string>`, `<functional>` …) instead of relying on
  transitive includes; `boost/asio.hpp` narrowed to `thread_pool`/`post`.
- `test/common/common_test.cpp`: `throw e;` (slicing) → `throw;`; duplicated
  catch bodies removed; leftover `std::cout` debugging replaced by an
  assertion; the electron test's skip message corrected; the
  `#define private public` hack documented (and `<any>` deliberately kept
  before it).

## Documentation

- Doxygen comments rewritten/added for every public type in `common/`
  (units, coordinate conventions, thread-safety, ownership, which functions
  are unused, and pointers to the known issues at the exact code sites).
- File-level comments for every source file that lacked one.
- Stale comments fixed: C++23 instructions in the build scripts (the scripts
  use gnu20), "Construct an electron" on `Ion`, "proton bombardment",
  "Immediately fires … once created" (it does not), "strictly greater" layer
  depths (non-decreasing is what is checked), angle-based `Velocity` docs,
  `getIsotopicMass`, etc.
- `requirements.txt` aligned with Docker/CI (`conan==2.7.0`);
  `config.h.in` copyright year aligned with the header check (2025);
  Dockerfile comment typo.
- New: `AGENTS.md`, `CLAUDE.md`, `CONTRIBUTING.md`, `docs/KNOWN_ISSUES.md`,
  `docs/ROADMAP.md`, `docs/TESTING.md`, this file, and
  `tools/determinism_harness/`.

## Not changed on purpose

- Every item in `docs/KNOWN_ISSUES.md` §1 (physics/numerics).
- Public `Simulation::Settings` constructor overloads (callers depend on
  them), file formats, output locations, GUI layout and settings keys.
- The File/Edit placeholder actions (kept as hooks for planned features).
- `LayerWindow` / `openLayerWindow` (unreachable but planned).
