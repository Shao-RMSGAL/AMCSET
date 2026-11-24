# Development guide

Environments, builds, running, validation and benchmarking. Rules and
conventions: `AGENTS.md`. Tests: `docs/TESTING.md`.

## 1. Build matrix

| Environment | Builds | How | Notes |
|---|---|---|---|
| Windows + MSVC + Conan (the workstation, CI) | core, CLI, GUI, tests | `tools\windows\build_and_test.bat` or `.\qbuild.ps1 -c -b -t` | Qt built by Conan once (hours); kept in `build\Release\generators` |
| Linux + GCC, no Qt | core, CLI, tests | plain CMake (§2) | the fastest loop; what most agent sessions need |
| Linux + Conan | everything | `./qbuild.sh -c -b -t` | Qt built from source on first install |
| Docker (Ubuntu 24.04, GCC 14) | everything | `./run-dev.sh`, then `./qbuild.sh -b -t` inside | forwards Wayland/X11 for the GUI |

Only Release is supported by the scripts (Conan installs Release
dependencies; `CMakePresets.json`'s `release` preset inherits
`conan-release`). Version: `project(amcset VERSION 0.1.0)` in
`CMakeLists.txt` → `common/include/config.h` (`AMCSET_VERSION`).

CMake options: `AMCSET_BUILD_GUI`, `AMCSET_BUILD_CLI`, `AMCSET_BUILD_TESTS`
(all ON), `ENABLE_ASAN` (OFF), `AMCSET_COMPILE_SHADERS` (OFF; GUI shaders
from `gui/shaders/baked`), `BUILD_SHARED_LIBS` (ON; `amcset-common` is
static regardless), `AMCSET_PGO` (empty; MSVC only: `generate` or `use`
profile-guided optimisation of `amcset.exe` with the profile in
`AMCSET_PGO_DIR`; it gained 2 %, `docs/SRIM_MODE.md` §8.1).

Release flags: GCC `-O3 -march=native -flto -funroll-loops`; MSVC `/O2 /GL`
with `/LTCG`, plus `/Zi` and `/DEBUG /OPT:REF /OPT:ICF` so that every
Release binary has a PDB for crash reports and the profiler (no effect on
speed). `-march=native` binaries run only on CPUs like the build machine's.

## 2. Linux without Qt

```bash
# Dependencies (Ubuntu 24.04): compiler, CMake, Ninja, Boost headers, GoogleTest
sudo apt install build-essential cmake ninja-build libboost-dev libgtest-dev python3
# abseil: the distro version lacks absl::log; build one (once)
git clone --depth 1 --branch 20250127.1 https://github.com/abseil/abseil-cpp.git
cmake -S abseil-cpp -B abseil-build -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_CXX_STANDARD=20 -DABSL_PROPAGATE_CXX_STD=ON -DCMAKE_INSTALL_PREFIX=$HOME/absl
cmake --build abseil-build -j && cmake --install abseil-build

# AMCSET
cmake -S . -B build-noqt -G Ninja -DCMAKE_BUILD_TYPE=Release -DAMCSET_BUILD_GUI=OFF \
      -DCMAKE_PREFIX_PATH=$HOME/absl
cmake --build build-noqt -j
ctest --test-dir build-noqt --output-on-failure
build-noqt/amcset template -o /tmp/fe.json && build-noqt/amcset run /tmp/fe.json -n 200 -o /tmp/out
```

Release uses `-O3 -march=native -flto`; a change in `trim_engine.cpp`
recompiles in 1–3 minutes on two cores. For debugging configure a second
tree with `-DCMAKE_BUILD_TYPE=Debug` (and `-DENABLE_ASAN=ON`).

## 3. Windows

Prerequisites: Visual Studio with the C++ workload (the workstation has
Visual Studio 18, MSVC 14.5), CMake, Ninja, Python 3, Conan 2 (`pip
install conan`, then `conan profile detect`).

```bat
set AMCSET_NO_PAUSE=1
tools\windows\build_and_test.bat            :: conan install only if build\Release\generators is missing
tools\windows\build_and_test.bat --conan    :: force conan install
```

It runs `conanbuild.bat`, `cmake --preset conan-release`, `cmake --build
--preset conan-release -- -k 0` (all errors at once), checks that
`amcset.exe` starts, then `ctest --preset conan-release`. Status: last line
of `build\build_and_test.status`; log: `build\build_and_test.log`.

Running programs: `call build\Release\generators\conanrun.bat` (PATH for
the Qt and abseil DLLs), then `build\Release\amcset.exe ...` or
`amcset-gui.exe`. The post-build `windeployqt` step copies the Qt DLLs next
to the GUI when allowed.

Known Windows hurdles:

| Symptom | Cause | Fix |
|---|---|---|
| `LNK1104: cannot open file amcset-gui.exe` | the GUI is running | close it, rebuild |
| new `amcset.exe` does not start, log names Smart App Control / Device Guard | unsigned new binary blocked | allow it in Windows Security (protection history) |
| `windeployqt` or `qsb.exe` blocked | application-control policy | GUI: use `conanrun.bat`; shaders: keep the baked `.qsb` |
| paths too long during `conan install` of Qt | Windows MAX_PATH | Qt Quick Controls styles are disabled in `conanfile.py`; keep the Conan home short |
| a file briefly locked or replaced | Nextcloud sync of the working copy | retry; keep scratch output in `build\` |

### 3.1 Single-file installer

`powershell -ExecutionPolicy Bypass -File tools\windows\make_installer.ps1`
(after `build_and_test.bat`) writes `build\installer\AMCSET-Setup.exe`
(about 33 MB, built with Windows' own IExpress, no extra tools). It packs
`amcset.exe`, `amcset-gui.exe`, the three Qt libraries and plugins the GUI
uses, the app-local Visual C++ runtime, `data\` and a few docs. It installs
per user to `%LOCALAPPDATA%\Programs\AMCSET` (no administrator rights), adds
Start Menu and Desktop shortcuts, an "Installed apps" entry and an
uninstaller. The files are unsigned: SmartScreen warns, and Smart App
Control may block them. The bundled SRIM tables keep their non-commercial
terms (`data/srim2013/README.md`).

## 4. Running

```bash
amcset template > fe.json                 # example scenario (Fe -> Fe, SRIM mode)
amcset run fe.json -n 5000 -o out         # TRIM files + csv/ + results.json + summary.json
amcset trim TRIM.IN -o out                # an existing TRIM input (TRIM.DAT next to it for codes 4-7)
amcset run fe.json --set "ion.energy=2 MeV" --set damage=full --fix-trim-bug struck_atom_direction
amcset stopping --ion He --target Fe --from "10 keV" --to "10 MeV"
amcset-gui [scenario.json | TRIM.IN | results.json]
```

Output directory layout: `docs/OUTPUT.md`. The bundled stopping tables are
found automatically (search order in `data/srim2013/README.md`); set
`AMCSET_DATA_DIR` for installed builds outside the source tree.

## 5. Debugging tips

- Reproduce a single ion: SRIM-mode results depend only on (seed, ion
  index, block of 16), so `-n` small with the same seed replays the first
  ions exactly, and `--set first_ion=B -n K` replays ions B .. B+K-1 of a
  longer run when B is a multiple of 16 (from the run's first ion);
  `output.collision_details` / `collision_recoils` write COLLISON.txt for
  the first N ions; `output.trajectories` keeps paths (GUI).
- Compare one case with TRIM quickly: §6 "single case".
- Energy bookkeeping: `summary.json` → `energy_percent` (ionization,
  phonons, energy to recoils); with the TRIM bugs reproduced, low-energy
  cascades exceed 100 % (`docs/SRIM_MODE.md` §7).
- **Diagnostic log** (`common/include/diagnostics.h`): `amcset ... --log
  run.log` (add `-v` for debug lines, also on stderr); the GUI always
  writes `%LOCALAPPDATA%\AMCSET\AMCSET\logs\amcset-gui.log` (Linux:
  `~/.local/share/AMCSET/AMCSET/logs/`; Help → Open log folder; Help →
  Detailed diagnostic log for debug lines). It records run-level events —
  command line, setup time, threads, blocks per worker, files written,
  warnings, errors — never single collisions. Use `AMCSET_LOG(info) << ...`
  in new code; abseil `LOG(ERROR)` messages go into it as well
  (abseil INFO/WARNING stay compiled out in Release because the native
  engine has `VLOG`s in its inner loops).
- **Crashes**: every front end installs the crash handlers (`docs/CLI.md`
  "Crash reports"). The report names the crashing thread's stack (file and
  line, from the PDB that Release builds now write next to each `.exe`),
  the SRIM-mode ion each worker was transporting and the options that
  rerun it, and the last log lines; Windows also writes a minidump (open it
  in Visual Studio with the PDB next to the `.exe`). `amcset diag-crash
  segv` checks the whole chain.
- **Profiling** (Windows, no administrator rights): `tools/perf/winsample.cpp`
  is a small sampling profiler (stack walks through DbgHelp): build it once
  from a Native Tools prompt (`cl /O2 /EHsc /std:c++20
  tools\perf\winsample.cpp dbghelp.lib winmm.lib psapi.lib`), then
  `winsample --top 40 -- build\Release\amcset.exe trim case\TRIM.IN -n 1000
  -j 1 -q -o out`. It prints the hottest innermost (inlined) functions,
  source lines, C-runtime leaves with their callers, and inclusive times.
  On Linux use `perf record -g`. Reading its output: samples inside
  `ucrtbase.dll` are named after the nearest exported symbol, so the C
  runtime's internal `pow` appears as `finite` (the "library leaf <-
  caller" list shows the calling line); profile with `-j 1` for the
  transport and with all threads to see the single-threaded phases (start-up,
  merging, writing the output files).
- **Speed and bit-identity of a change**: `tools/perf/bench.py` (13
  reference cases, §7.1).
- **Input robustness**: `tools/perf/fuzz_inputs.py` (§7.2).

## 6. Validation against TRIM

Reference data (on the workstation, never in the repository):
`C:\Users\natha\Applications\SRIM-2013-Pro\_amcset\results\trim\` (batches
1–2) and `results\trim3\` (batch 3: replicates, alloys, layers, gases,
TRIM.DAT, options): one folder per TRIM run with its TRIM.IN, TRIM.DAT and
all outputs, `_done.txt` = TRIM's wall time.

**Full comparison** (~1 h, uses all cores):

```bat
tools\windows\compare_with_trim.bat     :: -> build\trim_validation\trim3\analysis.txt / .json
python tools\trim_compare\doc_tables.py build\trim_validation\trim3\analysis.json --units
```

Move the previous `build\trim_validation\trim3` aside first (finished runs
are skipped otherwise). Interpretation: `docs/SRIM_MODE.md` §6.2 and
`docs/TESTING.md`. Publish new numbers as in `docs/HOWTO.md` §10.

**Single case** (seconds):

```bat
set C=%USERPROFILE%\Applications\SRIM-2013-Pro\_amcset\results\trim3\328_rep_fe_fe_1000keV_fc_s1
build\Release\amcset.exe trim %C%\TRIM.IN -s 101 -q --formats trim -o build\one\101
python tools\trim_compare\compare.py %C% build\one\101
```

One AMCSET run against one TRIM run shows only the size of a difference;
whether it is real needs several seeds (`batch.py --only REGEX`) and
`analyze.py`.

**Effect of a TRIM bug fix:** run the same cases with and without
`--fix-trim-bug NAME`, 3+ seeds each, and compare the means (that is how
the numbers in `docs/SRIM_MODE.md` §7 were made).

## 7. Benchmark

`tools\windows\benchmark_vs_trim.bat [SRIM folder]` runs 12 reference
inputs with a fixed ion count: TRIM.exe alone (batch mode, private copy of
the SRIM folder in `build\benchmark\trim_work`), then AMCSET on one thread
and on all threads. Results: `build\benchmark\benchmark.csv` (wall-clock
seconds, speed-ups). It needs an otherwise idle machine and opens TRIM's
window; it takes 30–60 minutes. The numbers in `docs/SRIM_MODE.md` §8 come
from the validation runs instead (TRIM ran ~14 instances at once).

### 7.1 Speed and bit-identity of a change (`tools/perf/bench.py`)

13 TRIM reference inputs (Fe → Fe 1 keV - 1 MeV, 316 SS, layers, Ar gas,
TRIM.DAT ions and recoils; ion counts chosen for about 3 s each on one
thread with 0.1.0) run from `results\trim3`:

```bat
call build\Release\generators\conanrun.bat
python tools\perf\bench.py run --exe build\perf\baseline\amcset.exe --tag before
python tools\perf\bench.py run --exe build\Release\amcset.exe --tag after
python tools\perf\bench.py compare build\perf\runs\before build\perf\runs\after
python tools\perf\bench.py ab --exe old\amcset.exe --tag old --exe build\Release\amcset.exe --tag new
python tools\perf\bench.py table build\perf\runs\old build\perf\runs\new
```

`compare` must report every file identical for a change of the exact mode
(`common/src/trim/CLAUDE.md` rule 7; lines with times and paths are
ignored). `ab` interleaves the executables case by case and keeps the
fastest of 3 runs, which removes most of the drift of a busy machine;
`--args=--fast` adds options. Timings are `summary.json`'s transport time
on one thread (`--threads N` for more). Keep a copy of the executable
before a change (`build\perf\baseline\` holds 0.1.0's).

`tools/perf/fast_vs_exact.py --exe ... [--seeds 6 --scale 0.5]` checks the
fast approximations statistically (z per `summary.json` quantity between
the modes, and the same-seed differences); `docs/SRIM_MODE.md` §8.2 has its
last result. `tools/perf/winsample.cpp` is the sampling profiler (§5).

Compilers measured on the workstation (`docs/SRIM_MODE.md` §8.1): MSVC
14.5 (the default; `/arch:AVX2` and PGO gave nothing worth keeping) and
GCC 15 in WSL Ubuntu (6 % faster exact, 14 % faster fast mode). Clang is not
supported (the generated isotope table is empty under clang).

### 7.2 Input robustness (`tools/perf/fuzz_inputs.py`)

```bat
call build\Release\generators\conanrun.bat
python tools\perf\fuzz_inputs.py --exe build\Release\amcset.exe --out build\fuzz --count 800 --seed 1 ^
       --trim-in <case>\TRIM.IN <case>\TRIM.IN ...
```

Feeds `amcset validate` truncated and mutated TRIM.IN files and scenario
JSON files with hostile values (NaN, infinities, 1e308, negative, huge
counts, wrong units, wrong types, flipped bytes), runs the ones that
validate with 3 ions (a quarter of them with the automatic target depth)
and reports crashes (exit code 3), hangs (time-outs), and inputs that
validate but fail to run (a gap in the validation unless a named file is
missing). 2400 inputs (seeds 1-3, seeded from seven TRIM reference
inputs) gave no crash and no hang; the two problems found, a NaN
scattering angle in solids of gas density and missing-file errors
reported only at run time, are fixed and tested.

## 8. CI

- `.github/workflows/check-headers.yml` (GitHub-hosted Ubuntu): runs
  `scripts/check_headers.py`.
- `.github/workflows/tests.yml` (self-hosted Windows runner): installs Conan
  2.7.0, `qbuild.ps1 -Conan`, `-Build`, then `ctest --test-dir build/test`.
  It points at Visual Studio 2022's dev shell; the workstation has Visual
  Studio 18 — update the path if the runner is this machine.

There is no Linux test job yet (`docs/ROADMAP.md`).

## 9. Git hygiene

- Branches: `develop` (work), `main` (releases), remote
  `git@github.com:Shao-RMSGAL/AMCSET.git`. Commit or push only when asked.
- Never commit: `build*/`, `output/`, `common/include/config.h`, Conan
  output, TRIM reference data or SR Module output, screenshots, local
  validation results. `compile_commands.json` and `.vscode/.conan_tools/`
  are tracked by mistake (`git rm --cached` them in a cleanup commit).
- Do commit: `data/srim2013/srim2013_stopping.amcset` (42 MB, allowed by a
  `.gitignore` exception), `common/generated/*.h`.
- Line endings: LF, except CRLF for `.bat`/`.ps1` (`.gitattributes`).
