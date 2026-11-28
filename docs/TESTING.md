# Testing

## Running the tests

All tests are registered with CTest.

```bash
# Linux (any build tree, e.g. build-noqt from AGENTS.md §3)
cmake --build build-noqt && ctest --test-dir build-noqt --output-on-failure
build-noqt/common-test --gtest_filter='Trim*'      # one suite directly

# Windows (MSVC + Conan): configure, build and test in one go
set AMCSET_NO_PAUSE=1
tools\windows\build_and_test.bat          # log: build\build_and_test.log
```

`build_and_test.bat` runs `conan install` only when the Conan generators do
not exist yet (or with `--conan`), then `cmake --preset conan-release`, the
build (`ninja -k 0`, so all errors are reported at once) and `ctest`. Its last
status line (`SUCCESS`, `BUILD FAILED`, `TESTS FAILED`) is also written to
`build\build_and_test.status`.

## Suites

289 CTest tests on Windows (288 on Linux, where the segfault crash-report test is
left out); 6 are skipped by design (see Notes).

| Executable / file | Suites | Covers |
|---|---|---|
| `common-test` `trim_test.cpp` | `TrimRandom`, `TrimKinchinPease`, `TrimConfig`, `TrimBugs`, `TrimIn`, `TrimDat`, `TrimEngine`, `TrimTables`, `TrimRegression`, `TrimAutoDepth`, `TrimMdEvents`, `TrimTrajectories`, `TrimDamageModels`, `TrimCrystal` | SRIM mode: L'Ecuyer generator (recurrence, the division-free modular reduction, uniformity, Gaussian moments, independent per-ion streams); Kinchin-Pease vs TRIM's COLLISON.txt values; configuration validation, element defaults, the bundled tables as default; the TRIM bug switches (names, JSON, each fix changes what it should); TRIM.IN parsing, round trip, rejected modes; TRIM.DAT ions and recoils; results independent of the thread count, seed dependence, continue = one long run; identical ion trajectories in both damage modes; energy conservation; thin targets; per-element Ed; cancellation; TRIM-format and GUI output files; datasets and their units; statistical regression of Fe → Fe against TRIM; `fast_math` against the exact transport and its record in TDATA/summary.json; `first_ion` reruns the same ions; `snapshot_tallies` derived from the TRIM tallies; the limits of `validate` (energies, masses, thicknesses below one atomic layer, densities, files) with their messages; malformed TRIM.IN counts and codes; TRIM.DAT lines that cannot be transported; sputtered atoms use their own layer's Es; a solid of gas density runs (no NaN); plausibility warnings (1 keV into 1 m, densities, gas flags, transmission, relativistic ions); the automatic target depth (estimate vs the deepest ions, the deepest test ion in the last margin from too thin and too thick starts, layers above, recoil sources, cancellation, transmitted ions counted during a run); trajectory points in every transport (one birth per track, parents that exist with the child's position, times that never decrease, the collision energies; atoms that never move leave no points); MD events (record mode leaves every tally and the event file identical at 1 and 4 threads; handoff conserves the energy to 1e-3; in a channeling crystal one crossing per atom and E_c, and the events of seven E_c nested; recoil-source events; the summary table and JSON; settings validated); damage models (Nordlund 2018 Table 1 constants; NRT, Kinchin-Pease, arc-dpa and rpa of one 10 keV Fe primary against the formulas by hand, its damage energy against Robinson's partition, booked in its starting bin; quick-calculation NRT = TRIM's vacancies; arc < NRT < KP, rpa > NRT; the dpa table summed over depth = the per-ion totals at the given fluence and density; arc_dpa/rpa left out for elements without constants; validation); the crystal transport (the orientation report: Bunge angles that give the rotation back, the directions along the depth and the beam, the atom spacing along <100>, <110> and <111> of diamond and Lindhard's angle, in summary.json, the crystals table and the lattice offset columns of the export files; polycrystals: large grains of one orientation channel like the single crystal, small random grains do not, a fiber texture keeps the channel, wide disordered boundaries remove it, the same result with 1 and 4 threads, energy and defects accounted for; grains of 300 Å of one orientation in registry give the single crystal's mean range within 4 standard errors, apart they give less; the per-grain results count every ion once where it entered, add up to the run's range tallies, do not depend on the thread count and show a channel grain's ions going deeper; names, limits and warnings; it takes none of TRIM's bugs whatever the switches say, with results identical to a scenario with every bug fixed; results independent of the thread count; energy conserved and every displaced atom accounted for; 30 keV B in Si 7° off with a 22° twist = the amorphous target, whatever the stopping division and the maximum impact parameter; ⟨001⟩ and ⟨110⟩ channels deeper, straighter and with fewer displacements, a 100 Å oxide in between; hotter crystals dechannel, a static lattice channels best; the dip of 1 MeV He in W ⟨001⟩ against Lindhard's ψ₁ = 2.10°: a few per cent on the axis, out of the dip at ψ₁, the shoulder at 3°; a point spot on a static lattice repeats one trajectory, spread spots do not; amorphous layers = TRIM within the statistics; alloy and interstitial atoms are struck, atoms are sputtered; with the cut-off at Ed the vacancies are NRT's within 20 %) |
| `common-test` `crystal_test.cpp` | `CrystalGeometry`, `CrystalLattice`, `CrystalStructures`, `CrystalOrientation`, `CrystalVibrations`, `CrystalFiles`, `CrystalSearch`, `CrystalMaterials`, `CrystalDatabase`, `CrystalGrains` | Crystals: the crystal database (every row builds its structure and names known sources, every element's density within 4 % of the element table, Gemmell's rows first and unchanged; search by name, elements, structure and note; composition matching with the preferred form, alloys and off-stoichiometry refused; names in scenarios); polycrystals (the grain of a point and its distance to the nearest boundary against a search over every seed and its lateral images, for Voronoi cells, columns, Laguerre cells and a map of three seeds; voxels and their faces; the rotations of cubic, hexagonal, tetragonal, orthorhombic and triclinic lattices and the misorientation with them, with Mackenzie's largest and mean angle for random cubic orientations; generated grains: their number, the same grains from the same description, random, fiber and mosaic textures, columns, limits; AMCSET's grain file with seeds and with voxels, and its error messages; a Neper tessellation with every orientation descriptor in Neper's conventions (passive unless marked, the matrix g, Kocks and Roe angles; the values of Neper's documented example) and Laguerre weights; EBSD maps in .ang and .ctf form with points that are not indexed, symmetry-equivalent orientations, the threshold and the phase; boundaries between grains: columns through the periodic repeat, Voronoi cells (points on the boundary, fewer than 2 % of touching pairs missed), voxel faces; the grain summary (number, sizes, Mackenzie's mean misorientation for random grains, zero for a mosaic without spread, the fiber's normal direction, columns); low-angle boundaries in registry (grains of one orientation become one lattice everywhere, misoriented grains are not joined, with a spread the boundaries the groups are built through are in registry); a lattice search that shares another's sites finds the same atoms number for number); the common crystals of Gemmell 1974 (the Debye model reproduces u₁(293 K) of all 19 elements of Table IV and the worked examples for Ge and Au; all 43 rows build valid cells with handbook densities and the right stoichiometry; the thermal amplitude falls back to the table only for a pure element; Lindhard's angle for 1 MeV He in Ge ⟨110⟩); the collision-partner search of the crystal transport (every atom ahead within the maximum impact parameter, against a brute-force search over all sites, in five cells with random orientations, vibrations, alloy sites and vacancies; thermal displacements with the right rms, cut at 3.5 σ, uncorrelated between neighbours; occupancy fractions; one seed = one frozen lattice; the largest hole of seven structures against textbook geometry); vectors and rotations; lattices from parameters and vectors (volumes, angles, reciprocal vectors, impossible cells); every built-in structure against textbook geometry (atoms per cell, coordination numbers and neighbour distances; every interstitial set: sites per cell and their neighbours); density from the lattice constant and back (Si, Cu); cell validation (occupancy, overlapping sites); orientation from Miller indices (cubic, tetragonal through the reciprocal lattice, hexagonal four-index notation), Bunge Euler angles (cube and Goss textures), the twist, tilt/twist for a beam direction (the beam then travels along it, random directions in cubic and triclinic cells), naming a direction; the Debye function (series, the tabulated value at x = 1, limits) and the Debye amplitude against the oscillator average over the Debye spectrum, zero-point and classical limits; cell files: symmetry operations, a CIF with 16 operations (rutile: sites and Ti-O bond lengths), disorder and site kinds, refused CIFs, POSCAR (VASP 4/5, Direct/Cartesian, volume scale), LAMMPS data (styles, triclinic box, types by comment, map or mass), extended XYZ, writers round trip in all four formats |
| `common-test` `compute_test.cpp` | `ComputeBackend`, `ComputeKernel`, `ComputeGpu` | The GPU backends (docs/GPU.md): backend names and the `compute` scenario key (round trip, short form, errors); the device list; the kernel source compiled as C++ (`reference`) gives the CPU transport's tallies bit for bit (every histogram, total, moment, spectrum, position and exit record) with full cascades, the quick calculation, a two-layer compound at 7°, monolayer sputtering at 60°, every TRIM bug fixed with a recoil limit, a depth window in a thin target from `first_ion` 1000, the scattering integral, TRIM.DAT ions with repeats, TRIM.DAT recoil events (full cascades and quick damage, TRIM's bugs on and off, one below Ed) and EXYZ.txt (its bug on and off; an interval too fine for the GPU runs on the CPU with a warning); the device's ring of block slots and its merge in block order are part of every such run; blocks that outgrow the recoil slots and blocks with kept trajectories run on the CPU with the same results; an unsupported transport runs on the CPU with a warning; the report lines both front ends show and the `compute` object of summary.json (none for a CPU run); a GPU run of fewer than 100,000 ions gets a note that the CPU is usually faster (not for `reference`); Kokkos on a CPU execution space = the CPU (skipped without Kokkos); mixed precision within the statistics of 1600 ions and not number for number; a GPU present on the machine (`AMCSET_TEST_GPU=opencl\|hip\|cuda\|kokkos`, optional `AMCSET_TEST_GPU_DEVICE`; skipped otherwise): reproducible, and range, straggling, vacancies and ionization of 2000 ions within the statistics of the CPU's |
| `common-test` `physics_test.cpp` | `Units`, `Elements`, `Stopping`, `StoppingTable`, `Scattering` | Units and CODATA constants; element data (TRIM's ATOMDATA values, I-values); relativistic kinematics; Lindhard-Scharff; Bethe-Bloch vs NIST PSTAR and the PDG minimum ionization of Fe; density effect; built-in stopping vs SRIM-2013 values; Bohr straggling; stopping tables (interpolation, velocity scaling, extrapolation, loading, tables parsed on demand with malformed/truncated files, Bragg's rule with compound correction, Bethe-Bloch continuation); screening normalisation; magic formula vs the scattering integral; finite angles for extreme ε and b (denormal potentials); Rutherford limit; head-on collisions; nuclear stopping vs the ZBL fit; conservation in the collision kinematics; FFP impact parameter |
| `common-test` `fast_math_test.cpp` | `Scattering`, `FastMath` | The optimised scattering code: phi and phi' from shared exponentials are bit-identical to separate calls; kinematics from sin/cos of theta/2 equal the sin^2 form; `MagicTable` reproduces the magic formula (relative error < 5e-5 over solids and gases) and declines outside its range; the inline sine/cosine (< 1e-11); `OctaveTable` on power laws (its error bound, no answer outside the table); the tabulated Kinchin-Pease partition, flight path and log10 against their formulas; the tabulated stopping power against the curve (< 4e-5) |
| `common-test` `diagnostics_test.cpp` | `Diagnostics`, `ParallelTasks` | log levels, the log file (disabled levels not even evaluated), the parallel task runner (every task runs, the first error in task order is rethrown) |
| `common-test` `output_test.cpp` | `OutputUnits`, `OutputFormat`, `OutputSpec`, `OutputWriters`, `OutputText`, `OutputStream` | output units and conversions, number formats, OutputSpec parsing and shaping (columns, units, filters, every/max_rows), CSV/JSON/Parquet writers (written files read back); the fast number text (`output/text.h`, `format_number_to`) is character for character `snprintf`'s over random values, ties, denormals and extremes, and the buffered writer keeps every byte in order; the streaming CSV/Parquet writer (header, pieces, several Parquet row groups, wrong columns refused) |
| `common-test` `tallies_test.cpp`, `tallies_more_test.cpp` | `TalliesUnit`, `TalliesIntegration` | Tallies/Histogram merging and binning, tallies after real simulations |
| `common-test` `common_test.cpp` | `CommonTest` | Native engine: Coordinate/Velocity/Quaternion, IsotopeData, Properties, Layer (number density), Volume, Settings validation, model selection, collision-point clock, completion counting, cancellation, output files, a full ion bombardment; physics fixes P1-P6 |
| `cli-test` `test/cli/cli_test.cpp` | `RunConfigUnits`, `RunConfig`, `Cli` | Quantities with units (parsing, wrong units, formatting round trips), complete scenarios, unknown-key rejection, JSON round trip, `--set` overrides, examples build engines, save/load, malformed numbers; command line: options, usage errors, option → configuration mapping, compositions, TRIM.IN detection, `stopping` output, `--fast`, `--log`, `first_ion`, `--auto-depth` and the `auto_depth` keys, the trajectory/MD event options and keys (JSON round trip, malformed values), `damage_efficiency` and `output.dpa_fluence` (units, round trip, missing constants), the layers' `crystal` (built-in structures with and without a lattice constant, alloys and interstitials, sites, cell files, overrides, round trip, the layer's composition and density from the cell, a crystal layer running as the amorphous one with a warning, 40 malformed inputs each naming its key), the crystal transport's keys (`physics.transport`, `physics.crystal`, `ion.spot`: round trip, overrides with list values, nothing written for scenarios that do not use them, errors naming the key), the `cell` command, energy lists, the LAMMPS stopping table |
| CTest `cli.*` (`test/CMakeLists.txt`) | — | End-to-end runs of `amcset`: version, template, elements, stopping, cell, grains (mosaic grains in registry form one group), devices, validate, run (SRIM mode; and with `--compute reference`), convert to TRIM.IN, trim, bad option; crash reports (`amcset diag-crash segv` on Windows and `thread-exception` everywhere must print a report with the stack, the worker's ion and the rerun options) |

### Optional tests against SRIM-2013 data

Two tests run only when SRIM data are named in the environment (they are
skipped otherwise; the bundled file works):

```bash
export AMCSET_SRIM_STOPPING_TABLE=$PWD/data/srim2013/srim2013_stopping.amcset
export AMCSET_SRIM_DATA=/path/to/SRIM-2013/Data                       # optional
build/Release/common-test --gtest_filter='*IfAvailable*'
```

- `StoppingTable.RealSrimTablesIfAvailable`: H in Fe at 100 keV against the
  SR Module value.
- `TrimRegression.AgainstTrimWithSrimTablesIfAvailable`: Fe → Fe at 10 and
  100 keV (KP) and 100 keV (full cascades), range and vacancies within 3-4 %
  of SRIM-2013 TRIM.

### Optional test against public grain files

`CrystalGrains.PublicSampleFiles` reads files that Neper, MTEX and orix
publish, from the folder `AMCSET_PUBLIC_GRAIN_DATA` names (skipped
otherwise; the files are not in the repository). On the workstation they are
in `build\public_grain_data\` (downloaded 2026-11-27):

| Path in the folder | From | What the test checks |
|---|---|---|
| `neper/ori/ref.tess` | Neper's test suite (`tests/T/ori`, made with `-ori "fibre(x,1,1,1)"`), github.com/neperfepx/neper | a ⟨111⟩ axis of each cell along sample x: the passive convention |
| `neper/ori_file/ref.tess` | `tests/T/ori_file` (`-ori "file(ori)"`, three Rodrigues vectors) | the cells keep the file's orientations |
| `neper/oridescriptor/ref.tess`, `neper/format_2d/ref.tess` | `tests/T/oridescriptor`, `tests/T/format_2d` | quaternions; a 2D tessellation as columns |
| `mtex/testdata_sqr.ctf`, `mtex/testdata_hex.ctf` | MTEX `data/EBSD`, made by its `testdata_generator.py` (2 x 2 grains of 15 x 15 points, ±1° noise; the MAD and BC columns hold each point's grain), github.com/mtex-toolbox/mtex | every point in its grain, square and hexagonal grids; the grains' orientations |
| `mtex/ACOM.ang` | MTEX `data/EBSD` | a map without a step in its header |
| `orix/sdss_ferrite_austenite.ang` | orix-data (CC BY 4.0, J. Hjelen, NTNU), github.com/pyxem/orix-data | a measured two-phase map: size, period, phase selection |

```bash
export AMCSET_PUBLIC_GRAIN_DATA=$PWD/build/public_grain_data
build/Release/common-test --gtest_filter='CrystalGrains.PublicSampleFiles'
```

### Validation against TRIM (manual, on the workstation)

The full comparison of the reference set is one command on Windows:
`tools\windows\compare_with_trim.bat` (batch.py + analyze.py; ~1 h;
output in `build\trim_validation\<set>\`), and
`tools/trim_compare/doc_tables.py analysis.json --units` turns its result
into the tables of `docs/VALIDATION.md` §2.2. Details:
`docs/DEVELOPMENT.md` §6. The tools individually:

`tools/trim_compare/compare.py TRIM_DIR AMCSET_DIR` compares a TRIM run with
an AMCSET run of the same TRIM.IN (`amcset trim`), quantity by quantity
(`--json` for machine-readable output, `--tolerance PCT` to use it as a
check). `tools/srim/` generates the reference cases and runs TRIM in batch
mode. Results: `docs/VALIDATION.md` §2.

For many cases at once, `tools/trim_compare/batch.py CASES_DIR OUT_DIR
--amcset EXE` runs AMCSET on every finished TRIM case with several seeds
(`--seeds`; TRIM seed replicates named `..._s<k>` are run once, with
`--replicate-seeds`), and `tools/trim_compare/analyze.py OUT_DIR CASES_DIR
[MORE_CASE_DIRS] -v` compares the AMCSET seed runs with all TRIM runs of the
same input: per quantity the difference, a pooled-variance t-statistic
(reported as the normal deviate z with the same tail probability), and
chi^2 per bin of the vacancy and ion-range depth profiles. With only AMCSET's
seed scatter as the error, |z| > 3 in a few of ~20 quantities x ~100 cases is
expected by chance; a difference is real when it persists with more seeds and
across related cases.

### GUI smoke test

`amcset-gui --screenshots DIR --quit [--timeout S] [file]` runs the scenario
(the example one without a file), saves images of every view and every plot
(`plot_<id>.png`) and quits, with
the log on stderr (exit code 0 = finished, 1 = the run failed, 2 = it could
not start, 3 = timeout). It also writes the exit code to `DIR/exit_code.txt`
once the window has closed, and leaves with code 10 + the code if closing
takes more than 60 s. On Windows, `tools\windows\gui_smoke_test.bat` does
this with the Qt runtime on PATH and writes to `build\gui_smoke\<n>\`.

### Notes

- Skipped by design: `CommonTest.ion_test` (placeholder), the two `*IfAvailable` tests unless
  `AMCSET_SRIM_STOPPING_TABLE` is set (they could use the bundled tables,
  `docs/ROADMAP.md` P1.3), `CrystalGrains.PublicSampleFiles` unless
  `AMCSET_PUBLIC_GRAIN_DATA` is set, `ComputeGpu.AgreesWithTheCpu` unless
  `AMCSET_TEST_GPU` names a GPU backend, and
  `ComputeKernel.KokkosOnTheHostMatchesTheCpuTransport` in builds without
  Kokkos (`-DAMCSET_WITH_KOKKOS=ON`; on the workstation `build\kokkos\amcset`).
- The fixture in `common_test.cpp` redirects `stderr` to `log.txt`.
- Several native tests run complete simulations and write to `./output/…`
  in the working directory (the source tree under ctest).
- `common_test.cpp` uses `#define private public` to read private **data
  members** only; never call private or protected member functions through
  it (MSVC encodes access in mangled names, so such calls do not link).
- Native simulations are seeded from `EntropySeeder` and are not
  reproducible; SRIM-mode runs are (fixed seed, results independent of the
  thread count), so `trim_test.cpp` can compare exact tallies.

## Determinism harness (`tools/determinism_harness`)

Purpose: prove that a refactor of `common/` did not change numerical
results. The harness constructs `Ion`/`Electron` objects directly with a
fixed-seed `mt19937` (bypassing `EntropySeeder`) and prints every collision
point at full precision. Coverage: Si/Fe/H/He/Cu ions with every screening
model (ThomasFermi, ZBL, Moliere, QuickBasic/Firsov, Custom) and every
electronic stopping model (ZBL, LSS, CustomPowerLaw), normal and oblique
incidence, recoil cascades, and 20 keV electrons; a three-layer
water/iron/water target.

```bash
# Before the change
cmake -S tools/determinism_harness -B build/harness -DCMAKE_BUILD_TYPE=Release
cmake --build build/harness
./build/harness/amcset-harness > /tmp/before.txt      # from the repo root

# After the change (rebuild, same command)
cmake --build build/harness
./build/harness/amcset-harness > /tmp/after.txt
cmp /tmp/before.txt /tmp/after.txt && echo "bit-identical"
```

To compare two checkouts, configure two build directories with
`-DAMCSET_SOURCE_DIR=/path/to/checkout`. The harness only needs Boost headers
and abseil (≥ 20230125 for `absl::log`); on Ubuntu 24.04 build abseil from
source (`20250127.1`, `-DCMAKE_CXX_STANDARD=20`) and point
`CMAKE_PREFIX_PATH` at it.

Caveats:

- Floating-point results depend on compiler, flags and standard library.
  Compare outputs produced with the **same** toolchain and flags. The
  project's Release flags (`-march=native -flto`) may legitimately produce
  different (but self-consistent) numbers from a plain `-O3` build.
- The harness deliberately bypasses `Simulation::run_simulation()` (thread
  pool, tallies, files). Changes to those paths need the unit tests, plus a
  manual comparison of `output/*.csv` for a single-threaded run if relevant.
- After an **intentional** physics change, regenerate the reference output
  and say so in the commit message.

Reference: at the end of the November 2026 cleanup the harness produced
57,570 lines; the output was byte-identical to the pre-cleanup code. The
physics fixes of 0.1.0 change the output intentionally (every scattering
angle changes); regenerate the reference from 0.1.0 before the next
refactor. Since 0.1.0 the harness compiles all of `common/src` recursively
(physics/ and trim/).

## Suggested additional tests

- Native engine physics: projected ranges against SRIM (the native engine is
  not yet validated), electron CSDA ranges and backscatter coefficients.
- Native energy bookkeeping (incident = electronic + nuclear + escaping).
- GUI: form ↔ scenario round trip (`buildRunConfig` / `applyRunConfig`) with
  Qt Test.
