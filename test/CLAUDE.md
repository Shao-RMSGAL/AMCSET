# test/ — GoogleTest suites and CLI end-to-end tests

Catalogue, optional SRIM-data tests, the determinism harness and the
validation workflow: `docs/TESTING.md`.

## Layout

| File | Executable | Suites |
|---|---|---|
| `common/trim_test.cpp` | `common-test` | `TrimRandom`, `TrimKinchinPease`, `TrimConfig`, `TrimBugs`, `TrimIn`, `TrimDat`, `TrimEngine`, `TrimTables`, `TrimRegression`, `TrimAutoDepth`, `TrimMdEvents`, `TrimTrajectories`, `TrimDamageModels`, `TrimCrystal`, `DataFiles` |
| `common/crystal_test.cpp` | `common-test` | `CrystalGeometry`, `CrystalLattice`, `CrystalStructures`, `CrystalOrientation`, `CrystalVibrations`, `CrystalFiles`, `CrystalSearch`, `CrystalMaterials`, `CrystalGrains` |
| `common/physics_test.cpp` | `common-test` | `Units`, `Elements`, `Stopping`, `StoppingTable`, `Scattering` |
| `common/output_test.cpp` | `common-test` | `OutputUnits`, `OutputFormat`, `OutputSpec`, `OutputWriters` |
| `common/fast_math_test.cpp` | `common-test` | `Scattering` (shared exponentials, half-angle kinematics, `MagicTable`), `FastMath` |
| `common/diagnostics_test.cpp` | `common-test` | `Diagnostics`, `ParallelTasks` |
| `common/tallies_test.cpp`, `tallies_more_test.cpp` | `common-test` | `TalliesUnit`, `TalliesIntegration` |
| `common/common_test.cpp` | `common-test` | `CommonTest` (native engine) |
| `cli/cli_test.cpp` | `cli-test` | `RunConfigUnits`, `RunConfig`, `Cli` |
| `CMakeLists.txt` | `amcset` | `cli.*` end-to-end runs (version, template, elements, stopping, cell, validate, run, convert, trim, bad option, crash reports via `diag-crash`) |

248 CTest tests in total (247 on Linux); 3 are skipped by design (two need SRIM data
through environment variables, one native test is a placeholder). CTest
runs everything with the repository root as working directory.

## Rules

- **Exact tallies are allowed in SRIM mode** (seeded, thread-count
  independent). If a change alters them, decide whether that is intended;
  if so, update the numbers in the same change and explain the cause. If not,
  the change is a regression.
- **Statistical tests** (`TrimRegression`, `TrimBugs.EachFixChangesWhatItShould`)
  use fixed seeds and tolerances wide enough for their sample size; do not
  tighten them until they pass by luck, and do not loosen them to hide a
  change — measure with more ions first.
- **Physics references:** test a model against an independent value (NIST
  PSTAR/ESTAR, SR Module output, a hand calculation in the comment), not
  against its own output.
- **Fast:** a single test should take well under a second unless it is a
  statistical regression; use few ions and small targets.
- **Files:** tests that write files use a temporary directory under the
  build or system temp directory, never the source tree (the native
  `CommonTest` suite still writes to `./output`; do not copy that).
- `common_test.cpp` uses `#define private public` (data members only, never
  private functions — MSVC). Do not spread it to new files; test through
  public interfaces.
- Every bug fix gets a test that fails without the fix.
