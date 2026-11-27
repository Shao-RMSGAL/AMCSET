# common/ — amcset-common (the simulation core)

Everything except the front ends: both engines, the shared physics, the
scenario model, the output tables and writers, element/isotope data. It is
a static library linked by the CLI, the GUI and the tests; it does not use
Qt and must not.

## Where things are

| Need | Look in |
|---|---|
| A quantity type, unit, constant | `include/amcset_units.h` |
| The engine interface | `include/simulation_engine.h` |
| Scenario keys, JSON, unit parsing, building an engine | `include/run_config.h`, `src/run_config.cpp` |
| Electronic stopping (SRIM tables, built-in, Bethe-Bloch, LSS), straggling | `include/physics/stopping.h`, `src/physics/stopping.cpp` |
| Screening, magic formula, scattering integral, kinematics | `include/physics/scattering.h`, `src/physics/scattering.cpp` |
| SRIM mode | `include/trim/`, `src/trim/` (own `CLAUDE.md`) |
| Crystals: unit cells, built-in structures, orientation, Debye amplitudes, cell files, the lattice's collision-partner search, polycrystals (`grain_map.h`: grains and their generation; `grain_io.h`: grain files, Neper, EBSD; `grain_summary.h`: boundaries, the summary, low-angle boundaries in registry) (not TRIM; followed by the crystal transport only) | `include/crystal/`, `src/crystal/`; a layer's crystal: `trim::apply_crystal`, `check_layer_crystal` (`src/trim/layer_crystal.cpp`); design: `docs/CRYSTAL_DESIGN.md` |
| Result tables, units of output columns, CSV/JSON/Parquet | `include/output/`, `src/output/` |
| Locating `data/` at run time | `include/data_files.h`, `src/data_files.cpp` |
| Elements (masses, densities, Es, I-values) | `include/element_data.h`, `assets/element_data.csv` → `generated/element_database.h` |
| Isotopes | `include/isotope_data.h`, `assets/isotope_data.txt` → `generated/isotope_database.h` |
| Native engine | `include/amcset_common.h`, `src/{particle,ion,electron,layer,volume,simulation,quaternion,entropy_seeder}.cpp` |
| Live statistics shared by both engines | `include/tallies.h` |

## Rules for this library

- **Units.** Use the types in `amcset_units.h` (`energy_quantity`,
  `length_quantity`, `stopping_power_quantity`, `number_density_quantity`,
  …). Construct with unit constants, read with `in_units`. Where a formula is
  published for given units (SRIM book: E in keV, N in Å⁻³), convert at the
  boundary and comment it. `amcset_utilities.h` has `using namespace
  boost::units;` inside `amcset::common`; in code that also uses
  `amcset::common`, qualify `degree`, `radian`, `dimensionless_quantity`.
- **Physics shared by both engines** belongs in `physics/`, not in an
  engine. Each model gets a unit test against a reference value
  (`test/common/physics_test.cpp`: NIST PSTAR/ESTAR, SRIM values, hand
  calculations).
- **SRIM tables.** `SrimStoppingDatabase::load` caches tables by path, size
  and mtime for the process; `data::bundled_srim_stopping_table()` finds the
  bundled file (search order in `data/srim2013/README.md`; the source-tree
  path is compiled in as `AMCSET_SOURCE_DATA_DIR`).
- **Scenario model.** A new setting: field in `trim::TrimConfig` (both
  engines) or `config::NativeOptions` (native only), parse + `reject_unknown`
  list + `to_object` in `run_config.cpp`, validation in
  `TrimConfig::validate`, docs in `docs/CLI.md`. `apply_override` and the
  JSON round trip then work by themselves; add a `RunConfig` test.
- **Output.** Tables are built once (`trim::make_dataset`) and shaped by
  `output::apply(table, OutputSpec)`; writers never compute physics.
- **Generated files.** Never edit `generated/*.h`. Regenerate with
  `scripts/generate_isotope_data.py` (run by CMake) or
  `scripts/generate_element_data.py`.
- **Crystals.** The crystal library (`crystal::`) depends on no engine.
  Its geometry is in Angstrom doubles (`_A` in names; fractional
  coordinates and directions are dimensionless); public setters and the
  scenario take quantities. A layer with a crystal has the unit cell's
  composition and density: change the cell, then call
  `trim::apply_crystal(layer)` (validation refuses a layer that disagrees
  with its cell). A built-in structure needs a test of its coordination
  numbers and distances against the textbook values
  (`CrystalStructures`); do not add material constants (lattice constants,
  Debye temperatures) without a cited source. The ones that are built in
  are Tables IV and V of Gemmell, Rev. Mod. Phys. 46, 129 (1974)
  (`crystal/materials.h`) and the crystal database
  (`crystal/crystal_database.h`), every row of which names its source
  (`database_sources()`, docs/CRYSTAL_DATABASE.md); add rows only with a
  source, and keep Gemmell's rows first and unchanged.
- **No front-end code** here: no Qt, no printing to stdout/stderr from
  library code (return warnings; throw exceptions with the key name).

## Native engine notes

Not validated; see `docs/KNOWN_ISSUES.md` §1. All its types are in
`amcset_common.h`; `Simulation::Settings` has three constructor overloads
(add a field to all three and to `print_settings()`); the static model
functions of `Simulation` are implemented in `ion.cpp`; native runs are not
reproducible (seeded by `EntropySeeder`). Refactors must keep the
determinism harness output bit-identical (`docs/TESTING.md`).
