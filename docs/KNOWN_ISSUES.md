# Known issues

Status as of AMCSET 0.1.0 (2026-11-23). The first cleanup sweep
(`docs/CLEANUP_2026-11.md`) was behaviour-preserving and only documented the
problems below; the physics fixes P1-P6 were made afterwards, in 0.1.0
(`CHANGELOG.md`). Severity:

- **High** — almost certainly wrong and materially changes results.
- **Medium** — likely wrong or physically questionable; needs a decision.
- **Low** — minor inaccuracy, limitation, or cosmetic.

Sections 1 and 2 concern the **native** engine (`common/src/ion.cpp`,
`electron.cpp`, `simulation.cpp`) and the GUI. The **SRIM mode** engine
(`common/src/trim/`) is a separate implementation; its deviations from TRIM
are listed in `docs/SRIM_MODE.md` §6.1 and summarised in §3 below.

Verify every fix against an independent reference (SRIM/TRIM, NIST
PSTAR/ESTAR, published data), add a test, and compare the determinism
harness output before/after (`docs/TESTING.md`).

---

## 0. Fixed in 0.1.0

| Id | Issue (details kept below for reference) | Fix | Test |
|---|---|---|---|
| P1 | 1.1 rotation axis not perpendicular to the direction | axis (−y, x, 0) | `CommonTest.rotation_turns_direction_by_scattering_angle` |
| P2 | 1.2 "CM" energy used the projectile mass | `cm_energy(reduced_mass(m1, m2), v)` | `CommonTest.center_of_mass_energy_uses_reduced_mass` |
| P3 | 1.3 magic-formula half angle used as the full CM angle | `laboratory_scattering_angle(2Θ/2, ...)` | `CommonTest.laboratory_angle_from_cm_angle` |
| P4 | 1.4 screening coefficient 0.2817 → 0.02817 (φ(0) = 1) | coefficient and derivative | `CommonTest.screening_functions_are_normalized`, `Scattering.UniversalScreeningFunction` |
| P5 | magic formula G: the SRIM book prints numerator and denominator exchanged | G = 1/(γ(√(1+A²) − A)) in both engines | `Scattering.MagicFormulaAgreesWithTheScatteringIntegral` |
| P6 | 1.7 electron model units (931 vs 931494 keV, table index Z vs Z−1, cm² vs m⁻³ mean free path, missing Z/A in the ionization prefactor, non-relativistic speed) | electron.cpp rewritten with Boost.Units: exact relativistic elastic transfer, Bethe collision stopping with n_e = N Z and the Joy-Luo J, Bethe-Heitler bremsstrahlung, relativistic speed | `CommonTest.electron_collision_stopping_matches_estar` |
| — | per-element displacement energies ignored | Ed per element and layer (`Layer::set_displacement_energies`; GUI column wired) | `TrimEngine.PerElementDisplacementEnergyMatters` |
| — | `tallies.json` always written to `./output` | written next to the histogram files | — |
| — | number density used N_A with ρ in kg/m³ | N = ρ / mean atom mass | `CommonTest.layer_test` |

Also new: per-element Ed/Eb/Es in the GUI and scenario files, user seeds and
thread-count-independent results (SRIM mode), and SRIM-table, Bethe-Bloch
and built-in stopping options for the native engine.

## 1. Physics and numerics of the native engine

### 1.1 [Fixed, P1] Scattering rotation axis is not perpendicular to the direction

`Particle::rotate_inplace` (`common/src/particle.cpp`). For a general
direction (x, y, z) with x ≠ 0 and y ≠ 0 the polar-rotation axis is
`(y, x, 0)`, whose dot product with the direction is `2xy`, so the realised
deflection is smaller than the requested scattering angle. The three special
cases (direction along an axis / in a coordinate plane) are correct, which is
why the first collision of a normally incident particle looks right.

Evidence (direction (0.6, 0.48, 0.64), reproduced in Python with the same
quaternion algebra):

| requested | realised (current) | realised with axis (−y, x, 0) |
|---|---|---|
| 5.73° | 3.79° | 5.73° |
| 28.6° | 18.8° | 28.6° |
| 57.3° | 37.0° | 57.3° |
| 114.6° | 67.7° | 114.6° |

Fix: use `Quaternion(0.0, -y, x, 0.0)` in the general branch. Affects both
ions and electrons (everything that scatters).

### 1.2 [Fixed, P2] "Center-of-mass" energy is the laboratory energy

`Ion::fire` (`common/src/ion.cpp`) calls `cm_energy(m_1, speed)`, i.e.
½·m₁·v² = E_lab. `cm_energy` expects the reduced mass, and
`Particle::reduced_mass` exists but is never called. E_cm should be
E·m₂/(m₁+m₂): for Si→Fe it is overestimated by 1.5×, for Fe→H by 57×. This
feeds the reduced energy, the distance of closest approach, the radius of
curvature and therefore every scattering angle.

Fix: `cm_energy(reduced_mass(m_1, m_2), speed)`.

### 1.3 [Fixed, P3] Magic-formula half angle used as the full CM angle

`Ion::magic_formula_scattering_angle` returns
`acos((B + Rc + Δ)/(R0 + Rc))`. In the ZBL magic formula this is **Θ/2**.
`binary_collision_nuclear_stopping` uses `Tmax·sin²(angle)` which is
consistent with Θ/2, but `laboratory_scattering_angle` treats the value as
the full Θ, so laboratory deflections (and the recoil angle derived from them)
are too small.

Fix: pass `2.0 * cm_scatter_angle` to `laboratory_scattering_angle`, or make
`magic_formula_scattering_angle` return Θ and change the energy transfer to
`Tmax·sin²(Θ/2)` (clearer). Fix together with 1.2.

### 1.4 [Fixed, P4] Default "Thomas-Fermi" screening function has a coefficient typo

`Simulation::thomas_fermi_screening` and its derivative (`ion.cpp`) use the
ZBL universal-function exponents with a last coefficient **0.2817** instead
of **0.02817**, so φ(0) = 1.2536 instead of 1. This is the **default** model
(`ScreeningModel::ThomasFermi`, also used by `QuickBasic` and `Custom`).

Related observations:
- `thomas_fermi_screening_length` is 0.25·(9π²/2)^{1/3}·a₀/(Z₁^0.23+Z₂^0.23)
  = 0.8853·a₀/(Z₁^0.23+Z₂^0.23), i.e. the **ZBL universal screening length**,
  numerically identical to `zbl_screening_length`. So the "Thomas-Fermi" and
  "ZBL" models currently differ only by the typo and 0.4028 vs 0.4029.
- `Ion::screening_length_for_model(…, false)` (used only by a test) returns
  the same expression.

Fixed: the coefficient. Still open: decide what "Thomas-Fermi" is meant to be (e.g. Molière's approximation
to the TF function with the Firsov length, which is the classical TF
choice), implement it, and correct the coefficient. Rename the enum value if
it is really "ZBL universal".

### 1.5 [Medium] Recoil kinematics

`Ion::fire` knock-on branch:
- The recoil's azimuth is sampled independently by its own
  `rotate_inplace` call; momentum conservation requires it to be opposite
  (φ + π) to the projectile's azimuth.
- `asin(sin(θ_lab)·sqrt(E'·m₁/(T·m₂)))` can receive an argument slightly
  above 1 from round-off, giving NaN. `rotate_inplace` then logs "NaNs
  detected" and resets the direction to +z.
- The recoil starts with the full transferred energy T; no lattice binding
  energy is subtracted, and a displacement is counted whenever
  T > `ion_displacement_energy_` (no Kinchin–Pease / NRT treatment).

### 1.6 [Medium] Free-flight path and impact parameter model

`Ion::free_flying_path_length` always returns the mean interatomic spacing
N^(−1/3); its mass/energy/screening arguments are unused, so there is no
high-energy (long free path) regime. `Ion::impact_parameter` switches
formulas at reduced energy 10 but still uses that fixed L.
`Settings::flying_distance_number_` is validated and shown in the GUI but not
used anywhere, and `Ion::collision_diameter` is unused.

### 1.7 [Fixed, P6] Electron model: unit errors

All in `common/src/electron.cpp` (energies in keV internally):

1. **Target rest energy** in `elastic_energy_loss` is `A·931.0` — the MeV
   value — while everything else is keV (511 = electron rest energy). It
   should be `A·931494`. Effect for 20 keV on A = 56: elastic loss at 1 rad is
   0.187 keV instead of 1.8e-4 keV (≈1000× too large).
2. **Table index off by one**: both parameter tables are stored from Z = 1
   at index 0 (the first screening row −184.39, 185.39, … is Salvat's
   hydrogen set), but `mott_differential_cross_section` indexes with Z
   (clamped to ≥ 1), so element Z uses the parameters of Z + 1.
3. **Cross-section units**: the Mott cross section is built from
   r_e = 2.817938e-13 **cm** (→ cm²), but `generate_flying_distance` combines it
   with an atom density in **m⁻³**; the mean free path is therefore 10⁴× too
   short and each step is dominated by the added interatomic spacing.
4. **Ionization loss**: `153.55e24·ρ/N_A/β²·[…]` has no Z/A factor; the
   standard form is 153.55·(Z/A)·ρ[g/cm³]/β² keV/cm. With ρ in kg/m³ and the
   result in keV/m that is 15.355·(Z/A)·ρ, so the current prefactor is
   ≈16.6·A/Z (≈36× for Fe) too large — *to be confirmed against the original
   source of this formula*. The bremsstrahlung term should be checked the
   same way.
5. Speeds and times use the non-relativistic v = sqrt(2E/m) (≈3 % error at
   20 keV, unphysical above ~250 keV).

Evidence: 20 keV electrons into Fe (determinism harness) lose ≈1 keV in the
first ≈3 Å and reach a maximum depth of 120–350 Å; the CSDA range of 20 keV
electrons in iron is of order 1 µm.

After the fix the collision stopping power of electrons in Al is within 2 %
of ESTAR at 0.1-1 MeV (1.496 vs 1.465 MeV cm²/g at 1 MeV); there is no
density-effect correction, so it is ~11 % high at 10 MeV. Electron transport as a whole (ranges,
backscatter) is still unvalidated, and `CommonTest.electron_bombardment_runs`
is still skipped (it needs `common/assets` in the working directory).

### 1.8 [Medium] Electronic stopping constants

`Simulation::zbl_electronic_stopping` multiplies k_L·√E by an undocumented
`c_2 = 1.59 Å³`. Its origin should be documented or the formula rewritten as
S_e[eV·Å²]·N[Å⁻³]·Δx[Å] with literature constants. The `LSS` model has not
been checked against a reference either.

### 1.9 [Medium] Tally semantics

- `Bombardment::record_energy_loss` feeds **every** call (electronic, nuclear
  *and* the per-step total) into `energy_loss_stats`, so "Mean loss (eV)" in
  the GUI / `energy_loss_mean` in tallies.json is the mean over a mix of
  three kinds of records, not a physically meaningful quantity.
- Range = depth of the incident particle's last *recorded* coordinate. For
  ions that is the position before the final free flight; sputtered
  (backscattered) and transmitted histories are included.
- `depth_histogram` mixes incident final depths and displacement depths; it
  is only exported to tallies.json.
- "Sputtered" counts any particle (including the backscattered incident
  ion/electron) that crosses z < 0. There is no surface binding energy.

### 1.10 [Low] Closest-approach solver

`Ion::closest_approach`: the Newton step uses the derivative at the previous
iterate (`dF(old_result)`) with the function value at the current one, which
is neither Newton nor secant; the relative tolerance is 1e-3. (The
non-terminating bisection branch was fixed in the cleanup — see
`docs/CLEANUP_2026-11.md`.)

---

## 2. Software issues and limitations

### Core library (native engine)

- **Memory growth** — all `Bombardment` objects are created up front and keep
  their RNG (≈2.5 kB), six histograms, the incident particle and its full
  coordinate list (stored twice: `Particle::coordinates_` and the
  `ParticleDataLeaf`) until the `Simulation` is destroyed. Millions of
  bombardments will exhaust memory. See `docs/ROADMAP.md`.
- **Electron tables are process-global** and loaded once; a later
  simulation with different table paths silently reuses the first tables.
  They are read relative to the working directory
  (`common/assets/*.csv`).
- **Seeding** — `EntropySeeder` falls back to clock XOR thread-id hash when
  `std::random_device::entropy()` is 0 (true on some standard libraries even
  when the device is usable). All bombardments are seeded on one thread in
  quick succession, so the fallback can produce correlated or identical
  seeds. Native runs are never reproducible (the scenario seed is ignored);
  SRIM-mode runs are.
- **Clang** — `IsotopeData::database` is an empty array under `__clang__`, so
  clang builds compile but every isotope lookup throws.
- `Histogram::merge` ignores histograms of a different bin count and never
  compares ranges.
- `Volume::get_layer(depth)` is a linear search per step (fine for few
  layers).
- `Simulation::run_simulation()` can only be called once (the pool is
  joined).

### GUI

- No automated GUI tests beyond the smoke test
  (`tools\windows\gui_smoke_test.bat`: runs a scenario, saves screenshots
  of every view); form ↔ scenario round trips are untested.
- Native engine: the ion/electron stop energies and the displacement
  energy default (`native.ion_stop_energy`, `electron_stop_energy`,
  `displacement_energy`) and the electron parameter tables can only be set
  in scenario files, not in the GUI.
- The native engine takes the ion's isotope from the mass rounded to a
  mass number (`isotope_mass_number`), so it uses the tabulated isotope
  mass, not the exact mass entered.
- `windeployqt` can be blocked by application-control policies (Device
  Guard); the post-build step then only warns and the GUI needs the Qt DLLs
  on PATH (`build\Release\generators\conanrun.bat`).
- The automatic-target-depth controls, the "ions escaped" question and the
  setup-warning indicator were added without a run of the GUI smoke test
  (they build, and the engine side is unit-tested; the dialogs are not).
- Stop during the automatic target depth takes effect between test runs:
  a test run of `auto_depth.test_ions` ions always finishes (typically well
  under a second; seconds for very energetic heavy ions).

### Build, CI and repository

- Release flags include `-march=native` (non-portable binaries) and are
  applied twice (CMAKE_CXX_FLAGS_RELEASE and the `compiler_flags` target).
  `ABSL_MIN_LOG_LEVEL=2` strips INFO/WARNING logs in Release.
- Test directories differ: `qbuild.sh` uses `build/Release/test`,
  `qbuild.ps1` and CI use `build/test`, `build_and_test.bat` uses
  `ctest --preset conan-release`. The batch file is the verified path on
  Windows.
- CI (`.github/workflows/tests.yml`) calls Visual Studio 2022's dev shell
  on a self-hosted runner; the workstation has Visual Studio 18. There is
  no Linux test job.
- `compile_commands.json` and `.vscode/.conan_tools/` are tracked although
  `.gitignore` lists them; `.vscode/.conan_tools/.devcontainer/devcontainer.json`
  declares `"mounts"` twice (the second wins). `.conan/windows-gcc` pins
  `compiler.version=16.2`.
- `test/common/common_test.cpp` relies on `#define private public` to read
  private data members (undefined behaviour; MSVC also encodes access in
  mangled names, so it must never be used to call private member
  functions).
- Tests: `CommonTest.electron_bombardment_runs` and `CommonTest.ion_test`
  are skipped; `CommonTest.electronic_stopping_energy_test` checks almost
  nothing. `StoppingTable.RealSrimTablesIfAvailable` and
  `TrimRegression.AgainstTrimWithSrimTablesIfAvailable` still skip unless
  `AMCSET_SRIM_STOPPING_TABLE` is set, although the bundled tables are now
  always available. The native engine's transport results (ranges, damage)
  are not validated.
- Qt 6.11.1 is pinned exactly; `conanfile.py` keeps Qt 3D and Qt Graphs
  enabled only so that the existing Conan package is reused (the GUI no
  longer uses them).
- File names are kept as `std::string` and opened as narrow paths: on
  Windows a path with characters outside the ANSI code page (e.g. a user
  folder in another script) cannot be opened. Fix: hold paths as
  `std::filesystem::path` from the front ends down (Qt `QString` →
  `toStdWString()`), or UTF-8 code page in the application manifest.
- The 42 MB `data/srim2013/srim2013_stopping.amcset` is a plain file in git
  (below GitHub's 50 MB warning size); Git LFS would keep clones small.
- `common/src/data_files.cpp` and `data_files.h` lacked the license header
  until 2026-11-23 (CI's header check would have failed); fixed.

---

## 3. SRIM mode

Details, numbers and diagnostics: `docs/SRIM_MODE.md` §6–8 and
`docs/TRIM2013_ALGORITHM.md`.

- **[Medium, open] Replacements in Xe gas.** `gas_xe_xe_1000keV_fc_gas`
  (Xe → Xe gas, 1 MeV, full cascades) has 4.6 % more replacements than
  TRIM's single 500-ion run (z = 5.0), with vacancies +0.7 % and recoil
  ionization −1.6 %; every other comparison (1031 of 1032) is within 3σ.
  Ruled out: TRIM's recoil-stack limit, Tmagic's iteration limits, the
  zero-energy defaults. Suspect: TRIM's single-precision recoil slots at
  depths of 10⁶ Å (not reproduced). Next steps: emulate single precision for
  the slot coordinates and energies behind an option and compare; a second
  TRIM run of this input to measure TRIM's own scatter (8 h).
- [Low] TRIM's recoil-stack limit (Data/TRIM.cfg, 200 000 in the reference
  installation) and its "Minor TRIM Error" message box are not reproduced;
  AMCSET's stack is unlimited. Only heavy-ion gas cascades come near it.
- [Low] `physics.transport: refined` (AMCSET's earlier transport) keeps the
  old deviations from TRIM; use it only where TRIM's coarse steps are not
  wanted.
- [Low] The built-in stopping model is within ~10 % of SRIM for H, ~25 %
  for He, but up to 2× high for slow heavy ions (Fe in Fe at 10–100 keV).
  The bundled SRIM-2013 tables are the default; without them a run warns
  and falls back to the built-in model.
- [Low] The performance figures in `docs/SRIM_MODE.md` §8.0 come from the
  validation runs (TRIM ran ~14 instances at once, AMCSET 0.1.0); the
  controlled benchmark (`tools\windows\benchmark_vs_trim.bat`) has not been
  run yet. (The speed-ups of §8.1-8.2 against 0.1.0 were measured on the
  idle workstation.)
- [Low] The fast approximations (`physics.fast_math`) agree with the exact
  calculation (`docs/SRIM_MODE.md` §8.2) but have not been compared with the
  TRIM reference runs themselves (`compare_with_trim.bat` with `--fast`).
- [Low] Trajectory and MD event files (`docs/OUTPUT.md`): times are an
  estimate (flight length / speed; TRIM has none); the files hold one
  run()'s ions (GUI *Continue* replaces them); only the trim2013 transport;
  `parent_event` is always -1 (no event descends from another by
  construction). The Parquet files were checked with an independent
  decoder of the format, not yet with pyarrow on the workstation.
- **[Medium] The crystal transport's stopping division is not
  calibrated.** `physics.transport: "crystal"` follows the lattice
  (docs/CLI.md "Crystals"). Its trends are checked, and a measured
  channeling dip (He in Ge) and measured boron profiles in Si are
  reproduced in peak, width and half-angle, but those data only bracket
  `physics.crystal.local_stopping_fraction` (0.25 to 1; default 0.5): the
  depth of a channeling tail carries that uncertainty (about 10 % for
  30 keV B in Si, up to 50 % for 0.5 keV B). Other species and targets are
  not compared. With TRIM's 1 eV cut-off its full-cascade vacancy count is
  11 % above TRIM's in the one case compared (not understood). Missing:
  damage build-up, TRIM.DAT sources, COLLISON.txt/trajectory/MD/dpa
  outputs. Polycrystals run from a scenario but have no GUI, and their
  Neper and EBSD readers are tested on hand-written files only: check the
  orientation convention with a known grain (docs/CLI.md "Polycrystals",
  Limits; what is left: docs/CRYSTAL_DESIGN.md §13). The GUI's Crystal tab has no
  automated test: it was clicked through by hand once (2026-11-24: the atom
  menu, the periodic table, sharing and occupancy, undo, turning and
  zooming, typed view and beam directions, the presets, the structure
  dialog, CIF export and import, the site table, bad input, a run), which
  left out Shift+drag (rolling the view), the prompt for the elements of a
  POSCAR or LAMMPS file, and cells of thousands of atoms. The other
  transports ignore crystals and say so. Next phases:
  `docs/CRYSTAL_DESIGN.md` §8.
- [Low] Cell files: a CIF must list its symmetry operations (space-group
  names alone are refused); only the first data block is read; no built-in
  lattice constants or Debye temperatures; the Debye amplitudes have not
  been compared with tabulated values (`docs/CRYSTAL_DESIGN.md` §10).
- [Low] Results of the same seed are identical on every thread count and
  run, but not between Windows (MSVC, Microsoft's C runtime) and Linux
  (GCC, glibc): `exp`, `log` and `pow` differ in the last bits.
- [Low] The automatic target depth sizes the last layer from the deepest of
  `auto_depth.test_ions` test ions; a much longer calculation can still
  transmit a few ions of the far tail (reported during the run). It does not
  resize the depth window except to reset one that ended beyond the target.
- Fixed in 0.1.0 (trim2013 transport, found by disassembling TRIM.exe):
  sputtering yields 4–10× below TRIM's, H ions with 2–10 % fewer vacancies,
  the ions/recoils split of IONIZ and PHONON, the sputtered-atom spectra.
  TRIM.DAT inputs (codes 4–7) are supported; the stopping tables are bundled
  uncompressed.
