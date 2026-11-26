# Changelog

## Unreleased

### Added

- **Random seed for every run.** `"seed": "random"` (CLI `--seed random`,
  GUI *New random seed for every run*): the engine draws a seed (system
  entropy, clock and a counter; 31 bits) when it is built, so every run
  differs. The seed used is printed (CLI), logged (GUI) and recorded in
  TDATA.txt, summary.json (`"seed"`, plus `"seed_random": true`) and the
  run table; entering it as the seed repeats the run. A continued GUI
  calculation keeps its seed. Tests `TrimEngine.RandomSeedIsDrawnRecordedAndRepeatable`,
  `Cli.RandomSeedKeyAndOption`.
- **A recoil atom (PKA) inside the target, in the CLI and the GUI.**
  `source.repeat` (new key) runs every particle of a source N times, each
  with its own random numbers; `--recoil EL,E,DEPTH[,CX,CY,CZ]` sets up one
  recoil atom run `-n` times, and the GUI's Ion tab has *One recoil atom in
  the target (PKA)* (element, energy, depth, direction, number of runs),
  whose cascades the path views show. The crystal transport now runs
  particle sources (before: refused): TRIM.DAT ions start at their points,
  a recoil atom at the nearest lattice atom of its element, leaving its site
  empty (`LatticeSearch::nearest_atom`). docs/CLI.md "A recoil inside the
  target". Tests `TrimDat.RepeatRunsEachParticleSeveralTimes`,
  `TrimCrystal.RecoilStartsFromItsSite`, `Cli.RecoilOptionAndRepeatKey`.

### Changed

- **Automatic target depth shrinks to the ions and rounds.** The sized
  target is accepted when the deepest test ion reaches
  `auto_depth.shrink_below` of it (new key, default 0.7, GUI *Shrink
  below*; before: the last 10 %, i.e. 0.9); shallower, the last layer is
  made thinner so that the deepest stops at 1 − margin, now also when the
  ions stop above the last layer (it is made about a monolayer thin;
  before, it was left as it was). Every thickness tried is rounded up to
  at most one decimal in the unit of its size. 5 MeV Fe → Fe, 1000 test
  ions: 2.07 µm in 2 test runs before, now 2.4 µm in 1 (the first guess,
  deepest ion at 78 %); bcc Fe (1 1 1), crystal transport: 4.29 µm → 4.3
  µm. Sized thicknesses therefore differ from 1.0.0's. Tests `TrimAutoDepth.*` (new:
  `ShrinkThresholdDecides`, `IonsStoppingAboveTheLastLayerMakeItThin`).

### Fixed

- **The refined transport refuses particle sources.** It started every
  particle at the surface as a beam ion, so a TRIM.DAT source silently ran
  the beam; validation now says to use `trim2013` or `crystal`.

- **MD events in the crystal transport: missing crossings** (a user's
  report). The steps with no atom within reach (`crystal_follow`'s free
  steps, most of a channeled path) and the flight in which a mover stops
  lowered its energy without the crossing check, so a level whose E_c was
  passed there got no `crossing` event and kept treating the atom as the
  BCA's: its later births were written as that level's `birth` events and
  the event sets of the levels no longer nested. In 40 ions of 5 MeV Fe →
  bcc Fe (1 1 1) with E_c = 0.5 … 50 keV, 1 to 8 ions per level from
  1 keV up had no projectile crossing (21 in all) and 694 events broke the
  nesting of adjacent levels; now none. The check now
  follows every energy loss. Record mode: the events change, no tally or
  random number does; handoff mode: such atoms are now handed to MD where
  they cross. The `trim2013` and refined transports are unchanged. Test
  `TrimMdEvents.CrystalCrossesEveryLevelOnceAndLevelsNest`.

## 1.0.0 — November 2026

### Release 1.0: packages, signing, platform integration

- **Version 1.0.0**; a logo (`resources/icons/`, drawn by
  `scripts/make_icons.py`): the window icon, the executables' icon, the
  packages' icons.
- **Packages** from `cmake --install` and CPack (`packaging/`): Windows MSI
  (WiX; Program Files, Start menu, `amcset` on the PATH), ZIP and MSIX
  (`tools/windows/make_msix.ps1`); Linux AppImage, `.deb` (Debian 13), `.rpm`
  (Fedora 42) and a CLI tarball for clusters. Install rules for both
  systems (`bin/`, `share/amcset/`, `share/doc/amcset/`; on Linux a desktop
  entry, AppStream metadata, hicolor icons and the man pages `amcset(1)`,
  `amcset-gui(1)`). docs/INSTALL.md.
- **Release workflow** `.github/workflows/release.yml`: builds and tests on
  Windows, Ubuntu 22.04, Debian 13 and Fedora 42 on GitHub's runners,
  signs the Windows files when Azure Artifact Signing is set up
  (`tools/windows/sign.ps1`), writes SHA256SUMS.txt and build attestations,
  and drafts the GitHub release. `conanfile.py` takes `-o with_qt=False`
  (Qt from Qt's own binaries). 1.0 is not signed: the release notes say so
  and docs/INSTALL.md "Running the unsigned Windows release" tells users
  what SmartScreen and Smart App Control do and what they can do.
  docs/RELEASING.md, including the signing
  options and what the maintainer has to set up.
- **Windows integration**: version information and icon in both
  executables; an application manifest (Windows 10/11, Per-Monitor v2 DPI,
  long paths, UTF-8 code page: file names outside the ANSI code page now
  work in scenarios and on the command line). docs/PLATFORM.md lists
  Microsoft's desktop-app requirements and how each is met.
- **Signals and shutdown.** The CLI stops a run on SIGINT, SIGTERM and
  SIGHUP (Windows: Ctrl+C, Ctrl+Break, closing the console, log-off,
  shutdown), writes the partial results, then ends by the same signal
  (status 128 + n; before: only Ctrl+C, exit 130); a second interrupt ends
  it at once. The GUI quits quietly on SIGTERM, SIGHUP and SIGINT (Linux),
  and when the desktop session ends it stops a calculation and keeps a
  changed scenario in a recovery folder instead of asking (Windows'
  rule: never block a log-off). The GUI parses numbers in the C locale
  whatever `LC_NUMERIC` says (Linux).
- **Validation study, docs/VALIDATION.md**: every comparison in one
  document (moved there: the TRIM.exe comparison and the damage models
  from docs/SRIM_MODE.md §6, the crystal measurements from docs/CLI.md),
  and new comparisons with measurements (`tools/validation/`): a scan of
  the SRIM-2013 stopping tables over every ion (Ni up to 39 % below its
  neighbours Co and Cu, Au 26 % below Pt and Hg, in every target), Ni
  ranges in Fe and Mo (Gurbich et al. 2025), 50 keV Sb in Si by SIMS (Titze
  et al. 2021: AMCSET, like SRIM, 22 % deep), and the ionization yield of
  keV Ge and Si recoils (25-32 % low in Ge, 28 % to 3 × high in Si; no
  stopping model fixes it). Further: electronic stopping against the
  38 268 measurements in elemental targets of the IAEA database (median
  0.996, 81 % within 10 %; `stopping_reference.py`) and NIST PSTAR/ASTAR
  (within 3 % above 0.1 MeV); measured Ni peak depths at 2-11.5 MeV in Ni,
  Fe and steels (Chauhan et al. 2026, Pena et al. 2025, Song et al. 2026:
  SRIM-2013 1.2-1.9 × deep, the correction brings six of ten within 9 %);
  350 measured ranges in Si, SiO2, C, B and SiC (Behar et al. 1995,
  Oosterhoff 1988, Zhang et al. 2009: agreement in Si, heavy ions 12-25 %
  short in light targets, as in SRIM); sputtering yields of Ar, H, D and He
  (Yamamura et al.: median 1.05, pair factors 0.68-1.9); published SRIM
  damage calculations (Stoller et al. 2013, Weber and Zhang 2019); and the
  crystal transport's channeling tails of 15 keV As in Si against tilt
  (Tasch et al. 1994: tails 9-20 × too high at 100 nm, perfect crystal).
  All of it, with the TRIM comparison and the damage models, is explained as
  a review with 35 figures in `docs/validation/` (nine chapters);
  `tools/validation/make_figures.py` (needs matplotlib) redraws them from
  the scripts' output, and `literature_cases.py` now also writes
  `results.json`.
- **Ni stopping correction** `physics.stopping_corrections.nickel_ions`
  (CLI `--correct-nickel`, also for `amcset stopping`; GUI Calculation tab;
  off by default): Ni projectiles and recoils take the geometric mean of
  Co's and Cu's SRIM-2013 electronic stopping at the same velocity. It
  shortens 5 MeV Ni in Mo by 1.275 × (measured: SRIM-2013 1.3 × too long)
  and 11.5 MeV Ni in Fe by 1.273 × (1.5 ×); measured peak depths at
  3.5-9 MeV then agree within 9 %, 2 MeV stays 1.2-1.8 × deep
  (docs/VALIDATION.md §4.1). With it off nothing changes.
  New TRIM.exe reference runs of Ni ions (`tools/srim/gen_cases4.py`, 11
  runs; `run_trim_cases3.ps1 -CaseSet/-OutSet`): AMCSET as TRIM agrees
  (166 of 169 quantities within 3σ, docs/VALIDATION.md §2.3); with the
  correction Ni ranges are 21 % shorter at 5-11.5 MeV, Co and Cu unchanged.
  `StoppingTable.NickelCorrectionIsTheNeighboursGeometricMean`,
  `Cli.NickelStoppingCorrection`.
- **Automatic target depth with the crystal transport** (it was refused):
  the first test run is made in a box 32 times as deep as the estimate,
  which knows no lattice, and the box is then trimmed to the deepest
  channeled ion (two test runs as a rule; docs/CLI.md "Checks and
  warnings"). `TrimCrystal.AutomaticDepthHoldsTheChanneledIons`.
- **GUI**: *Turn the crystal to this view* (Crystal tab) sets the crystal's
  orientation from the view of the unit cell, as *Point the beam along this
  view* sets the beam (`crystal::orientation_for_frame`: Miller indices
  where a low-index plane fits, else Euler angles;
  `CrystalOrientation.OrientationForAFrame`). The Output tab has switches
  *Write a trajectory file* and *Write an MD event file*; before, a file was
  switched on only by typing its name, and the group's tooltip still said
  "TRIM-2013 transport only".
- **Native engine, electrons**: the bundled electron tables are installed
  (`share/amcset/electrons/`) and found from any working directory; before,
  only `./common/assets` was tried, so electron runs failed outside the
  source tree. `CommonTest.electron_bombardment_runs` is no longer skipped.
- 247 CTest tests on Windows (246 on Linux):
  `Cli.InterruptsSetTheFlagAndGiveTheSignalsStatus`,
  `TrimCrystal.HeavySlowIonsKeepTheAmorphousRange`,
  `TrimCrystal.AutomaticDepthHoldsTheChanneledIons`,
  `CrystalOrientation.OrientationForAFrame`,
  `StoppingTable.NickelCorrectionIsTheNeighboursGeometricMean`,
  `Cli.NickelStoppingCorrection`,
  `DataFiles.ElectronTablesAreFoundOutsideTheSourceTree`.

### Since 0.1.0

A stability and performance sweep. The default (exact) SRIM-mode results are
bit-identical to 0.1.0 for every reference case (`tools/perf/bench.py
compare`: 377 of 377 output files before the dpa table was added), except
for the two fixes under "Behaviour changes". The dpa table adds `csv/dpa.csv`
and a `dpa` table and `damage_models` block to `results.json` and
`summary.json`; every TRIM-format file and every other table stays
identical (351 of 377 files identical, the other 26 are those two JSON files
of the 13 cases, differing only by the additions).

### Build

`-DAMCSET_MARCH=...` replaces the `-march=native` of GCC/Clang Release builds
(default unchanged); portable builds (the AppImage, the CLI tarball) use
`x86-64-v2`.

### Behaviour changes

- **Sputtering from lower layers** (SRIM mode, multilayer targets): an atom
  that reaches the surface from a deeper layer was judged by the surface
  binding energy at its element index in the *current* layer's list, i.e.
  another element's Es, or memory past the end of a shorter list. It now
  uses its own Es (`TrimEngine.SputteredAtomsUseTheirOwnSurfaceBindingEnergy`).
  No reference case changes (the layered case is identical).
- **Scattering at huge impact parameters**: where the screened potential
  underflows (reduced impact parameters above ~3500: only in solids of gas
  density), the magic formula returned NaN, which stopped the run; it now
  returns no deflection (`Scattering.ExtremeImpactParametersGiveFiniteAngles`,
  `TrimEngine.SolidOfGasDensityRuns`). Condensed targets never reach it.

### Performance

- SRIM-mode transport 1.35× faster on one thread (geometric mean of 13
  reference cases, 1.29-1.40×; whole process 1.67×) with identical
  results: the screening function and its derivative share their
  exponentials, the GUI's live tallies are no longer accumulated twice per
  event.
- Start-up ~0.8 s → ~0.05 s: the 42 MB stopping file is read in one piece
  and each table parsed on first use. Output files are written in parallel.
  The CLI no longer waits up to 200 ms after a run (a condition variable
  instead of polling). Fixed cost per run: 0.93 s → about 0.14 s.
- Output numbers through `std::to_chars` into large buffers instead of
  `snprintf` and a `std::string` each (identical files): for 2 million ions
  the writing after the transport fell from 2.4 s to 1.0 s. The depth-bin
  edges are precomputed and the random-number generator reduces without a
  division (identical results).
- **Fast approximations** (`physics.fast_math`, CLI `--fast`, GUI check box
  with a disclaimer; off by default): tabulated magic formula (relative
  error < 4e-5), polynomial sine/cosine of the azimuth, rotations with one
  reciprocal, cheaper powers; then also the stopping powers, the
  Kinchin-Pease partition and the flight path from bit-indexed tables (no
  log/exp/pow/cbrt per lookup). A further 2.08× with MSVC before the tables
  (2.81× over 0.1.0; GCC 15 another 14 %);
  statistically the same results as the exact mode (421 quantities,
  median |z| 0.65; `tools/perf/fast_vs_exact.py`), but not number for number,
  so not a step-by-step reproduction of TRIM.exe (`docs/SRIM_MODE.md` §8.2).
- Measured and not kept: `/arch:AVX2` (no gain), PGO (+2 %; the option
  `AMCSET_PGO` stays), polynomial `exp`/`log` (slower). Release builds now
  carry PDBs (`/Zi`, `/DEBUG /OPT:REF /OPT:ICF`).

### Stability and diagnostics

- **Crash reports** (`diagnostics.h`): on a crash (access violation, stack
  overflow, abort, an exception escaping a thread) the CLI and the GUI write
  the symbolised stack, the calculation, what each worker was transporting
  and the options that rerun that ion (`first_ion`, new scenario key), the
  last log lines, and on Windows a minidump; exit code 3. `amcset
  diag-crash ...` tests the chain. Diagnostic log: CLI `--log FILE`, GUI log
  folder (Help menu) with a detailed-log toggle.
- **Exceptions**: none reaches Qt (every slot through `guarded()`, an
  `Application::notify` override, overrides catch); worker exceptions name
  the ion and seed; the native engine's thread pool and the parallel output
  writer rethrow after all tasks finished (an escaping exception used to
  terminate the program); a failed thread start is reported, not fatal.
- NaN guards in the transport (a NaN energy or depth used to loop forever);
  the error now prints the last flight's values.
- `FormBinder` no longer stays stuck ignoring edits after an exception;
  rejected edits are reported; `std::clamp` with inverted bounds fixed.

### Input checking

- Physical limits with messages that suggest what to use (`docs/CLI.md`
  "Checks and warnings"): ion energy above the stop energy and at most
  100 GeV, isotope masses, angles below 90°, layers at least one atomic layer
  thick, densities up to 50 g/cm3, Ed up to 10 keV, Eb and Es up to 100 eV,
  stop energies up to 1 keV, bins, threads, existing table files and SRIM
  folders. The GUI refuses such values as they are entered.
- Malformed numbers are errors in scenario JSON, TRIM.IN and TRIM.DAT (no
  truncation, wrap-around of negative or huge counts, NaN or infinity); a
  negative element count in a TRIM.IN used to index past the arrays. Output
  file names cannot leave the output folder.
- **Plausibility warnings** (`validate`, before a run, the GUI status bar):
  densities far from the material's, gas flags that do not fit the density,
  targets far thicker than the ions travel (with a suggested thickness),
  targets that transmit most ions, relativistic ions.
- **Automatic target depth** (`auto_depth`, CLI `--auto-depth`, GUI
  *Target depth* group; off by default): quick test runs size the last layer
  so that the deepest of 1000 test ions stops in the last 10 % of the
  target (first guess from the stopping powers; 2 test runs, ~0.3 s, for
  every case tried from H to U and 1 keV to 1 GeV). During a run, ions
  transmitted through the back are reported (the run continues); the GUI
  offers to rerun with the automatic depth.

### Trajectory and MD event files (opt-in; not TRIM)

- **Trajectory file** (`output.trajectory_file`, CLI `--trajectory-file`,
  GUI Output tab): every recorded point to CSV or Parquet, in ion order for
  any number of threads, with ion, point index, track, parent track and
  point, kind (birth, collision, stop, exit), generation, Z, mass, position,
  direction, energy, energy before the collision, energy transferred,
  electronic loss and time; `points: births` and `min_energy` filters.
  `output.trajectories: "all"`. Trajectory points now carry these fields in
  the library too (`common::CollisionPoint`).
- **Time** along trajectories (flight length / speed, from the ion's entry;
  an estimate, TRIM has none), filled when trajectories or events are
  recorded.
- **Hybrid BCA-MD events** (`output.md_events`, CLI `--md-events`,
  `--md-ec`, `--md-emin`, `--md-handoff`): recoils born below a transition
  energy E_c and movers falling below it, for a list of E_c in one run;
  record mode (transport and tallies unchanged, tested) or handoff mode
  (the BCA stops following events; energy conserved to 10⁻³ with the TRIM
  bugs fixed). Totals in summary.json, the CLI output and the table
  `md_events_summary`. 5 MeV Fe → Fe, E_c = 1 keV: 3 090 events and 0.98 MeV
  to MD per ion.
- `--fix-trim-bug(s)` are applied before `--set`, so `--set
  physics.trim_bugs.NAME=true` can keep one of TRIM's behaviours (needed for
  Eb = 0 with the other bugs fixed).
- **`amcset stopping --lammps`**: the table for LAMMPS `fix electron/stopping`
  (metal units, one column per atom type, `--ntypes`).
- Streaming CSV/Parquet writer (`output/stream.h`; Parquet with several row
  groups) for outputs larger than the memory.
- Fixed: a target atom that never moves (TRIM's stale-energy primary with
  energy 0) no longer leaves trajectory points with zero energy and no
  direction; a recoil reflected at the surface is recorded where it is
  reflected. (Trajectories only; no tally changes.)

### Crystal transport: two corrections; why ranges in crystals are longer

Both corrections change crystal-transport results only (the other
transports are untouched); both made ranges too short.

- **Electronic straggling** of the crystal transport is now always a gamma
  deviate. Its flights are the ~1 Å steps between atoms, whose straggling
  exceeds the mean loss; the Gaussian cut at zero (the default
  `straggling_distribution`) raised the mean loss: 100 keV He in a random
  direction of Si came out 8 % short (0.92 of the amorphous range; now
  1.02). `straggling_distribution` now applies to `refined` only.
- **The path is charged once.** A flight was charged to the struck atoms'
  closest approach, but the mover turns x₁ = s tan(θ/2) before it and the
  next flight started there, charging that stretch again; in amorphous
  layers the next partner was also placed one monolayer from the deflection
  point instead of from the last closest approach. Heavy slow ions, whose x₁
  is a few per cent of a monolayer, had ranges 2-5 % short in amorphous
  layers: 30 keV W in amorphous W 50.3 Å against `trim2013`'s 52.9 Å
  (20000 ions; now 52.7 Å), Cu 106.4 against 109.0 Å (now 108.8 Å)
  (`TrimCrystal.HeavySlowIonsKeepTheAmorphousRange`).
- `tools/crystal/range_study.py`: ranges of B, As, He in Si, Fe in Fe, Cu in
  Cu, W and Xe in W along axes, planes, the 7°/22° "random" direction, 16
  random orientations, a static lattice, tilt scans, and 0.5 keV B against
  measured profiles in crystalline and preamorphized Si. The longer ranges
  in crystals are channeling: the transport without a lattice gives
  `trim2013`'s ranges (0.98-1.02), axes give 1.5 to 16 times the amorphous
  range, slow heavy ions channel even at 7° (their critical angle, ψ₂, is
  7.5-11°), and the measured profiles show the same. The study and its
  numbers: docs/CLI.md "Why ranges in a crystal are longer". summary.json
  and the `crystals` table now give the critical angle (ψ₂ below E′) besides
  ψ₁.
- The numbers of `channeling_checks.py`, `calibration.py` and
  `polycrystal_checks.py` in docs/CLI.md are from runs after both
  corrections.

### Trajectory and MD event files in every transport; crystal orientation output

- The trajectory file and the MD event file (with the cascade time) now
  work with the `refined` and `crystal` transports too, not only
  `trim2013`: tracks, parents, birth/collision/stop/exit points, times,
  record and handoff modes. Recording changes no tally (tested for each
  transport; the results with the files off are the same as before, number
  for number). docs/OUTPUT.md "Other transports and crystals".
- Crystal orientation in the outputs: summary.json `crystals` (rotation
  matrices to TRIM's axes and to the files' axes, Bunge angles, the crystal
  directions along the depth and the beam, the thermal amplitudes used, and
  the low-index axes nearest the beam with the atom spacing and Lindhard's
  ψ₁), a `crystals` result table, and in the trajectory and event files of
  the crystal transport the lattice offset of each ion (`lattice_dx_A`,
  `lattice_dy_A`) and the orientation in the Parquet metadata.
  docs/OUTPUT.md "Crystals".
- 240 CTest tests on Windows (239 on Linux):
  `TrimTrajectories.EveryTransportRecordsTreesAndEvents`,
  `TrimCrystal.OrientationIsReported`.

### Crystals, phase 3: polycrystals (scenario and transport; no GUI yet)

- **`layers.N.crystal.polycrystal`**: a layer of grains of one unit cell.
  Grains are generated (a Poisson-Voronoi tessellation with random, fiber
  or mosaic orientations, optionally columns) or read from AMCSET's grain
  file (seeds or voxels), a Neper `.tess` or an EBSD map (`.ang`, `.ctf`);
  `boundary_width` makes a disordered slab at every grain boundary.
  docs/CLI.md "Polycrystals".
- The crystal transport follows the lattice of the grain a moving atom is
  in (`crystal::GrainMap`, one lattice per grain). 30 keV B into Si along
  ⟨001⟩: 2641 A in the single crystal, 1136-1253 A in random grains of
  50-3000 A, 1743 A with a ⟨001⟩ fiber texture, 1120 A amorphous
  (`tools/crystal/polycrystal_checks.py`).
- Single-crystal results are unchanged; `trim2013` and `refined` are
  bit-identical (`bench.py compare`: 390 of 390 files).
- Not yet: the GUI for polycrystals, a check of the Neper and EBSD readers
  against real files (docs/CRYSTAL_DESIGN.md §13).
- 238 CTest tests on Windows (237 on Linux): `CrystalGrains` (7),
  `TrimCrystal.PolycrystalGrains`, `RunConfig.PolycrystalKeys`.

### Crystals: common crystals, measurements, the GUI's Crystal tab

- **Common crystals** (`crystal/materials.h`): Tables IV and V of Gemmell,
  Rev. Mod. Phys. 46, 129 (1974): structure and room-temperature lattice
  constant of 19 cubic elements and 24 cubic compounds, and the Debye
  temperatures of the elements. **Behaviour change of the crystal
  transport:** a crystal of one element that the table lists now vibrates
  with its tabulated Debye temperature when the scenario gives none (it
  stood still before); the run says so. Compounds and other elements are
  unchanged (static unless given).
- **Against measurements** (`tools/crystal/calibration.py`; docs/CLI.md
  "Against measurements"): the Debye model reproduces u₁(293 K) of all 19
  elements of Table IV; the dip of 1 MeV He along ⟨110⟩ of Ge has a
  half-angle of 0.99° (measured 0.95°) and a minimum yield of 0.022
  (measured 0.03); B in (100) Si at 15 keV: Rp 62.2 nm and ΔRp 33.1 nm
  against the measured 64.2 and 31.8 (the amorphous calculation: 60.5 and
  25.7); the 30 keV tail and the first decade of the 0.5 keV profile agree
  to a few per cent. The division of the electronic stopping
  (`physics.crystal.local_stopping_fraction`) is bracketed by these data
  (0.25 to 1), not fitted; the default stays 0.5 and the warning of the
  transport says what is established.
- **GUI: the Crystal tab.** Per layer: amorphous or single crystal; the
  unit cell from 43 common crystals, a built-in structure (sublattice
  elements, mixtures, interstitials, lattice constant or density) or a CIF,
  POSCAR, LAMMPS or XYZ file, with export; lattice parameters and a table
  of sites; a 3D view of the cell in which clicking an atom replaces its
  element, shares the site, sets the occupancy, changes the kind of site or
  removes it; the view direction can be made the beam direction (and a
  typed beam direction turns the view); orientation by Miller indices or
  Euler angles, twist; temperature and thermal vibrations; the settings of
  the crystal transport and the beam spot. The Target tab shows a crystal
  layer's density, elements and stoichiometry read-only.
- **The crystal transport takes none of TRIM.exe's bugs**, whatever
  `physics.trim_bugs` says (`TrimConfig::effective_bugs`): its results never
  depended on the switches, but results.json recorded them as reproduced
  and the GUI showed them ticked. A crystal run now records them as
  corrected, and its TRIM-format count tables hold the exact sums
  (`output.trim_single_precision` does not apply; this changes VACANCY.txt,
  NOVAC.txt and RANGE.txt of crystal runs by up to the 1.5 % of TRIM's
  rounding). The scenario keeps its switches for the other transports.
- **GUI: a switch for the crystal transport.** *Crystal transport (follows
  the lattice)* in the Calculation tab, and the same box at the top of the
  Crystal tab, turn it on and off; the transport list holds TRIM-2013 and
  Refined and applies again when the switch is off. With it on, the
  *TRIM.exe bugs* group shows every bug corrected and says why.
- 229 CTest tests on Windows (228 on Linux): `CrystalMaterials`,
  `TrimCrystal.TakesNoneOfTrimsBugs`,
  `RunConfig.CrystalRunsRecordEveryTrimBugAsFixed`. `trim2013` and
  `refined` are untouched (`bench.py compare`: 390 of 390 files identical).

### Crystals, phase 2: the crystal transport (not TRIM; not calibrated)

- **`physics.transport: "crystal"`** follows the lattice of layers that have
  a crystal, so that ions can channel: the binary-collision scheme of
  Robinson and Torrens, Phys. Rev. B 9, 5008 (1974) (collision partners on
  lattice sites within a maximum impact parameter, Gaussian thermal
  displacements, simultaneous collisions added vectorially, replacement
  when the projectile keeps less than min(cut-off, Ed)), with AMCSET's
  potential and stopping powers; the electronic stopping is divided between
  the path and the collisions (`physics.crystal.local_stopping_fraction`),
  normalised so that random directions lose the full stopping power.
  Amorphous layers take random partners. docs/CLI.md "Crystals",
  docs/CRYSTAL_DESIGN.md §11 (what follows the paper and what does not).
- **`ion.spot`**: where the ions meet the surface (rectangle, circle,
  Gaussian or a point); `physics.crystal` (maximum impact parameter,
  simultaneity distance, stopping division, thermal vibrations on/off).
- **Checked, not calibrated** (`tools/crystal/channeling_checks.py`):
  30 keV B in Si goes 2838 Å along ⟨001⟩, 4407 Å along ⟨110⟩ and 1150 Å at
  7° tilt / 22° twist, against 1128 Å in the amorphous target; hotter
  crystals and a covering oxide dechannel; a fast proton loses the stopping
  power to 0.2 % however it is divided; the channeling dip of 1 MeV He in
  W ⟨001⟩ has a minimum of 2 % and a half-width of 0.82 ψ₁; with the
  cut-off at Ed full cascades give 11 % more vacancies than NRT. Nothing is
  compared with measurements yet: absolute channeling ranges depend on the
  uncalibrated stopping division.
- The `trim2013` and `refined` transports are untouched (`bench.py
  compare`: 390 of 390 files identical) and still treat crystals as
  amorphous, with a warning.
- `--set` takes a list where the key holds one
  (`layers.0.crystal.orientation.surface=[1,1,0]`); the GUI's transport list
  has the crystal transport.
- 224 CTest tests on Windows (223 on Linux): `CrystalSearch` (the partner
  search against a brute-force search, thermal displacements and occupancy
  statistics, the largest hole of each structure), `TrimCrystal` (thread
  independence, energy conservation and defect bookkeeping, random
  direction = amorphous, channels, temperature, the He/W dip against
  Lindhard's angle, the beam spot, amorphous layers, alloys, interstitial
  atoms and sputtering, the cut-off at Ed against NRT),
  `RunConfig.CrystalTransportKeys`.

### Crystals, phase 1: the description

- A layer can carry a **crystal** (`layers.N.crystal`): a unit cell from a
  built-in structure (sc, bcc, fcc, hcp, diamond, zinc blende, wurtzite,
  rock salt, CsCl, fluorite, perovskite, L1_2, L1_0; substitutional alloys
  by site fractions; octahedral/tetrahedral interstitial sites), from a
  cell file (CIF with symmetry operations, VASP POSCAR, LAMMPS data,
  extended XYZ) or written out; an orientation (surface (h k l) and a
  reference direction, four-index hexagonal notation, or Bunge Euler
  angles, plus a twist); a temperature and Debye temperatures or rms
  amplitudes. The crystal sets the layer's composition and density.
  docs/CLI.md "Crystals".
- The `trim2013` and `refined` transports treat a crystal layer as
  amorphous and warn about it; a test checks that it gives the amorphous
  layer's results number for number. (The lattice is followed by the
  crystal transport of phase 2, above.)
- `amcset cell`: summary of a unit cell (parameters, density, sites,
  neighbour shells, directions of the lab axes), conversion between the
  four formats, the scenario's `crystal` object (`--json`), and the tilt
  and twist that send the beam along a crystal direction (`--beam`).
- New `common/*/crystal/` library (no engine dependency); `temperature`
  quantities (`K`, `C`).
- Scenario files write number lists in the shortest form (`[0.5,0.0,1.0]`
  instead of `[5E-1,0E0,1E0]`); both forms read back the same.
- 209 CTest tests on Windows (208 on Linux): `CrystalGeometry`,
  `CrystalLattice`, `CrystalStructures` (coordination numbers and distances
  of every structure and interstitial set), `CrystalOrientation`,
  `CrystalVibrations`, `CrystalFiles`, `RunConfig.Crystal*`,
  `Cli.CellCommand`, `cli.cell*`. All earlier outputs are unchanged
  (`bench.py compare`: 390 of 390 files identical).

### dpa by damage model (SRIM mode; not TRIM)

- New table `dpa` (all formats, GUI plots *dpa* and *rpa*) and
  `summary.json` block `damage_models`: displacements per atom by depth at
  a fluence (`output.dpa_fluence`, default 1e15 ions/cm², GUI Output tab)
  from TRIM's vacancies and from four damage functions of every primary's
  damage energy: NRT, Kinchin-Pease, arc-dpa and rpa (Nordlund et al.,
  Nat. Commun. 9, 1084 (2018)). arc-dpa/rpa constants built in for Fe, Cu,
  Ni, Pd, Pt, W; other elements through
  `layers.N.elements.M.damage_efficiency`. trim2013 transport only.
  docs/OUTPUT.md "dpa".
- Checked against published SRIM-2013 post-processing (Mitsi et al. 2023:
  arc/NRT of 30 ion-target pairs within 0.43 % rms) and compared with MD
  and 6 K experiments (docs/VALIDATION.md §5): for surviving defects in
  metals arc-dpa is the most accurate; TRIM's full-cascade vacancies are
  about 2× NRT and 4-19× arc-dpa; TRIM's replacements are 21-275× below
  the mixing of MD (rpa).
- Tools: `tools/damage_models/pka_damage_functions.py` (damage functions
  per PKA and model), `arc_dpa_vs_srim_postprocessing.py`.
- 181 CTest tests on Windows (180 on Linux): `TrimDamageModels` (published
  constants, the formulas for one primary, table scaling with fluence and
  density, elements without constants) and
  `RunConfig.DamageEfficiencyAndDpaFluence`.

### Tests and tools

- 164 CTest tests on Windows (163 on Linux): crash reports, diagnostics,
  the fast approximations, the limits and malformed inputs, plausibility
  warnings, the automatic depth.
- `tools/perf/`: `bench.py` (speed and bit-identity of a change),
  `fast_vs_exact.py`, `fuzz_inputs.py` (2400 mutated inputs: no crash, no
  hang), `winsample.cpp` (sampling profiler for Windows).

## 0.1.0 — November 2026

### Physics fixes (native engine)

Each fix changes simulation results; see `docs/KNOWN_ISSUES.md` §0 and the
tests named there.

- P1 Scattering rotation axis perpendicular to the direction (was (y, x, 0)).
- P2 Centre-of-mass energy from the reduced mass (was the projectile mass).
- P3 Laboratory angle from the full CM angle (the magic formula gives Θ/2).
- P4 Screening-function coefficient 0.02817 (was 0.2817, φ(0) = 1.25).
- P5 Magic-formula G = 1/(γ(√(1+A²) − A)) (the SRIM book prints it inverted).
- P6 Electron model rewritten with Boost.Units: target rest energy, parameter
  table index Z − 1, cross-section/density units of the mean free path, Z/A
  in the collision stopping (Bethe with n_e = N Z and the Joy-Luo J),
  Bethe-Heitler bremsstrahlung, relativistic speeds.
- Number densities without the stray Avogadro factor; per-element
  displacement energies are used.

### New

- **SRIM mode**: a TRIM-compatible engine (`common/src/trim/`) validated
  against SRIM-2013 TRIM on 78 cases (Fe → Fe 1 keV - 5 MeV; H and He in
  transition metals from Ti to W; Cr, Co and Cu self-ions; Fe in Fe-Cr and
  316 SS; Ar → Cu); `docs/SRIM_MODE.md`. TRIM's inputs (TRIM.IN),
  damage models (quick Kinchin-Pease with TRIM's damage-energy threshold,
  full cascades, monolayer steps), per-element/per-layer Ed, Eb, Es, compound
  corrections, gas targets, all of TRIM's tallies and output files, plus
  recoil spectra and `summary.json`. Reproducible (seed) and independent of
  the thread count.
- **Stopping**: SRIM-2013 stopping tables as an option (generated locally by
  `tools/srim`), relativistic Bethe-Bloch (effective charge, density effect)
  for high-energy particles up to 100 GeV, a built-in Lindhard-Scharff /
  Bethe-Bloch model, TRIM's electronic straggling. Available to both
  engines.
- **Scattering**: SRIM-2013 pair potentials (POTCOEF03/SNUC03), the classical
  scattering integral (Gauss-Mehler) as an alternative to the magic formula.
- **CLI** `amcset`: run, trim, convert, stopping, elements, template,
  validate (`docs/CLI.md`); JSON scenario files with units, shared with the
  GUI.
- **GUI**: engine selection (SRIM mode / native), damage model, seed,
  stopping source, SRIM table and data paths; Ed, Eb and Es columns per
  element; File > Open/Save for scenarios (JSON or TRIM.IN); incident
  energies up to 100 GeV; SRIM-mode runs write TRIM-format files.
- **Element database** (`element_data.h`): TRIM's ATOMDATA values (masses,
  densities, heats of sublimation) and NIST mean excitation energies.
- **Units**: `amcset_units.h` with every quantity type, unit and constant
  used by the code.
- **Tests**: physics, SRIM mode, TRIM.IN, run configuration, CLI and
  end-to-end CLI tests (101 CTest tests); optional tests against SRIM data.
- **Tools**: `tools/trim_compare` (TRIM vs AMCSET comparison),
  `tools/srim` (SR Module extraction, TRIM batch runs),
  `tools/windows/build_and_test.bat`.

### SRIM mode follows TRIM.exe (`physics.transport: trim2013`, default)

TRIM.exe (SRIM-2013.00) was disassembled to find the causes of the remaining
differences; `docs/TRIM2013_ALGORITHM.md` records its algorithm. The new
default transport follows it step by step; the previous SRIM-mode transport
is kept as `physics.transport: refined`.

- SRIM mode uses the SRIM-2013 stopping tables by default: the tables of
  every ion in every element, solid and gas (16928), are in
  data/srim2013 and found next to the executable, in the source tree, or
  through AMCSET_DATA_DIR. TRIM.IN inputs, scenario files without a
  stopping model, the CLI and the GUI use them; without the file a run
  warns and uses the built-in stopping. A loaded table is kept for the
  process.
- Bugs of TRIM.exe can be corrected one by one (physics.trim_bugs,
  --fix-trim-bug/--fix-trim-bugs, GUI "TRIM.exe bugs" with tooltips on the
  effect of each correction): stale primary energy, energy created in slow
  recoil collisions, struck-atom direction, replacements counted as
  interstitials, sub-threshold sputtered atoms counted as vacancies, ion
  transfers counted twice, replayed cascade random numbers, zero-energy
  defaults that skip the first atom, rounded layer Z, EXYZ.txt's skipped
  first step, doubled TRIM.DAT recoil events. The default reproduces TRIM
  (docs/SRIM_MODE.md §7).
- Validation against 161 TRIM-2013 runs of 75 inputs (with the bundled
  tables, new gas cascades and TRIM.DAT inputs): 1031 of 1032 comparisons
  within 3σ; replacements in Xe gas at 1 MeV (+4.6 %, z = 5.0) are open
  (docs/VALIDATION.md §2.2).
- Performance (docs/SRIM_MODE.md §8): the validation cases took 4.7 h of the
  PC with TRIM.exe and 5.3 min with AMCSET (53×; 30× per case, geometric
  mean). tools/windows/benchmark_vs_trim.bat measures TRIM alone against
  AMCSET on one and on all threads.
- tools/windows/benchmark_vs_trim.bat: wall-clock time of TRIM.exe and
  AMCSET (one thread, all threads) on the same inputs.
- Ion flights: TRIM-85 free-flight-path rule with the layer's mean Z and M,
  no 5 % energy-loss limit, random flights in full-cascade mode, boundary
  shortening; collisions at the energy after the flight's electronic loss;
  TRIM's uniform straggling deviate.
- Recoils: monolayer steps, TRIM's slot stack, sub-threshold atoms near the
  surface set in motion, reflection at the surface barrier, TRIM's
  replacement and interstitial bookkeeping, TRIM's Kinchin-Pease constants.
- Atoms struck by recoils leave in the direction TRIM gives them (its DIRCOS
  call with sin x = 0); this brings sputtering yields to TRIM's (Ar -> Cu
  1 keV: 3.54 vs 3.55 atoms/ion; Fe -> Fe 1 keV: 1.98 vs 1.98).
- The mover's energy after a recoil collision is limited to >= 0 as in TRIM
  (recoil phonons near the surface now match TRIM's).
- TRIM's use of random numbers: one generator, the same draws in the same
  order, a random starting point within a monolayer, and each cascade run on
  a copy of the ion's generator (TRIM saves and restores the generator state
  around every cascade). This correlates consecutive cascades as in TRIM and
  removes a ~1 % excess of recoil energy deposition in surface-sputtering
  mode.
- No pair potentials by default (`physics.pair_potentials: none`): TRIM.exe
  never reads POTCOEF03.DAT.
- SPUTTER.txt lists sputtered atoms as TRIM does (energy inside the surface,
  direction before the barrier); the header says so.
- `tools/trim_compare/compare.py` compares sputtering yields and reads TRIM's
  surface-sputtering (cascade code 3) VACANCY.txt correctly.
- TRIM replaces a displacement, surface or lattice binding energy of 0 by
  20, 4 or 2 eV for every target atom except the first (TDATA.txt shows it);
  AMCSET does the same, in TRIM's atom order, and warns. With Eb = 0 in a
  TRIM.IN, alloys and layered targets ran with Eb = 2 eV for their second and
  later elements in TRIM: this was the cause of the 2-6 % excess of
  replacements and of the recoil-phonon differences in Fe-Cr, Cr/Fe and
  316 SS.
- Quick calculation: E2RECOIL loses Eb once per displaced atom, not ν Eb;
  E2RECOIL's second column is TRIM's "Absorbed by Recoils" (one column; per
  atom only for full cascades, none for surface sputtering).
- TRIM.exe sums its count tables (VACANCY, NOVAC, RANGE) in single precision,
  which shifts them by up to ±1.5 % from the counts (+0.6 % for Fe → Fe at
  2 MeV); AMCSET's TRIM-format files reproduce it
  (`output.trim_single_precision`), the CSV/JSON/Parquet tables are exact.
- TRIM's atom order is kept from a TRIM.IN (`trim_atom`), for the column
  order and the energy replacement above.
- Quick calculation: "displacements per ion" in results.json and the GUI
  equal the vacancies (it was the number of the ion's knock-ons, now
  `knock_ons_per_ion`); no replacements or interstitials rows.
- GUI: an unattended run (--run/--screenshots) logs the end of its event
  loop and leaves with code 10 + its code if closing takes over 60 s.
- The SR Module job list (tools/srim) extracts gas-phase tables for every
  element, as TRIM uses gas stopping for all elements of a gas layer (C and
  O in CO2): He -> CO2 2 MeV range now within 0.03 % of TRIM (was +1.9 %).
- tools/trim_compare: TRIM's backscattered/transmitted counts from
  TDATA.txt (RANGE.txt does not always have them); NOVAC.txt of very wide
  bins (TRIM prints six decimals) no longer counts as zero replacements;
  z of very large t without a floating-point overflow.
- TRIM.DAT recoils follow TRIM's codes as TRIM.exe runs them: code 6 is a
  full cascade from each atom (it was read as the quick calculation: half
  of TRIM's vacancies), code 7 moves each atom like an ion with monolayer
  steps (its own ionization and phonons are the ion's, as in TRIM); code 6
  also books the events in RANGE.txt as TRIM does.
- EXYZ.txt: as in TRIM, the record at the start sets the next one two
  energy steps down (TRIM does not record the first crossing after the
  start), and TRIM's column widths.
- A sputtered atom set in motion below Ed now counts as a displacement as
  well as a vacancy in the totals, as in TRIM: the header displacements of
  low-energy runs (Ar -> Cu 1 keV: 14.9 -> 17.8 per ion, TRIM 18) and
  COLLISON.txt's VAC = DISP - REPLAC.
- COLLISON.txt with TRIM's "Ion + Recoils" (Collisions = 2,
  `output.collision_recoils`): every displacement and replacement of each
  cascade, and the damage summary after each ion; tables `ion_damage` and
  `cascade_atoms`.
- `tools/trim_compare/batch.py` and `analyze.py`: AMCSET seed runs against
  all TRIM runs of the same input, with t-statistics and profile chi^2.

### TRIM features and outputs

- TRIM.DAT particle lists (TRIM's cascade codes 4-7): ions with their own
  energy, start and direction, and recoil atoms without an ion (damage from
  neutrons, electrons, photons), from TRIM.IN + TRIM.DAT, the scenario
  (`source`), `--particles` or the GUI. BACKSCAT/TRANSMIT/SPUTTER/TRIMOUT
  files can be used as particle lists.
- EXYZ.txt (`output.exyz_interval`, `--exyz`), TRIMOUT.txt, and the
  surface-arrival energy histogram behind TRIM's sputtering plots.
- Results as tables with units (`common/include/output`, `common/src/output`): CSV, JSON (results.json)
  and Apache Parquet writers (own Parquet writer, no Arrow dependency).
  Output settings choose formats, tables, columns and their order, row
  ranges (with units), zero rows, decimation, units per dimension or column
  (Boost.Units conversions) and number formats (general, scientific,
  engineering, fixed, shortest); `--formats`, `--output-units`,
  `--number-format`, `--output-spec`, `output.files` (docs/OUTPUT.md).
- Calculations can be paused, and continued with more ions (also with
  another ion or energy on the same target); continuing gives the same
  results as one longer run.

### GUI rewritten

- New Qt 6 Widgets application in `gui/src/{app,model,run,setup,views}`
  (docs/GUI.md): undoable scenario editing with unit-aware fields, layer and
  element tables, compound formulas, periodic table, TRIM.IN import/export,
  results.json viewer.
- Particle paths drawn through QRhi (Vulkan/Direct3D/Metal/OpenGL) from an
  append-only GPU buffer: four views (depth-Y, depth-Z, beam, 3D) stay
  interactive with millions of segments. Replaces the QPainter canvases,
  whose redraw time grew with the number of paths, and the Qt 3D view.
- TRIM's plots (distributions, damage, energy loss, sputtering, spectra,
  EXYZ) live during the run, custom plots of any table, tables browser,
  TRIM's calculation-window statistics.
- Removed: Qt Graphs, Qt 3D and gflags dependencies of the GUI, the stub
  menu actions (print, cut/copy/paste), the Fontconfig workaround and the
  unused layer window.

### Build

- CMake options `AMCSET_BUILD_GUI`, `AMCSET_BUILD_CLI`, `AMCSET_BUILD_TESTS`;
  the core, CLI and tests build without Qt.
- `windeployqt` failures (e.g. blocked by an application-control policy) no
  longer fail the GUI build.
- Version 0.1.0 in CMake, Conan and `amcset version`.

### Documentation and tools

- Agent-oriented documentation for handing the project to coding agents:
  `AGENTS.md` (rules, repository map, build, conventions, definition of
  done), `CLAUDE.md` (workstation specifics) and a `CLAUDE.md` per module
  (`common/`, `common/src/trim/`, `gui/`, `cli/`, `test/`, `tools/`);
  new `docs/ARCHITECTURE.md`, `docs/HOWTO.md`, `docs/DEVELOPMENT.md`,
  `docs/GLOSSARY.md` and the index `docs/README.md`; `docs/KNOWN_ISSUES.md`,
  `docs/ROADMAP.md`, `docs/TESTING.md` and `tools/srim/README.md` brought up
  to date.
- `tools/trim_compare/doc_tables.py` prints the validation tables of
  `docs/VALIDATION.md` §2.2 from `analyze.py`'s `analysis.json`.
- `tools/srim/gen_cases3.py`, `run_trim_cases3.*` and `run_held_cases.*`:
  the third batch of TRIM reference cases.
- License header added to `common/src/data_files.cpp` and
  `common/include/data_files.h` (the CI header check failed on them).
