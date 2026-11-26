# Crystalline targets — design

Status (2026-11-24): **phases 1 and 2 of §8 are implemented**: the crystal
library, the scenario's `layers.N.crystal`, the `amcset cell` command
(phase 1, §10), and the crystal transport for single crystals with a static
lattice and thermal vibrations, with the beam spot (phase 2, §11). Of the
later phases two parts are in (§12): the GUI's Crystal tab for single
crystals (phase 4 without the grain preview) and a first comparison with
measurements (phase 6: a channeling dip and boron profiles in silicon, which
bracket the stopping division without fixing it). The user documentation is
docs/CLI.md "Crystals" and docs/GUI.md. Polycrystals (phase 3) are built in
the library, the scenario and the transport, without a GUI (§13, which also
lists what is left). **Still a proposal**: the dose-dependent lattice
(phase 5), a fitted stopping division and the wider validation (phase 6).
§10-§13 list where the built parts differ from the proposal below. As further parts are built, their user-facing
description moves to the owning document and this file keeps the design
rationale.

## 1. What is asked

- Targets aware of crystallographic orientation, **off by default**.
- Crystal information: which elements sit on lattice sites (substitutional,
  with occupancies) and which on interstitial sites; orientation; single
  crystal or polycrystal.
- Polycrystals from generated Voronoi grains, grain files, Neper
  tessellations and EBSD maps; ions entering over a surface area, not at
  one point.
- A unit-cell setup tab: common structures built in; CIF, VASP POSCAR,
  LAMMPS data and extended XYZ files; a 3D view of the cell where atoms are
  clicked and replaced; the beam direction set from the view ("the beam
  looks the way I look") or typed in.
- Physics: a static lattice as in MARLOWE by default; optionally a
  dose-dependent lattice (damage build-up dechannels later ions, as in
  Crystal-TRIM).
- Tests of crystal behaviour, channeling in particular, on Si implants, W
  and Fe, and Cu/Au (RBS-style angular scans).

## 2. Where it fits

**The crystal is a property of the target, not of one engine.** A layer
gets an optional `crystal` description in `config::RunConfig` (rule 6: one
scenario model). Amorphous layers stay exactly as they are.

**The transport is new, and it is not TRIM.** `trim2013` must keep
reproducing TRIM.exe (AGENTS.md rule 2), so crystals get a third transport,
`physics.transport: "crystal"`, next to `trim2013` and `refined`. It reuses
everything around the transport that the SRIM-mode framework already has:
`Mover`, per-ion random streams, blocks of 16 merged in order (thread-count
independence), `TrimTallies` and every output table, the dpa table, the
trajectory and MD-event files, the GUI's live plots and path views. In
amorphous layers it behaves like `refined` (random collision partners).

**Why not the native engine first.** The native engine is unvalidated
(docs/KNOWN_ISSUES.md §1) and has none of the output machinery; building
crystals there would mean building it twice. The lattice geometry, the
neighbour search, the file readers and the grain maps go into a library of
their own (`common/include/crystal/`, `common/src/crystal/`) with no
dependency on either engine, so the native engine (or another transport)
can use them later.

"SRIM mode" is then a misnomer for the framework. The proposal keeps the
name in the code and the docs say plainly that the crystal transport is
AMCSET's own physics; renaming the user-facing mode (e.g. "BCA") is a
separate, optional change.

## 3. Crystal description (`common/include/crystal/`)

### 3.1 Unit cell

- Lattice vectors **a**, **b**, **c** (Å), entered as a, b, c, α, β, γ or as
  vectors. A built-in structure with no lattice constant takes it from the
  layer density (a³ = n_cell · M̄ / ρ), so density and lattice agree; given
  a lattice constant, the density is derived from it instead and the layer
  shows the derived value.
- Basis sites, fractional coordinates. Each site has a **kind** (`lattice`
  or `interstitial`) and an **occupancy**: a list of (element, fraction),
  with vacancy = the rest. Substitutional alloys are lattice sites with
  several elements (random occupation, drawn once per site); interstitial
  solutes are interstitial sites with a small fraction (e.g. C on the
  octahedral sites of bcc Fe).
- Per element: Ed, Eb, Es (as now), and the thermal vibration amplitude:
  a Debye temperature (with the target temperature) or a given rms
  displacement.
- Built-in structures: sc, bcc, fcc, hcp, diamond, zinc blende, wurtzite,
  rock salt, CsCl, fluorite, perovskite, corundum, rutile, and the
  octahedral/tetrahedral interstitial sites of fcc, bcc and hcp. The
  element-specific lattice constants and Debye temperatures need a cited
  table before they become defaults (§9).

### 3.2 Files

| Format | Read | Notes |
|---|---|---|
| CIF | cell, space-group symmetry operations (`_symmetry_equiv_pos_as_xyz` / `_space_group_symop_operation_xyz`), atom sites with occupancies | the symmetry operations expand the asymmetric unit; duplicate sites merged within a tolerance; Hall/Hermann-Mauguin symbols alone (no operations) not supported at first |
| VASP POSCAR/CONTCAR | lattice and scale, species line (VASP 5), Direct/Cartesian | VASP 4 files (no species line) need the elements given |
| LAMMPS data | triclinic box (xy xz yz), Masses, Atoms (atomic/charge/full styles) | atom types mapped to elements by mass or by a given map; used as the unit cell (it may be a supercell) |
| Extended XYZ | `Lattice="..."`, `Properties=species:S:1:pos:R:3` | ASE/OVITO convention |

Writers for the same formats (except CIF symmetry) make round trips
testable and let users take an edited cell back to DFT/MD.

### 3.3 Orientation

Crystal → laboratory rotation; the laboratory frame stays TRIM's (x =
depth). Entered as any of:

- surface normal [hkl] and an in-plane direction [uvw] (Miller-Bravais
  four-index for hexagonal cells);
- Bunge Euler angles (φ1, Φ, φ2);
- beam direction in crystal coordinates, or the usual implantation
  **tilt** and **twist** about a named axis (e.g. Si(100), 7° tilt,
  22° twist);
- the GUI's view direction (§6).

### 3.4 Polycrystals

A grain map answers "which grain, with which orientation, is at (x, y, z)".
The lateral directions are periodic with a period the user sets (a
representative patch of the surface).

| Source | What it gives |
|---|---|
| Generated Voronoi | seeds from the run seed; mean grain size; orientations random (uniform on SO(3)), a fiber texture (axis ∥ normal, spread) or one orientation with a mosaic spread; optional columnar grains |
| Grain file (documented JSON/CSV) | seed points + Euler angles per grain, or a voxel grid of grain ids + one orientation per id |
| Neper `.tess` | cells, seeds and weights (Laguerre), orientations |
| EBSD `.ang` (TSL) / `.ctf` (Oxford) | a surface map of orientations (and phases); grains are columns into the depth |

Grain boundaries: an atom belongs to the grain its site lies in; an
optional boundary width is treated as amorphous (random partners), the
simplest model with a physical meaning (disordered boundary region).

### 3.5 Beam over an area

`ion.spot`: `point` (now), `rectangle`, `circle`, `gaussian`, with the
size. Each ion starts at a random point of the spot (from its own
stream). For a single crystal a spot of at least one surface unit cell is
needed anyway: ions starting at one point would all see the same rows.

## 4. Crystal transport (static lattice, MARLOWE-like)

Following the published algorithm of MARLOWE (M. T. Robinson and
I. M. Torrens, Phys. Rev. B 9, 5008 (1974); Robinson's later papers) and
Crystal-TRIM (M. Posselt). Their source codes are not public (MARLOWE 15b
is distributed through RSICC/NEA as PSR-137, Crystal-TRIM by HZDR on
request), and AMCSET is GPL: the implementation follows the papers and
manuals, never code.

Per flight step of a mover:

1. **Where am I.** Layer; grain (polycrystal); position and direction in
   that grain's crystal frame.
2. **Collision partners.** Lattice sites (and occupied interstitial sites)
   ahead of the mover within a maximum impact parameter p_max, found by
   walking the unit cells along the path (cell lists precomputed per unit
   cell), ordered by their projection on the path. The atom struck last is
   excluded (no immediate re-collision); sites emptied by this ion's
   cascade are skipped.
3. **Thermal vibrations.** Each site met gets an uncorrelated Gaussian
   displacement (Debye model rms amplitude), drawn when first met and kept
   for the rest of that ion's history.
4. **Simultaneous collisions.** Partners whose path coordinates lie within
   a tolerance are struck together: momentum transfers add, as in MARLOWE.
5. **Scattering.** The shared `physics::` scattering (magic formula or the
   scattering integral, ZBL or other screening), plus the asymptotic
   trajectory corrections (time integral, path displacement) that matter
   for channeled trajectories.
6. **Electronic stopping.** A nonlocal part (the existing stopping models
   and SRIM tables, over the path) and a local part depending on the
   distance of closest approach (Oen–Robinson form), mixed by a fraction;
   channeled ions lose less because they avoid close collisions. The
   mixing fraction and any channel corrections are the parameters that
   decide channeling depths; they are calibrated against measured profiles
   (§7), never on one case.
7. **Recoils.** A struck atom with T > Ed leaves its site (vacancy), is
   followed through the same lattice and stops as an interstitial.
   Optional close-pair recombination: a stopped interstitial within a
   capture radius of a vacancy fills it (MARLOWE-style pairing), giving a
   count of surviving Frenkel pairs next to TRIM-style vacancies and the
   dpa models.
8. **Amorphous regions** (amorphous layers, grain-boundary slabs, damaged
   regions of §5): random partners as in `refined`.

Determinism, threads, cancellation, pause, continue, trajectories and MD
events behave as in the other transports.

Performance: crystal transports cost more per step than random partners
(the neighbour search); the cell lists and the per-ion displacement cache
are the main measures. Expect several times `refined`'s cost; measured
before any claim.

## 5. Dose-dependent lattice (option)

Crystal-TRIM's approach, in outline: the fluence is split into steps; after
each step a depth profile of damage (from the transport's surviving defects
or a damage model) is updated, with saturation and a recombination factor;
in the next step a collision at depth x finds a disordered (random) partner
with a probability given by the local damage fraction. Lateral homogeneity
is assumed (depth profile only). Inputs: fluence, number of steps, the
damage-to-disorder parameters. Outputs: ranges and damage per step.
Deterministic (steps run in order; blocks within a step as now).

## 6. GUI

- **Crystal tab** (setup dock, next to Target): per layer amorphous /
  single crystal / polycrystal; structure list; *Import…* (CIF, POSCAR,
  LAMMPS data, extended XYZ) and *Export…*; lattice parameters; a table of
  sites (kind, x, y, z, occupancy); Debye temperatures and temperature;
  orientation fields (normal [hkl] + in-plane [uvw], Euler angles,
  tilt/twist).
- **3D unit-cell view** (QRhi, like the path views): atoms as spheres
  coloured by element, cell edges, n×n×n cells, orbit/zoom. Click an atom →
  menu: change element (periodic table), set occupancy, make interstitial,
  delete; edits are undoable (`Scenario::modify`).
- **Beam direction from the view**: *Beam along view* sets the beam
  direction to the current view direction (crystal coordinates shown as
  [uvw] and as tilt/twist); typing a direction turns the view to it. The
  surface plane and the beam arrow are drawn.
- **Polycrystal**: grain generation settings and file import; a preview
  slice coloured by orientation.
- Plots: as now; an angular-scan result plot (below).

CLI: the same keys in the scenario; `amcset cell` (print/convert a unit
cell file); `amcset scan` (runs over tilt angles, writing a
yield-versus-angle table for channeling studies).

## 7. Validation plan

Unit tests (exact, no statistics):

- built-in cells: nearest-neighbour distances and coordination numbers;
  densities from lattice constants and back;
- readers: small embedded CIF (with symmetry), POSCAR, LAMMPS data, XYZ →
  the same cell; writers round trip;
- orientation: Miller → matrix → Euler → matrix; tilt/twist;
- neighbour search against brute force over random paths;
- thermal displacements: rms against the Debye amplitude;
- Voronoi/voxel/Neper/EBSD grain lookup against brute force;
- results independent of the thread count.

Physics tests (statistical, fixed seeds):

- **Random limit**: a crystal in a random (non-channeling) direction with
  vibrations gives ranges close to the amorphous transport.
- **Channeling**: Si ⟨100⟩ at 0° gives deeper B ranges with a channeling
  tail than at 7° tilt; tails shrink with temperature.
- **Critical angle**: He in Cu, Au, W along low-index axes: the
  close-encounter yield against tilt has a half-width of the order of
  Lindhard's ψ₁ = (2 Z₁ Z₂ e² / (E d))^½ (d: atomic spacing along the row).
- **Dose dependence**: tails shrink with fluence.

Literature comparisons (calibration of the electronic-stopping parameters
on part of the data, checked on the rest): B, P, As in Si (SIMS profiles,
channeled and random); He and self-ions in W and Fe (channeling tails, MD
range calculations); RBS channeling minimum yields and critical angles in
Cu and Au.

## 8. Phases

| Phase | Content | Ends with |
|---|---|---|
| 1 | crystal library: unit cell, built-ins, orientation, readers/writers, scenario keys, CLI `cell` | unit tests; scenario round trip |
| 2 | crystal transport, single crystal, static lattice, beam spot | channeling and random-limit tests; first Si/W comparisons |
| 3 | polycrystals: Voronoi, grain files, Neper, EBSD; grain-boundary width | grain lookup tests; polycrystal runs |
| 4 | GUI: crystal tab, 3D unit-cell editor, beam-from-view, grain preview | smoke test screenshots |
| 5 | dose-dependent lattice; `amcset scan` | dose tests |
| 6 | calibration and validation against the literature; docs | docs/CRYSTAL_MODE.md with the comparisons |

Each phase builds, passes all tests, leaves `trim2013` and `refined`
bit-identical (`tools/perf/bench.py compare`) and is documented.

## 9. Open points

(See also §10 for what phase 1 settled.)

- **Literature access.** The algorithms and validation data are in papers
  behind paywalls (Robinson & Torrens 1974; Posselt's Crystal-TRIM papers;
  SIMS channeling profiles; RBS channeling data). The user offered TAMU
  access; the papers would be read in the browser pane after the user signs
  in.
- **MARLOWE / Crystal-TRIM themselves.** Not downloadable. If the user
  obtains MARLOWE 15b through RSICC (licence) or Crystal-TRIM from HZDR,
  their manuals and outputs can be used for comparison runs; their source
  is never copied into AMCSET.
- **Default material data.** Lattice constants and Debye temperatures per
  element need a cited table before they become defaults; until then they
  are inputs (lattice constants can come from the density).
- **Native engine.** Uses the crystal library later, if wanted.

## 10. Phase 1 as built (2026-11-24)

Code: `common/include/crystal/` and `common/src/crystal/` (`geometry.h`,
`unit_cell`, `orientation`, `crystal`, `cell_io`; no dependency on an
engine), `common/src/trim/layer_crystal.cpp` (`apply_crystal`,
`check_layer_crystal`), the crystal part of `run_config.cpp`, `print_cell`
in `cli/src/cli.cpp`. Tests: `test/common/crystal_test.cpp` (21 tests) and
the `RunConfig.Crystal*` and `Cli.CellCommand` tests.

Differences from the proposal:

- **One stored form.** A crystal is kept and saved as its cell written out
  (lattice vectors and sites). `structure` and `file` are ways to write a
  scenario, not stored: a saved scenario never depends on a cell file, and
  the GUI will edit one representation.
- **The crystal owns composition and density** (§3.1 said "agree"): the
  layer's stoichiometry and density are set from the cell
  (`trim::apply_crystal`) and `TrimConfig::validate` checks that they still
  are. A density in the scenario within 2 % of the cell's is accepted and
  replaced; beyond that it is an error.
- **Built-in structures**: sc, bcc, fcc, hcp, diamond, zincblende, wurtzite,
  rocksalt, cesium_chloride, fluorite, perovskite, L1_2, L1_0. Corundum and
  rutile (proposed) are left to cell files: they need internal parameters
  that would have to be cited.
- **Material data.** Phase 1 had none. Since §12 the cubic crystals of
  Gemmell's Tables IV and V are built in (lattice constants; Debye
  temperatures of the elements), with the citation; for anything else the
  lattice constant is an input or comes from the density, and Debye
  temperatures are inputs.
- **CIF** needs the symmetry operations in the file (no space-group
  tables); **cell writers** exist for all four formats.
- **The beam direction** is not a scenario key: `ion.angle` and
  `orientation.twist` stay explicit, and `amcset cell --beam` computes them
  (the GUI's "beam along view" will do the same).
- **The beam spot** (§3.5) moved to phase 2, with the transport that needs
  it: it changes where ions start.
- **Thermal amplitudes**: the Debye-model formula is tested against the
  oscillator average over the Debye spectrum (its definition) and its
  limits; the amplitudes tabulated in the channeling literature (Gemmell,
  Rev. Mod. Phys. 46, 129 (1974)) have not been compared yet (the journal
  page could not be read in this session).
- Every existing output is unchanged: `tools/perf/bench.py compare`,
  390 of 390 files identical.

## 11. Phase 2 as built (2026-11-24): the crystal transport

Source: M. T. Robinson and I. M. Torrens, "Computer simulation of
atomic-displacement cascades in solids in the binary-collision
approximation", Phys. Rev. B 9, 5008 (1974), read in this session (the
paper, not the MARLOWE program). Code: `crystal::LatticeSearch`
(`common/*/crystal/lattice_search.*`, engine-independent) and
`Transporter::run_ion_crystal`, `crystal_follow`, `crystal_strike`,
`crystal_path_loss` in `common/src/trim/trim_engine.cpp`. Checks:
`tools/crystal/channeling_checks.py`; results in docs/CLI.md "The crystal
transport".

What follows the paper:

- **Straight segments between deflection points** (Sec. II). The deflection
  point lies x₁ before the closest approach on the incoming line and the
  struck atom starts x₂ from its position (Eq. 11: x₁ = [(1 + f)τ +
  (fA − 1) s tan(θ/2)] / [f(1 + A)], x₂ = s tan(θ/2) − x₁).
- **The partner search** (Sec. IV, Eqs. 20-21): a site qualifies when its
  path coordinate ξ is positive and its impact parameter s is below s_max.
  Atoms already struck are not struck again while they are still ahead (the
  paper's ξ_min; here by the sites' identities).
- **Simultaneous collisions** (Sec. IV): the atoms within a distance of the
  first along the path are each treated as if alone and the deflections
  added vectorially; if the transfers exceed the projectile's energy they
  are all scaled so that it just vanishes.
- **Thermal vibrations** (Sec. IV): static, uncorrelated Gaussian
  displacements along each axis, from the Debye model.
- **Displacement** (Sec. V): displaced if T > Ed, leaving with T − Eb; the
  projectile replaces the target when it keeps less than the lesser of its
  cut-off energy and Ed. With the cut-off at Ed this is Kinchin and Pease's
  model; AMCSET's default cut-off is TRIM's 1 eV.

What differs from the paper, and why:

| Paper (MARLOWE) | AMCSET | Reason |
|---|---|---|
| Molière potential | ZBL universal potential, the magic formula (or the scattering integral) | the shared `physics::` code; one potential for crystal and amorphous layers |
| Time integral τ by four-point Gauss-Mehler quadrature | its hard-core value τ = s tan(θ/2): x₁ = s tan(θ/2), x₂ = 0 (struck atoms start from their own positions) | what TRIM-type codes use; exact for glancing collisions, which decide channeling. A quadrature for τ is a possible refinement |
| Inelastic loss by a modified Firsov model, local only | AMCSET's stopping power (SRIM tables or models), divided between the path and the collisions; the collision part falls off as exp(−c R/a) with the apsis R, normalised so that straight paths at all impact parameters lose S_e (the form of Oen and Robinson's local loss with the stopping power in place of their k√E; the 0.045/πa² of that form is c²/2πa² for c = 0.3) | TRIM-quality stopping in random directions; a local part for channeling. The fraction and c are defaults, not fits |
| Neighbour lists around a lattice site, orthorhombic cells | all sites of the cells within reach, any cell | general unit cells and cell files |
| Thermal displacements drawn as the cascade develops, by table look-up | from a hash of the site and the ion: one frozen lattice per ion, the same from wherever a site is seen; cut at 3.5 σ per axis | no storage, no dependence on the order of visits or the thread count |
| Follows the fastest atom of the cascade first | a stack (last displaced first) | the order matters only through the sites vacated meanwhile |
| Separation-dependent Frenkel-pair recombination after the cascade | not done | phase 5/6, with the dpa models |
| Channeled trajectories may be truncated | followed to the end | simpler; costs time only |

Other choices: s_max defaults to the radius of the largest hole of the
lattice (the paper asks that no part of the crystal be free of collisions);
each layer has its own lattice, anchored at its top; the ion's entry point
(the beam spot) shifts the lattice instead of the ion, so that results stay
relative to the entry point; the ion starts in front of the surface so that
every surface atom lies ahead of it; amorphous layers take one random
partner per monolayer of path (counted from the last partner's closest
approach, not from the deflection point x₁ before it); the electronic loss
along the path is charged once (the stretch from a deflection point to the
closest approach is charged with the flight before it, not again with the
next), and its straggling is a gamma deviate (exact mean; a Gaussian cut at
zero raises the mean loss of the short steps).

What the checks showed (numbers in docs/CLI.md):

- A random direction reproduces the amorphous target (range within 2 %,
  quick-damage vacancies within 3 %), independent of how the stopping is
  divided; a fast proton on a straight path loses the stopping power to
  0.2 % whatever the division and s_max.
- Channels are deeper and cleaner, ⟨110⟩ more than ⟨001⟩ in Si; a 100 Å
  oxide, a higher temperature and a larger thermal amplitude each reduce
  the channeled range; a static lattice channels best.
- The channeling dip of 1 MeV He in W ⟨001⟩ has a minimum yield of 2 % and
  a half-width of 0.82 ψ₁. (A measured dip, 1 MeV He in Ge ⟨110⟩, is
  compared in §12.)
- With the cut-off at Ed, full cascades give 15 % more vacancies than NRT
  (293.9 against 256.0, 30 keV B in Si); with TRIM's 1 eV cut-off 436.7,
  which is 15 % above TRIM's own full cascades in amorphous Si (379.8). The
  second difference is not understood in detail (it does not come from the
  simultaneous collisions: sequential collisions give the same count); the
  quick-damage counts agree.
- Slow heavy ions channel in many directions: 20 keV Fe in bcc Fe at 7° tilt
  goes four times deeper than in amorphous Fe, and the average over random
  orientations (188 Å) is twice the amorphous range (89 Å), while the
  least-channeling orientation matches it (83 Å).

Open for the next phases: a fit of the stopping division to low-dose
profiles with stated tilt and twist (§12 only brackets it); P and As in Si,
channeling tails in W and Fe, RBS dips in Cu and Au; the time integral;
Frenkel-pair recombination; a search that prunes by direction (the
transport is 10-30× slower than `trim2013`); polycrystals and the
dose-dependent lattice (§8).

## 12. Common crystals, measurements and the GUI tab, as built

**Common crystals** (`crystal/materials.h`). Tables IV and V of Gemmell,
Rev. Mod. Phys. 46, 129 (1974): 19 cubic elements (structure, lattice
constant at room temperature, Debye temperature, u₁ at 293 K) and 24 cubic
compounds (structure and lattice constant). `thermal_amplitude()` is what
the transport uses for an element: the crystal's own entry; else, in a
crystal of one element that the table lists, the Debye model with the
tabulated Debye temperature (the run says so); else none. Compounds get no
default: the table has no Debye temperatures for them, and an element's
value does not carry over to a compound. `lindhard_angle()` is Gemmell's
Eq. (A5).

**Measurements** (`tools/crystal/calibration.py`, results and discussion in
docs/CLI.md "Against measurements"). The Debye amplitudes reproduce Table
IV. The dip of 1 MeV He along ⟨110⟩ of Ge comes out with a half-angle of
0.99° (measured 0.95°) and a minimum yield of 0.022 (measured 0.03); the
yield is the number of atoms the ion displaces between 200 and 1200 Å, a
close-collision yield like backscattering. Boron profiles in (100) Si at
0.5, 15 and 30 keV are reproduced in peak and width for every stopping
division, and bracket the division (0.25 to 1) without fixing it: tilt and
twist are mostly not stated, the higher doses carry damage, and the deep
SIMS tail at 0.5 keV is beyond any division. The default 0.5 is kept, and
the transport's warning says what is and is not established.

**GUI** (`gui/src/setup/crystal_panel.*`, `gui/src/views/cell_view.*`;
docs/GUI.md). Differences from §6:

- The unit-cell view is painted with QPainter, not QRhi: a cell has tens of
  atoms, and a parallel projection drawn back to front needs no depth
  buffer. More than 6000 atoms (a supercell from a file) are not drawn.
- Clicking an atom opens a menu (replace the element, share the site with
  another element, set the occupancy, lattice or interstitial, remove), for
  that site or for all sites like it. Sites are also edited in a table.
  There is no dragging of atoms and no adding an atom by clicking in space
  (*Add site...* takes fractional coordinates).
- *Point the beam along this view* sets `ion.angle` and the twist with
  `tilt_twist_for_beam`; a typed beam direction does the same and turns the
  view. The view does not follow the beam by itself.
- Presets are Gemmell's 43 crystals; structures are built from a dialog
  (the lattice constant, or the density); the four file formats import and
  export.
- No polycrystal settings and no grain preview yet (§13).

## 13. Phase 3 (polycrystals), as built so far (2026-11-24)

**On hold** (the maintainer's decision, 2026-11-25): work concentrates on
single crystals and amorphous targets. What is built stays, tested and
documented; the list of what is left (below) waits.

**Built** (library, scenario and transport; tested):

- `crystal/grain_map.h`: `GrainMap` (which grain is at a point of a layer,
  and how far the nearest grain boundary is): seeds of a Voronoi or
  Laguerre tessellation found through buckets (periodic laterally; several
  seeds may be one grain), columns, or a grid of voxels. `Polycrystal` (the
  description in `Crystal::polycrystal`), `generate_voronoi` (random,
  fiber and mosaic textures; its own seed), `lattice_rotations` and
  `misorientation_rad` (cubic, hexagonal, tetragonal, orthorhombic).
- `crystal/grain_io.h`: AMCSET's grain file (seeds with Bunge angles, or
  voxels), Neper `.tess` (seeds, weights, `*ori`), EBSD `.ang` and `.ctf`
  (points grouped into grains by misorientation; unindexed points filled).
- Scenario: `layers.N.crystal.polycrystal` (docs/CLI.md "Polycrystals");
  file sources are read when the scenario is parsed.
- Transport (`crystal_follow`): one `LatticeSearch` per grain; the mover's
  grain decides the lattice, its origin and its configuration; near a
  boundary only the sites of the mover's grain exist; a boundary slab of
  `boundary_width` takes random partners; sites are told apart by lattice
  (`SiteKey`), which also separates the sites of two crystal layers.
- Tests `CrystalGrains` (7), `TrimCrystal.PolycrystalGrains`,
  `RunConfig.PolycrystalKeys`; `tools/crystal/polycrystal_checks.py`.

**Differences from §3.4.** The grain file is a keyword text file, not
JSON/CSV. The generator's seed is the polycrystal's own, not the run's. EBSD
grains take the orientation of their first point (no averaging), and
unindexed points take the nearest grain instead of becoming amorphous.

**Not built yet.**

- GUI: polycrystal settings in the Crystal tab and a preview slice coloured
  by orientation (the rest of phase 4). The tab shows a polycrystal layer's
  unit cell and keeps its `polycrystal` when the cell is replaced, but
  cannot create or edit one.
- A warning when `ion.spot` is smaller than the grains' period.
- A check of the Neper and EBSD readers against real files (orientation
  conventions, axes); TSL hexagonal grids are read as points but only a
  square grid is tested.
- A summary of the grains for the user (`amcset cell` or a `grains`
  command: number of grains, size distribution, texture), and the grains in
  results.json.
- Grains in registry across low-angle boundaries (a mosaic crystal
  dechannels too much: docs/CLI.md "Polycrystals", Limits).
- Then phase 5 (dose-dependent lattice, `amcset scan`) and the rest of
  phase 6 (a fit of the stopping division; W, Fe, Cu/Au).
