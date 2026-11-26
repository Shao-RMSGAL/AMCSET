# The AMCSET desktop application

`amcset-gui [scenario.json | TRIM.IN | results.json]`

The window has the setup on the left, the results in the middle and the
running statistics on the right (both side panels can be moved, floated or
hidden from the View menu).

## Setup

- **Ion**: element (or the periodic table, *PT*), mass (default: most
  abundant isotope, as TRIM), energy and angle; values accept units
  ("2 MeV", "250 keV"). The particle source is either the ion beam or a list
  in TRIM.DAT format: ions with their own energy, start and direction
  (TRIM's cascade codes 4 and 5) or recoil atoms set in motion without an
  ion, for damage by neutrons, electrons or photons: with full cascades
  (code 6) each atom starts a cascade, with monolayer steps (code 7) it
  moves like an ion.
  BACKSCAT.txt, TRANSMIT.txt, SPUTTER.txt and TRIMOUT.txt can be loaded as
  lists too.
- **Target**: layers from the surface (name, thickness, density, gas phase,
  Bragg compound correction) and the elements of the selected layer
  (stoichiometry, mass, Ed, Eb, Es). Cells accept units ("1.5 um",
  "7.87 g/cm3"). *Compound...* adds a layer from a chemical formula
  ("Fe2O3", "Ca5(PO4)3OH", "Fe0.91Cr0.09") or a list of common compounds
  with their densities. Ticking *Gas* switches an element's default density
  between its solid and gas values, as TRIM does.
- **Crystal** (not TRIM): the crystal of a layer. Choose the layer and
  *Single crystal*; a pure element of the table of common crystals becomes
  that crystal, anything else is asked for. The unit cell comes from
  *Common crystals...* (43 cubic crystals with their lattice constants and,
  for the elements, Debye temperatures: Gemmell, Rev. Mod. Phys. 46, 129
  (1974)), from *Structure...* (a built-in structure with the elements of
  its sublattices, substitutional mixtures such as "Fe 0.9, Cr 0.1",
  interstitial atoms, and the lattice constant or the density) or from
  *Import...* (CIF, VASP POSCAR, LAMMPS data, extended XYZ); *Export...*
  writes it. The lattice parameters and the table of sites (fractional
  coordinates, lattice or interstitial, occupants) can be edited. The
  layer's elements, stoichiometry and density are then the unit cell's (the
  Target tab shows them read-only), and a change that gives the layer other
  elements also renames it after them ("GaAs", "UO2"; the name can be
  changed in the Target tab).
  The cell is drawn in three dimensions: drag to turn, Shift+drag to roll,
  wheel to zoom; shared sites are pie charts, interstitial sites smaller,
  empty sites rings. **Click an atom** to replace its element (periodic
  table), share the site with another element, set its occupancy, make it a
  lattice or an interstitial site, or remove it, for that site or for all
  sites like it. The projection is parallel, so looking along a crystal
  direction shows its channels: type a direction in *View along*, or use
  *Beam*, *Normal* and *Snap*. **Point the beam along this view** sets the
  beam's angle (Ion tab) and the crystal's twist so that the ions travel
  the way you look; typing a direction in *Beam along* does the same and
  turns the view. **Turn the crystal to this view** sets the crystal's
  orientation instead: the direction you look along becomes the surface
  normal (you look into the target) and the top of the view the lateral Y
  axis, towards which the beam tilts; the beam's angle stays. A view along
  the normal of a low-index plane (indices up to 6, within 0.05°: type a
  direction or use *Snap*) gives that plane as the surface in Miller
  indices, with a low-index reference direction or else the twist; any
  other view gives Bunge Euler angles, which are exact. The corner of the
  view shows the beam and the depth axis.
  Below: the orientation (surface plane and a direction in the surface, or
  Bunge Euler angles; twist), the temperature and each element's thermal
  vibrations (Debye temperature or rms amplitude; *In use* shows what the
  transport takes, including the tabulated value of a pure element), and
  the settings of the crystal transport (`physics.crystal`, `ion.spot`).
  Polycrystals (`layers.N.crystal.polycrystal`, docs/CLI.md "Polycrystals")
  cannot be set up here yet: a scenario that has one opens and runs, the
  tab lists the layer as a polycrystal and shows its unit cell, and
  replacing the cell keeps the grains' description.
  A crystal is followed by the crystal transport only: the box *Crystal
  transport (follows the lattice)* at the top of the tab switches it on and
  off (the same switch as in the Calculation tab), and the tab says so
  while it is off.
- **Calculation**: engine (SRIM mode or the native engine), TRIM's damage
  calculation (quick Kinchin-Pease, full cascades, monolayer steps), number
  of ions, seed (or *New random seed for every run*: the drawn seed is
  shown in the log and recorded in the outputs; continuing keeps it),
  threads; for SRIM mode the transport of amorphous targets
  (TRIM-2013 or refined) and the switch *Crystal transport (follows the
  lattice)*, which turns the crystal transport on and off: on, layers that
  have a crystal are followed atom by atom (not TRIM, and the division of
  its electronic stopping is not calibrated, docs/CLI.md "Crystals"), the
  list above it waits, and the
  *TRIM.exe bugs* below are shown all corrected, because that transport
  takes none of them; off, the transport chosen in the list is used again
  with the bug switches as they were. Then the electronic stopping (SRIM tables, built-in, Bethe-Bloch,
  Lindhard-Scharff), *Correct the stopping of Ni ions*
  (`physics.stopping_corrections.nickel_ions`, off by default, not TRIM:
  docs/VALIDATION.md §9), straggling, final energy, the largest cascade and
  *Fast approximations* (`physics.fast_math`, off by default: about 2x
  faster through a tabulated scattering angle; the tooltip warns that the
  results agree with the exact calculation statistically but not number
  for number, so the option may not faithfully replicate TRIM step by
  step; see docs/SRIM_MODE.md §8.2); for the native engine its own
  options. The group *TRIM.exe bugs* (SRIM mode,
  trim2013 transport only) has one check box per known TRIM.exe bug
  (checked = reproduce it, as TRIM does; each tooltip says what the bug does
  and that unchecking it departs from TRIM's results), a box for TRIM's
  single-precision table sums, and *Reproduce all* / *Fix all*. The list
  and the evidence are in docs/SRIM_MODE.md §7 and
  docs/TRIM2013_ALGORITHM.md. The group *Target depth* (SRIM mode) has
  *Automatic target depth* (`auto_depth.enabled`, off by default: before
  each run, quick test runs size the last layer so that the deepest test ion
  stops in its last part; the sized thickness then appears on the Target
  tab as an undoable edit), the number of *Test ions* (1000) and the
  *Margin* (10 %), and *Warn when ions leave through the back* (on by
  default). docs/CLI.md "Checks and warnings" describes the procedure,
  also for the crystal transport (a deep box first, then trimmed to the
  deepest channeled ion).
- **Output**: folder and formats (TRIM text, CSV, JSON, Parquet), the depth
  window and number of bins, the optional lists (ion end points, atoms
  leaving the target, collision details for the first N ions, the EXYZ
  energy step), the number of ion paths drawn (or *Record every ion's
  path*) and *dpa at fluence* (ions/cm², `output.dpa_fluence`, the fluence
  of the dpa plot and table; the arc-dpa/rpa constants are edited in the
  scenario file only). The group *Trajectory and MD event files (not TRIM outputs)*
  (SRIM mode, every transport) has a switch for each file: *Write a
  trajectory file* (then `trajectories.parquet` in the output folder unless
  another name is typed; every recorded point, or the births only, above a
  minimum energy; if no paths are recorded it records every ion's) and
  *Write an MD event file* (`events.parquet`; the list of transition
  energies E_c, E_min, record or hand-off mode; needs full cascades). The
  settings of a file are greyed out while its switch is off; each field's
  tooltip says what it does
  and docs/OUTPUT.md describes the files. The MD totals appear in the
  Tables tab (`md_events_summary`). *Tables, units,
  numbers...* opens the output settings (docs/OUTPUT.md) with a live
  preview of each table.

Every edit is checked as it is made: a value outside what SRIM mode can
simulate (an energy above 100 GeV, a layer thinner than one atomic layer,
a negative or non-finite number, a mass that is not an isotope of the
element, ...) is refused with a message that says what is accepted and
what to use instead, and the field shows the previous value again (the
limits are listed in docs/CLI.md "Checks and warnings"). Values that are
allowed but look unintended (a density far from the material's, a gas flag
that does not fit the density, a target far thicker than the ions travel,
a relativistic ion) are counted in the status bar (*n setup warnings*,
the list in its tooltip) and written once to the Messages tab.

Every change can be undone (Edit > Undo). Scenarios are saved as JSON
(File > Save) and TRIM inputs can be imported and exported (File > Import
TRIM.IN, Export TRIM.IN, which also writes TRIM.DAT for particle lists).

## Running

*Run* (F5) starts the calculation; *Pause* (F6) holds it; *Stop* ends it
early (the results cover the finished ions). *Continue* (F7) adds the given
number of ions to the finished calculation with the current ion settings:
more statistics, or, as in TRIM, a second implant of another ion, energy or
angle on the same target superimposed on the first. Because every ion has
its own random-number stream, 1000 ions followed by 1000 more give exactly
the results of one 2000-ion run.

With the automatic target depth on, *Run* first shows "Sizing the
target..." while the test runs go on (on the worker thread; *Stop* cancels
them), writes what was done to the Messages tab and then starts the
calculation on the sized target. When ions leave through the back of the
target during a SRIM-mode run (and the warning is on), the Messages tab
and a non-modal message say how many; the run goes on. The message offers
to stop and run again with the automatic target depth (or, if it was
already on, with ten times as many test ions), to continue, and to turn the
warning off for the scenario. The end of the run logs the number of
transmitted ions.

**Diagnostic log and crash reports.** The GUI keeps a log of run-level
events (start, warnings, errors, files written; Qt's own warnings) in
`%LOCALAPPDATA%\AMCSET\AMCSET\logs\amcset-gui.log` (Linux
`~/.local/share/AMCSET/AMCSET/logs/`; the previous log is kept as
`amcset-gui.old.log` beyond 5 MB). *Help > Detailed diagnostic log* adds
debug lines (it never logs single collisions, so it costs no speed). If the
program crashes, a crash report (stack trace, the calculation, the ion each
worker was transporting and how to rerun it, the last log lines) and on
Windows a minidump are written to the same folder; *Help > Open log folder*
opens it. Please attach both to a bug report.

## Results

- **Paths**: the ion and recoil paths in four views (depth-Y, depth-Z, the
  beam view and 3D) or one of them. Drag to pan (2D) or rotate (3D;
  right-drag pans), use the wheel to zoom, double-click to fit. Colours by
  element (the ion in red) or by recoil generation; the incident paths,
  recoil paths up to a chosen generation and the opacity (low values show
  where many paths overlap) can be set. The paths are drawn on the GPU
  (Vulkan, Direct3D, Metal or OpenGL) from a buffer that only grows, so a
  frame costs the same however many paths came before it, and millions of
  segments stay interactive. View > Path memory limit sets how many
  segments are kept (1 million take about 72 MB). If the graphics API does
  not work on a machine, set the environment variable `AMCSET_RHI` to
  `opengl`, `vulkan`, `d3d11` or `d3d12`.
- **Plots**: TRIM's plots, live while the calculation runs: ion
  distribution, ions and recoils, lateral spread, ionization, phonons,
  energy to recoils, energy absorbed by each element, vacancies, collision
  events, replacements, the integral and differential sputtering plots
  (with the surface binding energy marked), the energy spectrum of displaced
  atoms, the EXYZ energy-depth points and the final ion positions, dpa by
  damage model (TRIM vacancies, NRT, Kinchin-Pease, arc-dpa) and replacements
  per atom (rpa; SRIM mode, TRIM-2013 transport, docs/OUTPUT.md "dpa"), plus a
  custom plot of any columns of any table. Drag a rectangle to zoom, wheel to
  zoom the depth axis, double-click to reset; the context menu switches log
  scales, copies the data or saves the image. Lengths and energies can be
  shown in any unit.
- **Tables**: every result table, raw or as the output settings shape it,
  with copy, CSV and Parquet export.
- **Statistics**: TRIM's calculation-window numbers (ranges and straggling,
  backscattering, transmission, sputtering per element, damage per ion, the
  energy balance), with the progress, rate and remaining time.

File > Open Results reads a `results.json` written by an earlier run.

## Quitting, signals and where the GUI keeps its files

Closing the window asks to stop a running calculation and to save a changed
scenario. When the desktop session ends (log-off, shutdown, an installer's
restart), the GUI does not ask: it stops the calculation, keeps a changed
scenario as `recovery/scenario-<date>-<time>.json` in its data folder, and
the next start names it in the Messages tab. On Linux, SIGTERM, SIGHUP and
SIGINT (Ctrl+C in the terminal that started it) close it the same way, with
the exit status 128 + the signal number.

| | Windows | Linux |
|---|---|---|
| Settings (window layout, recent files, folders) | registry, `HKEY_CURRENT_USER\Software\AMCSET\AMCSET` | `~/.config/AMCSET/AMCSET.conf` |
| Log and crash reports (`logs/`), `recovery/` | `%LOCALAPPDATA%\AMCSET\AMCSET\` | `~/.local/share/AMCSET/AMCSET/` |

## Code layout

`gui/src/` is organised by role:

| Folder | Contents |
| --- | --- |
| `app/` | `MainWindow`: actions, menus, docks; connects the parts |
| `model/` | `Scenario` (the edited RunConfig with undo), table models of layers and elements, unit-aware text parsing |
| `run/` | `RunController` (engine on a worker thread, live snapshots, pause/stop/continue), `TrackStore` (path segments), datasets of the native engine and of results.json |
| `setup/` | the setup panels (ion, target, crystal, calculation, output), `FormBinder` (two-way binding of widgets to the scenario), the periodic table and the output settings dialog |
| `views/` | `TrackView` (QRhi path views) and `TrackPanel`, `PlotWidget` and `PlotsPanel`, `TablesPanel`, `StatsPanel`, `CellView` (the unit cell in three dimensions, QPainter) |

`gui/shaders/` holds the GLSL sources of the path shaders and their compiled
`.qsb` files in `baked/`. The build uses the baked files by default, because
some Windows machines block running `qsb.exe`; configure with
`-DAMCSET_COMPILE_SHADERS=ON` to compile them at build time instead. After
editing a shader, regenerate the baked files with `gui/shaders/bake.sh`
(set `QSB` to the path of Qt's `qsb` tool).

For automated checks the application accepts `--run` (start at once),
`--screenshots DIR` (save images of every view, `plot_<id>.png` of every
plot, `setup<n>.png` of every setup tab and `setup<n>_full.png` of the
whole of each tab that scrolls, and the Crystal tab and its cell view as
`crystal_tab.png` and `crystal_cell.png`, when the run ends) and
`--quit`.
