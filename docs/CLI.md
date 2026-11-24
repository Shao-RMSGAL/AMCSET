# The `amcset` command-line program

`amcset` runs everything the GUI can run, without the GUI: both engines
(SRIM mode and the native engine), TRIM batch inputs, stopping-power tables,
element data and scenario files. Scenario files are JSON and are the same
files the GUI opens and saves (File > Open / Save).

```
amcset run <scenario.json> [options]      Run a scenario (SRIM mode or native)
amcset trim <TRIM.IN> [options]           Run a TRIM batch input in SRIM mode
amcset convert <in> <out>                 Convert TRIM.IN <-> scenario JSON
amcset stopping [options]                 Stopping powers and ranges (like SR Module)
amcset elements [symbol|Z ...]            Element data (masses, densities, Es, I)
amcset cell <file> | --structure NAME     A crystal's unit cell: summary, conversion, scenario object
amcset template [srim|native] [-o FILE]   Print (or save) an example scenario
amcset validate <scenario.json|TRIM.IN>   Check a scenario without running it
amcset version | help
```

`amcset help` prints the full option list. Physical values always carry
units: `"1 MeV"`, `"1.2 um"`, `"7.87 g/cm3"`, `"40 eV"`, `"7 deg"`. Accepted
units: energy `eV keV MeV GeV`; length `A nm um mm cm m`; mass `amu u g kg`;
density `g/cm3 kg/m3`; angle `deg rad`; temperature `K C`.

## Quick start

```bash
amcset template > fe.json                      # Fe (1 MeV) into Fe, SRIM mode
amcset run fe.json -n 5000 -o out_fe           # TRIM-format files in out_fe/
amcset run fe.json --damage full --set "ion.energy=2 MeV" -o out_fe_fc
amcset trim "C:\SRIM-2013\TRIM.IN" -o out_trim # run an existing TRIM input
amcset convert fe.json TRIM.IN                 # ... or write one for TRIM
amcset stopping --ion He --target "Fe:0.9,Cr:0.1" --density "7.8 g/cm3" \
               --from "10 keV" --to "100 MeV" --units "MeV/(mg/cm2)"
```

SRIM mode uses the SRIM-2013 stopping tables that come with AMCSET
(`data/srim2013`, `docs/SRIM_MODE.md` §4) unless told otherwise:

```bash
amcset run fe.json --stopping builtin                       # no data files
amcset run fe.json --stopping-table my_tables.amcset        # other tables
amcset trim TRIM.IN -o out --fix-trim-bugs                  # correct TRIM's bugs
```

## Options of `run`, `trim` and `validate`

| Option | Meaning |
|---|---|
| `-o, --output DIR` | Output directory (default `output`) |
| `-n, --ions N` | Number of incident ions |
| `-j, --threads N` | Worker threads (0 or omitted: all cores). Results do not depend on it. |
| `-s, --seed S` | Random-number seed (SRIM mode; TRIM's default 716381) |
| `--mode srim\|native` | Engine (scenario files only) |
| `--damage kp\|full\|monolayer` | TRIM's damage calculation: quick Kinchin-Pease, full cascades, monolayer steps |
| `--stopping MODEL` | `srim`, `builtin`, `bethe-bloch` or `lss` |
| `--stopping-table PATH` | SRIM stopping-table file, or a folder of SR Module outputs (default: the SRIM-2013 tables in `data/srim2013`) |
| `--fix-trim-bugs` | Correct every bug of TRIM.exe that SRIM mode reproduces (results then depart from TRIM's; `docs/SRIM_MODE.md` §7). Applied before `--set`, so `--set physics.trim_bugs.NAME=true` keeps one of them (e.g. `zero_energy_defaults_skip_first_atom`, the only way the first atom keeps Ed/Eb/Es = 0) |
| `--fix-trim-bug NAME` | Correct one of them (repeatable): `stale_primary_energy`, `recoil_energy_created`, `struck_atom_direction`, `replacement_interstitial`, `subthreshold_sputter_counted`, `ion_phonons_counted_twice`, `cascade_random_replay`, `zero_energy_defaults_skip_first_atom`, `integer_layer_z`, `exyz_skips_first_step`, `recoil_event_counted_twice`, `single_precision_tables` |
| `--srim-data DIR` | SRIM `Data` folder: SNUC03 nuclear stopping (`stopping` command), and the POTCOEF03 pair potentials when `physics.pair_potentials` is `recoils` or `all` (TRIM.exe itself never reads them) |
| `--set KEY=VALUE` | Override any scenario key (repeatable), e.g. `--set "layers.0.elements.0.Ed=25 eV"` |
| `--save-config FILE` | Write the final scenario (after options and overrides) |
| `--json` | Print the summary as JSON |
| `--particles FILE` | Particle source: a TRIM.DAT list (ions, or recoils with `--source recoils`), or BACKSCAT/TRANSMIT/SPUTTER/TRIMOUT.txt (recoils) |
| `--source beam\|ions\|recoils` | Kind of particle source (TRIM's cascade codes 1-3, 4-5, 6-7) |
| `--exyz ENERGY` | Write EXYZ.txt: the ion's position every ENERGY of energy loss |
| `--formats LIST` | Output formats: any of `trim,csv,json,parquet` (default `trim,csv,json`) |
| `--output-units LIST` | Output units per dimension, e.g. `length=nm,energy=keV` or just `nm,keV` |
| `--number-format F` | `general`, `scientific`, `engineering`, `fixed` or `shortest`, optionally `:digits` (e.g. `scientific:4`) |
| `--output-spec FILE` | JSON file with the full output settings: tables, columns, rows, units, formats (docs/OUTPUT.md) |
| `--fast` | Fast approximations (`physics.fast_math`): about 2x faster SRIM mode, statistically the same results but not number for number; leave it off for TRIM comparisons (`docs/SRIM_MODE.md` §8.2) |
| `--trajectory-file F` | Write every recorded trajectory point to F (`.csv` or `.parquet`; relative names go into the output directory): ion, point, track, parent, kind, position, direction, energies, time. Records every ion unless `output.trajectories` gives a number. Not a TRIM output; `docs/OUTPUT.md` "Trajectory file" |
| `--trajectory-points P` | `all` (default) or `births` (the first point of every track only) |
| `--md-events F` | Hybrid BCA-MD recoil events to F (`.csv` or `.parquet`); needs `--damage full` (or `monolayer`) and `--md-ec`. `docs/OUTPUT.md` "MD event file" |
| `--md-ec LIST` | Transition energies E_c, ascending: `"1 keV"`, `"0.5,1,2,5 keV"` or `"500 eV,1 keV"` (repeatable) |
| `--md-emin E` | Recoils below E are not events (default `Ed`: the struck element's displacement energy) |
| `--md-handoff` | Handoff mode: the BCA stops following events (one E_c). Not TRIM: the tables then cover the BCA part only. Use `--fix-trim-bugs` for MD coupling |
| `--auto-depth` | Automatic target depth (`auto_depth.enabled`): before the run, quick test runs size the last layer so that the deepest test ion stops in its last 10 % (see "Checks and warnings" below) |
| `--log FILE` | Write the diagnostic log (configuration, threads, files written, warnings, errors) to FILE; crash reports go to the same folder |
| `-q, --quiet` / `-v, --verbose` | No progress bar / engine log messages and the diagnostic log (with debug lines) on stderr |

Ctrl+C stops a run after the ions in progress; the partial results are
written. The exit status is 0 on success, 1 on errors, 2 on usage errors,
130 after an interrupted run and 3 after a crash.

### Checks and warnings

`validate`, `run` and `trim` refuse values outside what the physics of
SRIM mode can simulate, with a message that names the key, the value, the
allowed range and what to use instead (the GUI shows the same messages when
a field is edited):

| Value | Accepted |
|---|---|
| Ion energy | above the stop energy (`physics.final_energy`, 1 eV), at most 100 GeV (where the stopping tables and their Bethe-Bloch continuation end) |
| Ion and atom masses | 0.9 Z to 3 Z + 10 amu (an isotope of the element) |
| Beam angle | 0 to below 90° to the surface normal |
| Layer thickness | at least one atomic layer, N^-1/3 (e.g. 2.28 Å for iron) |
| Density | positive, at most 50 g/cm3 |
| Compound correction | above 0, at most 10 (1 = Bragg's rule) |
| Ed | 0 to 10 keV; Eb and Es: 0 to 100 eV |
| Final and stop energies | at most 1 keV (TRIM: 1 eV) |
| Depth bins, threads | 1 to 10 000; at most 4096 (0 = all cores) |
| Files | `stopping.table` must exist, `srim_data` must be a folder |

Numbers must be whole where a count is expected, finite, and in range (no
wrap-around of a negative or huge count), also in TRIM.IN and TRIM.DAT
files; a TRIM.DAT line that cannot be transported is skipped with a warning.

Values that are allowed but probably not meant are printed as warnings
(`validate` and before a run; the GUI shows them in the status bar):
a density far from the material's (single elements: outside half to twice
the element's; compounds: a third to three times the mixture of the
elements' densities), a gas flag that does not fit the density, a target
more than 20 times as thick as the ions travel (1 keV ions in 1 m: the
depth tables would be meaningless; the message suggests a thickness), a
target so thin that most ions are transmitted (a note: fine for
transmission studies), and relativistic projectiles (classical collision
kinematics). The reach is the ions' path length from the stopping powers
(electronic + ZBL nuclear, through the layers in order).

**Automatic target depth** (`--auto-depth`, `auto_depth.enabled`): quick
test runs (Kinchin-Pease, fast approximations, `auto_depth.test_ions` ions,
default 1000) resize the **last** layer until no test ion is transmitted and
the deepest stops in the last `auto_depth.margin` (default 10 %) of the
target: the layer grows (doubles) while ions are transmitted and shrinks to
deepest/(1 - margin) when the deepest ion stops shallower. The first guess is
R (1.2 + f_n), R the path length and f_n the nuclear share of the stopping
along it (nuclear collisions straggle the range most); from 1 keV to 1 GeV
and from H to U it took 2 test runs (about 0.3 s). The new thickness is
printed (`Note: ...`) and used for the run (`--save-config` writes it); a
depth window that ended beyond the new target is reset to the whole target.
It is not used for recoil sources (TRIM.DAT atoms start at given depths) or
when the ions stop above the last layer. Leave it off to study transmission
through a given target.

**Transmitted ions.** During a SRIM-mode beam run, the first ion that
leaves through the back prints a warning that suggests `--auto-depth` (the
run continues), and the end of the run says how many escaped;
`--set auto_depth.warn_transmitted=false` turns this off.

**Crash reports.** If the program crashes (invalid memory access, stack
overflow, an exception escaping a thread, `abort`), it writes a report to
stderr and to `<log folder>/amcset-crash-<time>-<pid>.txt` (without
`--log`: `%TEMP%\amcset-crashes` or `/tmp/amcset-crashes`): the symbolised
stack trace, the calculation (engine, seed, ions), what every worker thread
was doing — for SRIM mode the ion and the command-line options that rerun
exactly that ion (`--set first_ion=B -n K`; the blocks of 16 ions start
there) — and the last log lines. On Windows a minidump (`.dmp`, for
Visual Studio or WinDbg) is written next to it (`AMCSET_MINIDUMP=0` turns
it off). `amcset diag-crash segv|abort|terminate|thread-exception|stack-overflow`
crashes on purpose, to check the reports on a machine.

## Output (SRIM mode)

The output directory receives TRIM's file set in TRIM's formats (format
`trim`), so that existing SRIM post-processing scripts work unchanged, and
the same results as tables in CSV (`csv/`), JSON (`results.json`) and
Parquet (`parquet/`) files. Which formats, tables, columns, rows, units and
number formats are written is set with `--formats`, `--output-units`,
`--number-format`, `--output-spec` or the scenario's `output.files`; see
[OUTPUT.md](OUTPUT.md). The TRIM-format files:

| File | Contents |
|---|---|
| `RANGE.txt` | Range moments and the ion / recoil depth distributions ((atoms/cm3)/(atoms/cm2)) |
| `VACANCY.txt` | Vacancies per Å per ion (KP: by ions / by recoils; full: knock-ons + per element) |
| `NOVAC.txt` | Replacement collisions (full cascades) |
| `IONIZ.txt`, `PHONON.txt` | Energy to ionization / phonons by ions and recoils (eV/Å/ion) |
| `E2RECOIL.txt` | Energy given to recoils and absorbed per element |
| `LATERAL.txt` | Lateral spread vs depth |
| `TDATA.txt` | Calculation parameters and totals |
| `BACKSCAT.txt`, `TRANSMIT.txt`, `SPUTTER.txt` | Exiting ions/atoms (energy, position, direction) |
| `TRIMOUT.txt` | All exiting ions and atoms (backscattered, transmitted, sputtered) in one file |
| `RANGE_3D.txt` | Final position of every ion |
| `COLLISON.txt` | Primary collisions and cascade sums of the first ions (`output.collision_details`) |
| `EXYZ.txt` | The ion's energy, position, electronic stopping and last recoil energy at every energy step (`output.exyz_interval`) |
| `PKA_SPECTRUM.txt` | Energy spectra of primary and all recoils (AMCSET addition) |
| `summary.json` | All totals and moments in one machine-readable file (always written) |

The native engine writes its histogram CSV files and `tallies.json`.

A TRIM.IN with cascade type 4-7 is run with the TRIM.DAT in the same folder:
`amcset trim TRIM.IN` (as TRIM, every line of TRIM.DAT is one particle and
the ion number of TRIM.IN is ignored).

## Scenario files

`amcset template` prints a complete scenario with every key; omitted keys
take their defaults and unknown keys are an error (so typos do not pass
silently). Relative paths are resolved against the scenario file's folder.

```json
{
  "mode": "srim",
  "title": "Fe (1 MeV) into Fe",
  "ion": {"element": "Fe", "energy": "1 MeV", "angle": "0 deg"},
  "layers": [
    {"name": "Fe", "thickness": "1.2 um", "density": "7.8658 g/cm3",
     "elements": [{"element": "Fe", "stoichiometry": 1.0,
                   "Ed": "40 eV", "Eb": "0 eV", "Es": "4.34 eV"}]}
  ],
  "ions": 1000,
  "damage": "kinchin-pease",
  "stopping": {"model": "srim", "table": "srim2013_stopping.amcset"},
  "srim_data": "C:/SRIM-2013/Data",
  "output": {"directory": "output"}
}
```

Keys:

- `ion`: `element` (symbol, name) or `z`; `mass` (default: most abundant
  isotope, as TRIM); `energy`; `angle` (to the surface normal). `spot`
  (written only when set; crystal transport only): where the ions meet the
  surface, see "Crystals".
- `layers[]`: `name`, `thickness`, `density`, `gas`, `compound_correction`
  (multiplies Bragg's-rule stopping, TRIM's "Compound Corr"), `elements[]`
  with `element`/`z`, `stoichiometry` (or `fraction`), `mass` (default:
  standard atomic weight), `Ed` (displacement), `Eb` (lattice binding), `Es`
  (surface binding). Defaults: Ed 25 eV, Eb 3 eV, Es = heat of sublimation
  (TRIM's defaults). `crystal` (optional): the layer's unit cell,
  orientation, temperature and thermal vibrations; it then decides the
  layer's composition and density, and `elements` may be left out. Followed
  by the crystal transport only; see "Crystals" below.
- `ions`, `seed`, `threads`, `damage` (`kinchin-pease`, `full-cascade`,
  `monolayer`). `first_ion` (default 0, written only when set): the number
  of the first ion of this run; every ion has its own random-number stream,
  so `first_ion = B` with `ions = K` repeats ions B .. B+K-1 of a longer run
  exactly when B is the start of a block of 16 (crash reports print B).
- `source` (optional): `{"type": "ions" | "recoils", "file": "TRIM.DAT"}` or
  with `"particles": [{"element": "Fe", "energy": "50 keV", "x": "1000 A",
  "y": "0 A", "z_position": "0 A", "direction": [1, 0, 0]}, ...]`. Ions take
  their species from `ion` and `ion.energy` is the highest energy allowed
  (as in TRIM); recoils are target atoms set in motion without an ion
  (`full-cascade`: TRIM's code 6, each atom starts a cascade; `monolayer`:
  code 7, each atom moves like an ion in monolayer steps).
  Every particle is run unless `ions` limits the number.
- `stopping`: `model` (default `srim`) and `table` (default: the SRIM-2013
  tables that come with AMCSET); `srim_data`.
- `physics`: the SRIM-mode model switches described in `docs/SRIM_MODE.md` §3:
  `transport` (`trim2013`, the default, follows TRIM.exe step by step;
  `refined` is AMCSET's earlier transport; `crystal` follows the lattice of
  layers that have a crystal, with its settings in `physics.crystal`, see
  "Crystals" below; not TRIM), pair potentials, scattering
  integral, straggling, final and stop energies, and the options of the
  refined transport. The defaults reproduce TRIM. `trim_bugs`: `"trim"`
  (the default: reproduce TRIM.exe's bugs), `"fixed"`, or an object of
  flags by name (`true` = as TRIM), e.g.
  `{"struck_atom_direction": false}`; see `docs/SRIM_MODE.md` §7. The
  crystal transport takes none of them, whatever is set.
  `fast_math` (default `false`): the fast approximations of
  `docs/SRIM_MODE.md` §8.2 (tabulated scattering angle, stopping powers and
  power functions; cheaper sine, cosine and rotations) — about 2x faster, statistically the same results, but not
  number for number, so not a step-by-step reproduction of TRIM.exe.
- `output.trajectories`: also `"all"` (every ion). `output.trajectory_file`
  (written only when set): `{"file": "traj.parquet", "points": "all" |
  "births", "min_energy": "0 eV"}`. `output.md_events` (written only when
  set): `{"file": "events.csv", "E_c": "1 keV" | ["0.5 keV", "1 keV"],
  "E_min": "Ed" | "40 eV", "mode": "record" | "handoff"}`. Both need the
  trim2013 transport; md_events needs full cascades or monolayer steps, 1 to
  32 E_c in increasing order, E_min below the smallest E_c and one E_c in
  handoff mode. Every column and rule: `docs/OUTPUT.md`.
- `auto_depth` (written only when not the defaults): `enabled` (default
  `false`), `test_ions` (1 to 10 000 000, default 1000), `margin` (0.01 to
  0.5, default 0.1), `warn_transmitted` (default `true`); see "Checks and
  warnings" above.
- `output`: `directory`, `depth_window` (TRIM's plot window, `["0 A","0 A"]`
  = whole target), `bins` (100 like TRIM), `range_3d`, `exit_lists`,
  `max_exit_records`, `collision_details` (number of ions in COLLISON.txt),
  `collision_recoils` (COLLISON.txt also lists every displacement and
  replacement of their cascades, TRIM's Collisions = 2),
  `trajectories` (ions whose cascades are kept for the GUI views),
  `exyz_interval` (EXYZ.txt energy step, e.g. `"10 keV"`),
  `trim_single_precision` (default `true`: the count tables of the
  TRIM-format files are summed in single precision, as TRIM.exe does; see
  docs/SRIM_MODE.md §5; always exact with the crystal transport), `files`
  (the output settings of docs/OUTPUT.md).
- Target elements take an optional `trim_atom`: their number in TRIM's list
  of target atoms (kept from a TRIM.IN; it orders the per-element columns and
  decides TRIM's replacement of zero energies, docs/SRIM_MODE.md §2.5).
- Target elements take an optional `damage_efficiency`: the arc-dpa and rpa
  constants of the `dpa` table, all four required,
  `{"arc_b": -0.568, "arc_c": 0.286, "rpa_b": "1018 eV", "rpa_c": 0.95}`
  (b_arc ≤ 0, 0 < c_arc ≤ 1, b_rpa > 0, c_rpa > 0). Without it Fe, Cu, Ni,
  Pd, Pt and W use the values of Nordlund et al. (2018); other elements have
  no arc-dpa/rpa. `output.dpa_fluence` (default `1e15`): the fluence the
  `dpa` table is given for, ions/cm² as a number or a string with `/cm2`,
  `ions/cm2`, `cm-2`, `/m2`, `ions/m2` or `m-2`; above 0 and at most 1e30.
  See docs/OUTPUT.md "dpa".
- `native`: options of the native engine (electron beam, screening and
  electronic-stopping models, thresholds, histogram range and bins, electron
  parameter tables).

`--set` addresses keys with dots and array indices: `ion.energy`,
`layers.1.density`, `layers.0.elements.2.Es`, `physics.final_energy`,
`output.bins`, `native.screening_model`. A list is given as JSON where the
key holds a list: `--set "layers.0.crystal.orientation.surface=[1,1,0]"`.

## Crystals

A layer can carry a crystal: a unit cell, its orientation, its temperature
and the thermal vibrations of its atoms. **Only the crystal transport
(`physics.transport: "crystal"`) follows the lattice**, so that ions can
channel; the `trim2013` (default) and `refined` transports treat every layer
as amorphous, as TRIM does, with the unit cell's composition and density,
and warn that the crystal is not used. Crystals are not TRIM: a TRIM.IN has
none, and `amcset convert` to TRIM.IN drops them.

**State (November 2026).** Single crystals, a static lattice with thermal
vibrations, one configuration of the lattice per ion; the GUI's *Crystal*
tab edits them (docs/GUI.md). Polycrystals (grains generated or read from
files, "Polycrystals" below) run from a scenario; the GUI cannot edit them
yet. Not there yet: damage that builds up from ion to ion, TRIM.DAT particle
sources. The transport reproduces a measured channeling dip and measured
boron profiles in silicon within the spread of those data ("Against
measurements" below), but **the division of the electronic stopping between
the path and the collisions (`physics.crystal.local_stopping_fraction`) is
not calibrated**: the measurements at hand allow 0.25 to 1, and the depth of
a channeling tail depends on it (about 10 % over that range for 30 keV B in
Si, up to 50 % for 0.5 keV B). Peak positions and trends (direction,
temperature, a covering oxide, the width of the channeling dip) depend on it
much less. "The crystal transport" below says what it does and how it was
checked.

### `layers.N.crystal`

Give exactly one of a built-in structure, a cell file, or a lattice with
sites:

```json
"layers": [{
  "thickness": "1 um",
  "crystal": {
    "structure": "diamond", "a": "5.431 A", "sublattices": ["Si"],
    "orientation": {"surface": [0, 0, 1], "reference": [1, 1, 0], "twist": "22 deg"},
    "temperature": "300 K",
    "vibrations": [{"element": "Si", "debye_temperature": "645 K"}]
  }
}]
```

| Key | Meaning |
|---|---|
| `structure` | A built-in structure: `sc`, `bcc`, `fcc`, `hcp`, `diamond`, `zincblende`, `wurtzite`, `rocksalt`, `cesium_chloride`, `fluorite`, `perovskite`, `L1_2`, `L1_0` (`amcset cell --list` shows their sublattices). |
| `sublattices` | With `structure`: the elements of each sublattice, in the structure's order. One entry per sublattice: `"Fe"`, `{"Fe": 0.9, "Cr": 0.1}` (a substitutional alloy: each site holds Fe with probability 0.9), or `[{"element": "Fe", "fraction": 0.9}, ...]`. Fractions lie in (0, 1] and add up to at most 1; the rest of a site is vacant. |
| `interstitials` | With `structure`: atoms on its built-in interstitial sites, e.g. `[{"sites": "octahedral", "occupants": {"C": 0.01}}]`. `octahedral` and `tetrahedral` for bcc, fcc and hcp; `tetrahedral` for diamond, zinc blende and rock salt; `octahedral` for fluorite and L1_2. |
| `a` | With `structure`: the lattice constant. **Without it the cell takes its size from the layer's `density`** (or, for one element, the element's density). |
| `c_over_a`, `u` | With `structure`: c/a of hexagonal and tetragonal cells (default: the ideal 1.633 for hexagonal, 1 for tetragonal); the internal parameter of wurtzite (default 0.375). |
| `file`, `format` | A cell file, relative to the scenario: CIF, VASP POSCAR/CONTCAR, LAMMPS data, extended XYZ. `format` (`cif`, `poscar`, `lammps-data`, `extxyz`) when the name does not tell. See "Cell files" below. |
| `elements`, `lammps_types`, `atom_style` | With `file`: the elements of a POSCAR without a species line (`["Ga", "As"]`); LAMMPS atom types as elements (`{"1": "Fe", "2": "Cr"}`); the LAMMPS atom style when the file does not say. |
| `lattice`, `sites` | The cell written out. `lattice`: `{"a": "3 A", "b": ..., "c": ..., "alpha": "90 deg", "beta": ..., "gamma": ...}` (b and c default to a, the angles to 90 deg) or `{"vectors_A": [[ax, ay, az], [bx, by, bz], [cx, cy, cz]]}` in Angstrom. `sites`: a list of `{"label", "position": [x, y, z], "kind", "element", "fraction"}` or with `"occupants"` (as in `sublattices`; `[]` is an empty site). Positions are fractional and are wrapped into the cell. `kind` is `lattice` (default; `substitutional` means the same) or `interstitial`. |
| `name` | A title. |
| `orientation` | `{"surface": [h, k, l], "reference": [u, v, w], "twist": "0 deg"}` or `{"euler_deg": [phi1, Phi, phi2], "twist": ...}`; see "Orientation". Default: surface (0 0 1). |
| `temperature` | The crystal's temperature, `"300 K"` (default) or `"27 C"`; it sets the thermal vibrations. 0 to 10000 K. |
| `vibrations` | Per element, `{"element": "Si", "debye_temperature": "543 K"}` (the rms amplitude then follows from the Debye model, the temperature and the atom's mass) or `{"element": "Si", "rms": "0.075 A"}` (the rms displacement along one axis). Without an entry, a crystal of one element takes that element's Debye temperature from Gemmell's table if it is there (C, Al, Si, V, Cr, Fe, Ni, Cu, Ge, Nb, Mo, Rh, Pd, Ag, Ta, W, Pt, Au, Pb; the run says so); any other atom stands still. |

**The crystal decides the layer's composition and density.** The layer's
`elements` may be left out (they are then the unit cell's, with default
energies) or listed to set Ed, Eb, Es and masses; an element that is not in
the cell is an error. Their stoichiometry becomes the number of atoms per
unit cell (stoichiometries written in the scenario must agree with it). The
layer's density becomes the unit cell's; a `density` in the scenario is
accepted if it agrees within 2 % (and replaced by the cell's), otherwise it
is an error naming both values. A gas layer cannot be a crystal.

Saved scenarios hold the cell written out (`lattice.vectors_A` and `sites`),
whatever it was made from, so that a scenario does not depend on other
files.

### Orientation

The laboratory frame is TRIM's: x is the depth, y and z lie in the surface,
and the beam tilt `ion.angle` turns the beam from x towards y.

- `surface`: the Miller indices (h k l) of the surface plane. Its normal is
  the depth axis; for non-cubic cells this is the plane's normal (from the
  reciprocal lattice), not the direction [h k l].
- `reference`: a direction [u v w] whose part in the surface points along
  the lab y axis, the direction the beam tilts towards. Without it: the
  lattice vector a (or b when a is along the normal).
- Hexagonal cells also take four indices: (h k i l) and [u v t w].
- `euler_deg`: Bunge Euler angles (phi1, Phi, phi2) instead of the indices,
  as texture and EBSD programs give them, relative to a sample frame with
  x_s and y_s in the surface and z_s the outward normal; lab x = -z_s,
  lab y = x_s, lab z = -y_s.
- `twist`: a rotation of the crystal about the surface normal (right-handed
  about the inward normal), applied last; the wafer's rotation on its holder.

To send the beam along a crystal direction, `amcset cell ... --beam "u v w"`
prints the `ion.angle` and `orientation.twist` to use; `--tilt ANGLE` prints
the crystal direction a tilted beam travels along.

### Cell files

| Format | What is read | Limits |
|---|---|---|
| CIF | The first data block: cell lengths and angles, the symmetry operations (`_symmetry_equiv_pos_as_xyz` or `_space_group_symop_operation_xyz`) applied to the atom sites, type symbols or labels (`Fe3+`, `O2-`, `Si1`), occupancies, fractional or Cartesian coordinates. Sites of several elements at one position become one site with several occupants. | A CIF that names a space group other than P 1 **without listing its operations** is refused (export it with the operations, or in P 1). Anisotropic displacement parameters, magnetic moments and the like are ignored. Split sites closer than 0.3 A are refused. |
| POSCAR / CONTCAR | Scale (or cell volume if negative), lattice, species line, counts, Direct or Cartesian coordinates, selective-dynamics flags (ignored). | VASP 4 files have no species line: give `elements`. |
| LAMMPS data | The box with tilt factors, `Masses`, `Atoms` in the styles atomic, charge, full, molecular, bond and angle. Atom types become elements by `lammps_types`, by a `# Fe` comment in `Masses`, or by the mass (within 0.1 amu of a standard atomic weight). | The whole box is the unit cell (it may be a supercell). With 6 or 9 columns and no `# style` comment after `Atoms`, give `atom_style`. |
| Extended XYZ | `Lattice="..."` and `Properties=` with `species` (or `Z`) and `pos`. | A plain XYZ file has no cell and is refused. |

Only CIF stores partial occupancy; the other writers refuse a cell that has
it. AMCSET's CIFs are in P 1 and mark interstitial sites with
`_amcset_site_kind`.

### The crystal transport

```json
{
  "ion": {"element": "B", "energy": "30 keV", "angle": "0 deg"},
  "layers": [{"thickness": "3 um",
              "crystal": {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"],
                          "orientation": {"surface": [0, 0, 1], "reference": [1, 1, 0]},
                          "vibrations": [{"element": "Si", "debye_temperature": "645 K"}]}}],
  "ions": 2000,
  "physics": {"transport": "crystal"}
}
```

This sends 30 keV boron down the ⟨001⟩ channels of silicon (mean range
2838 Å against 1128 Å in amorphous silicon); with `"angle": "7 deg"` and
`"twist": "22 deg"` it is the usual "random" implant (1150 Å). On the
command line: `--set physics.transport=crystal`.

What it does, after M. T. Robinson and I. M. Torrens, Phys. Rev. B 9, 5008
(1974) (the MARLOWE program):

- An ion or recoil moves on straight segments. Its next collision partners
  are the lattice atoms ahead of it whose impact parameter is below a
  maximum; each atom is displaced from its site by its thermal vibration (a
  Gaussian along each axis, uncorrelated between atoms, frozen for the whole
  history of one ion). Sites shared by several elements or partly vacant
  are filled at random, also once per ion.
- Atoms whose closest approaches nearly coincide along the path are struck
  together: each collision is computed as if alone and the momentum changes
  are added. This is what steers a channeled ion smoothly between rows.
- The scattering is AMCSET's (ZBL universal potential, magic formula or
  `physics.scattering: "integral"`), the electronic stopping AMCSET's
  stopping power, divided: a fraction is lost in the collisions, falling off
  exponentially with the distance of closest approach, the rest along the
  path. A channeled ion stays far from the atoms and loses less. The
  division is normalised so that a random direction loses the full stopping
  power.
- Struck atoms with more than Ed leave their site and are followed through
  the same lattice (`damage: "full-cascade"`), or the quick Kinchin-Pease
  estimate is booked (`"kinchin-pease"`). A projectile that keeps less than
  the smaller of its cut-off energy and Ed after displacing an atom of its
  own element takes that atom's site (a replacement).
- Amorphous layers (no `crystal`) take one random partner per monolayer
  step, so a target can mix both, e.g. an oxide on a crystal.

| Key | Meaning |
|---|---|
| `physics.transport` | `"crystal"`. The other transports ignore crystals. |
| `physics.crystal.max_impact_parameter` | The largest impact parameter of a collision. `"auto"` (default): for each crystal, the radius of the largest sphere that fits between its atoms (a√3/4 = 2.35 Å in Si, a/2 in fcc, a√5/4 in bcc), the smallest value for which no part of the lattice is free of collisions. Larger values change little. |
| `physics.crystal.simultaneous_distance` | Atoms whose closest approaches lie within this distance along the path are struck together (default `"0.5 A"`). 0 treats every collision in turn. |
| `physics.crystal.local_stopping_fraction` | The fraction of the electronic stopping lost in the collisions (default 0.5; 0 = all along the path, 1 = all in the collisions). **This decides the depth of channeling tails and is not calibrated**: the measurements under "Against measurements" allow 0.25 to 1. |
| `physics.crystal.local_stopping_decay` | The loss in a collision falls off as exp(−decay · R/a), R the distance of closest approach, a the screening length (default 0.3). |
| `physics.crystal.thermal_vibrations` | `false`: atoms stand exactly on their sites, whatever the crystals say (default `true`). |
| `ion.spot` | Where the ions meet the surface: `{"shape": "rectangle", "width": "1 um", "height": "1 um"}` (the default), `"circle"` (`"diameter"`), `"gaussian"` (`"fwhm"`, or `width`/`height` as full widths at half maximum) or `"point"`. Each ion enters at a random point of the spot; with `"point"` every ion meets the lattice at the same place (with a static lattice and no straggling: the same trajectory, repeated). Positions in the results stay relative to each ion's entry point. Only the crystal transport uses it. |

Atoms need their thermal vibrations (`crystal.vibrations`). A crystal of one
element that Gemmell's table lists takes its tabulated Debye temperature
when none is given, and the run says so; any other atom without them stands
still, channeling is then stronger than in a real crystal, and the run warns
about it. The built-in material data are Tables IV and V of D. S. Gemmell,
Rev. Mod. Phys. 46, 129 (1974): structure and room-temperature lattice
constant of 19 cubic elements and 24 cubic compounds, and the Debye
temperatures of the elements (`crystal::reference_crystals()`,
`common/include/crystal/materials.h`; the GUI offers them as *Common
crystals*). Nothing else is built in.

**No TRIM bugs.** The crystal transport is not TRIM and takes none of
TRIM.exe's bugs (docs/SRIM_MODE.md §7): `physics.trim_bugs` and
`output.trim_single_precision` do not apply to it, whatever they say. A
scenario keeps those switches for the other transports; a crystal run
records them as corrected in results.json, and its TRIM-format tables hold
the exact sums. TRIM's replacement of zero energies (Ed, Es, Eb = 0 become
20, 4 and 2 eV) belongs to the `trim2013` transport and is not applied
either: here a zero stays a zero.

**Cut-off energy and the number of vacancies.** With TRIM's cut-off (atoms
followed down to 1 eV, the default) few collisions count as replacements
and full cascades give TRIM-like vacancy numbers, well above NRT. With the
cut-off at Ed (`physics.recoil_stop_energy` and `physics.ion_stop_energy`
set to Ed) the displacement model is Kinchin and Pease's, as in the paper,
and the count is close to NRT: 30 keV B into Si (Ed = 15 eV), 284.9
vacancies per ion against NRT's 256.0, and 420.1 with the 1 eV cut-off
(TRIM's own full cascades in amorphous Si: 379.8).

**Checks** (`tools/crystal/channeling_checks.py`; SRIM-2013 stopping; the
numbers are in `build/crystal_checks/results.txt` on the workstation). They
test trends and internal consistency, not agreement with experiments.

| Check | Result |
|---|---|
| 30 keV B into Si at 300 K: mean range | ⟨001⟩ 2838 Å, ⟨110⟩ 4407 Å, 7° tilt / 22° twist 1150 Å; amorphous 1128 Å (`trim2013`), 1120 Å (this transport without a crystal) |
| ... vacancies per ion (quick damage) | ⟨001⟩ 43, ⟨110⟩ 14, tilted 156, amorphous 152 |
| A 100 Å oxide on Si ⟨001⟩ | 1522 Å: the oxide scatters part of the beam out of the channels |
| Temperature (⟨001⟩ range) | static lattice 3714 Å, 20 K 3220 Å, 300 K 2838 Å, 1200 K 1972 Å |
| Division of the stopping (⟨001⟩ range) | all on the path 2598 Å, half 2838 Å, all in the collisions 3132 Å; in the tilted direction 1154, 1150 and 1139 Å (no dependence) |
| Numerical parameters (⟨001⟩ range) | maximum impact parameter 3 Å instead of 2.35 Å: 2855 Å; sequential collisions: 2839 Å; simultaneity distance 1.5 Å: 2906 Å |
| 1 MeV H through 3000 Å of Si, random direction: energy to electrons | 1.2539 % (all on the path), 1.2555 % (all in the collisions), 1.2524 % (`trim2013`): the division is normalised |
| 1 MeV He through 600 Å of W along ⟨001⟩: atoms displaced by the ion against the tilt | 2.0 % of the amorphous value on the axis, half of it at about 1.73° (0.82 of Lindhard's ψ₁ = 2.10°), 1.25× at 2.5-3° (the shoulder) |
| Random crystal orientations (16), least-channeling orientation against amorphous range | 20 keV Fe in Fe 83 against 89 Å; 30 keV B in Si 1136 against 1129 Å; the averages are deeper (188 and 1237 Å) because some orientations channel |

**Against measurements** (`tools/crystal/calibration.py`; SRIM-2013
stopping; a perfect crystal for every ion; the numbers are in
`build/crystal_calibration/results.txt` on the workstation).

- *Thermal amplitudes.* The Debye model reproduces the tabulated u₁(293 K)
  of all 19 elements of Gemmell's Table IV from their Debye temperatures
  (test `CrystalMaterials`).
- *Channeling dip: 1 MeV He along ⟨110⟩ of Ge at 293 K* (atoms displaced by
  the ion at depths of 200-1200 Å, relative to amorphous Ge). Minimum yield
  0.022 and half-angle 0.99° = 0.80 ψ₁, against the measured 0.03 and 0.95°
  (Picraux et al., as quoted by Gemmell) and Barrett's formulas, 0.025 and
  0.93°. The shoulder reaches 1.57 at 1.7-2.0°.
- *Boron in (100) Si, as-implanted SIMS profiles*, for five divisions of
  the stopping:

| Measurement | Measured | Crystal transport, fraction 0 / 0.25 / 0.5 / 0.75 / 1 | `trim2013` (amorphous) |
|---|---|---|---|
| 15 keV, 7° tilt, 3·10¹⁵ cm⁻² (Jung 2010, Table 2): Rp (nm) | 64.2 | 62.1 / 62.4 / 62.2 / 62.2 / 62.5 | 60.5 |
| ... ΔRp (nm) | 31.8 | 32.1 / 32.7 / 33.1 / 33.5 / 34.1 | 25.7 |
| ... skewness | 0.98 | 0.59 / 0.68 / 0.76 / 0.85 / 0.94 | 0.00 |
| ... kurtosis | 5.56 | 3.49 / 3.78 / 4.07 / 4.40 / 4.80 | 2.52 |
| 30 keV, 1.5·10¹⁴ cm⁻² (Shao 2003, Fig. 33): depth of 10¹⁸ cm⁻³ (nm) | about 210 | 221 / 220 / 220 / 221 / 221 | 201 |
| ... depth of 10¹⁷ cm⁻³ (nm) | about 290 | 277 / 294 / 293 / 312 / 321 | 226 |
| 0.5 keV, 10¹⁴ cm⁻² (Schroer 2000, Fig. 1): depth of 10¹⁹ cm⁻³ (nm) | 14.2 | 14.5 / 14.6 / 14.9 / 15.1 / 15.3 | 9.7 |
| ... depth of 10¹⁸ cm⁻³ (nm) | 25.2 | 18.6 / 19.4 / 20.2 / 21.0 / 21.8 | 12.0 |
| ... depth of 10¹⁷ cm⁻³ (nm) | 35 | 20.5 / 21.8 / 23.3 / 25.5 / 27.9 | 13.8 |
| ... depth of 10¹⁶ cm⁻³ (nm) | 57 | 21.9 / 23.5 / 26.4 / 29.9 / 35.5 | 15.1 |

  The papers: Jung, Trans. Electr. Electron. Mater. 11, 120 (2010); Shao et
  al., Mater. Sci. Eng. R 42, 65 (2003); Schroer et al., Mater. Sci. Eng. B
  71, 219 (2000). The 15 keV moments are those of the profile to 450 nm, as
  in the paper, which gives 59.7 nm, 30.1 nm, 0.72 and 4.62 for UT-Marlowe;
  it does not state the twist (22° is assumed; 45° gives Rp 64.9 nm, and 0°,
  a planar channel, 90.7 nm). The 30 keV and 0.5 keV values were read off
  the published figures (a few per cent of the depth), which state neither
  tilt nor twist: 7°/22° is assumed at 30 keV and 0° at 0.5 keV (7°/22°
  gives 10.6 nm for the first 0.5 keV level; 10.1 nm was measured at
  10¹⁵ cm⁻²).

What this shows. The peak and the width are reproduced whatever the
division: Rp to 3 %, ΔRp to 1-7 %, the first decade of the 0.5 keV profile to
2-8 %, where the amorphous calculation misses the width by 19 % and the
0.5 keV depth by 32 %. The division itself is not fixed by these data: the
30 keV tail prefers 0.25-0.5; the skewness and kurtosis at 15 keV and the
deep 0.5 keV tail prefer 1. The unknowns of the measurements are as large as
its effect: tilt and twist; the damage of the earlier ions at 10¹⁵ cm⁻² and
above, which dechannels and which the perfect crystal of the simulation
lacks; and the deep SIMS tail at 0.5 keV, measured to 57 nm where no
division gives more than 36 nm. The default therefore stays 0.5. Fitting it
needs low-dose implants with stated tilt and twist.

Other limits: one frozen lattice per ion (the atoms an ion displaces stay
displaced for that ion only); a site emptied in a history is not refilled by
the atom that replaces it; at an interface the atoms of the next layer are
met only after crossing into it; local electronic losses use the layer's
mean stopping per atom (exact for one element); no COLLISON.txt, trajectory
file, MD events or dpa table (those need `trim2013`); `auto_depth` and
TRIM.DAT sources are refused. It is slower than `trim2013`: 2000 ions of
30 keV B into Si took 0.35 s against 0.011 s with quick damage and 1.7 s
against 0.16 s with full cascades (all cores of the workstation, which was
not idle; `elapsed_seconds` of the check runs).

### Polycrystals

A crystal with a `polycrystal` object is a layer of grains: the same unit
cell in many orientations. The crystal transport finds the grain a moving
atom is in and follows that grain's lattice; the atoms of the next grain
are met after crossing into it.

```json
"crystal": {
  "structure": "bcc", "a": "2.8665 A", "sublattices": ["Fe"],
  "polycrystal": {"grain_size": "50 nm", "texture": "random", "boundary_width": "5 A"}
}
```

| Key | Meaning |
|---|---|
| `source` | `"voronoi"` (the default: grains are generated), `"file"` (AMCSET's grain file), `"neper"` (a Neper `.tess`), `"ebsd"` (a TSL `.ang` or Oxford `.ctf` map). With `file` and no `source` the file's name decides. |
| `grain_size` | Generated grains: the mean grain diameter (default `"500 A"`): that of the sphere of the mean grain volume, or of the circle of the mean area for columns. 10 A or more. |
| `columnar` | Generated grains: `true` makes them columns through the layer. |
| `texture`, `spread` | Generated grains: `"random"` (the default: uniformly distributed orientations), `"fiber"` (the crystal's `orientation.surface` plane in every grain, turned at random about the normal; `spread` is the rms tilt of that axis) or `"mosaic"` (the crystal's `orientation`, turned by an rms angle `spread` about random axes). |
| `seed` | Of the grain generator (default 1). Not the run's seed: the same grains for every run. |
| `period` | The lateral period of the grains, one length or `[y, z]`: a representative patch of the surface that repeats. Default: eight grain sizes for generated grains, the extent of the file's map otherwise. |
| `boundary_width` | Width of a slab at every grain boundary that is treated as amorphous (random collision partners): a disordered boundary. Default 0; at most 100 A. |
| `file` | The file of the grains, relative to the scenario. |
| `scale` | Neper: the length of one unit of the tessellation (default `"1 um"`). |
| `misorientation`, `phase` | EBSD: neighbouring points closer in orientation than this (default `"5 deg"`, up to the lattice's own rotations) are one grain; `phase` selects the phase of the map that is this crystal (default 0: every indexed point). |

**The beam must cover the grains.** Each ion enters at a random point of
`ion.spot` (default 1 um square); with a spot smaller than the `period` only
part of the grains is met. Set the spot to the period.

**Sources.** Generated grains are a Poisson-Voronoi tessellation (random
seed points; each grain's lattice passes through its seed). The grain file
is a text file of seed points with Bunge angles, or of voxels; its format,
and how Neper tessellations (seeds with Laguerre weights, `*ori`) and EBSD
maps (points grouped into grains, which become columns into the depth) are
read, is described in `common/include/crystal/grain_io.h`. Orientations in
files are Bunge angles relative to the sample (x and y in the surface, z the
outward normal), as `orientation.euler_deg`; the crystal's own `orientation`
is then not used.

**Checks** (`tools/crystal/polycrystal_checks.py`; 30 keV B into Si along
the normal of (0 0 1), 4000 ions each; `build/polycrystal_checks/results.txt`
on the workstation). Trends, not measurements.

| Target | Mean range | Vacancies per ion |
|---|---|---|
| Single crystal ⟨001⟩; 7°/22°; amorphous | 2641 A; 1134 A; 1120 A | 52; 157; 152 |
| Random grains of 50 / 100 / 200 / 500 / 1000 / 3000 A | 1136 / 1153 / 1176 / 1199 / 1218 / 1253 A | 152 / 150 / 149 / 149 / 148 / 146 |
| Random columns of 200 A | 1300 A | 141 |
| Fiber texture (⟨001⟩ along the normal), grains of 500 A: spread 0 / 2° / 5° | 1743 / 1454 / 1272 A | 103 / 126 / 141 |
| Mosaic, grains of 500 A: spread 0 / 0.5° / 2° | 1740 / 1722 / 1514 A | 104 / 105 / 121 |
| Fiber, grains of 100 A, boundary width 0 / 5 / 10 / 20 / 40 A | 1366 / 1235 / 1188 / 1160 / 1125 A | 131 / 142 / 147 / 149 / 152 |
| 400 grains from a grain file; as many generated | 1181 A; 1176 A | 150; 149 |

Small random grains are the amorphous target, larger ones let some ions
channel; a texture that keeps the channel along the beam keeps part of the
channeling; wide disordered boundaries remove it.

**Limits.** The lattices of neighbouring grains are not in registry: even
grains of one orientation (mosaic without spread) dechannel part of the beam
at every boundary (1740 A against 2641 A for the single crystal above),
where real subgrain boundaries keep the lattice nearly continuous. The
readers of Neper and EBSD files were written from the formats and are
tested on small hand-written files, not on files produced by Neper or by an
EBSD system: check the orientation convention with a known grain before
relying on them. A tessellation is repeated laterally, so the cells at the
sides of a non-periodic Neper tessellation differ from Neper's. For voxels
the distance to a boundary counts faces only. The GUI shows a polycrystal's
unit cell but cannot edit or preview its grains.

### `cell`

```bash
amcset cell --list                                           # the built-in structures
amcset cell --structure bcc --sublattice Fe --a "2.8665 A"    # a summary
amcset cell --structure bcc --sublattice "Fe:0.9,Cr:0.1" \
            --interstitial octahedral=C:0.01 --density "7.8 g/cm3" --json   # for a scenario
amcset cell Si.cif -o POSCAR --to-format poscar               # read, check, convert
amcset cell Si.cif --surface "0 0 1" --reference "1 1 0" --beam "1 1 2"
```

| Option | Meaning |
|---|---|
| `<file>` | A cell file (formats above); `--from-format F` when the name does not tell; `--elements Ga,As` for a POSCAR without species |
| `--structure NAME` | A built-in structure, with one `--sublattice SPEC` per sublattice (`Fe` or `"Fe:0.9,Cr:0.1"`: site fractions, not normalised) and `--interstitial SITES=SPEC` |
| `--a LENGTH`, `--c-over-a X` | Lattice constant and c/a; without `--a`, `--density D` (or the element's density) sizes the cell |
| `--surface "h k l"`, `--reference "u v w"`, `--twist ANGLE` | The orientation |
| `--beam "u v w"` | Prints the tilt (`ion.angle`) and twist for a beam along this direction |
| `--tilt ANGLE` | Prints the crystal direction of a beam with this tilt |
| `-o FILE`, `--to-format F` | Writes the cell (format from the name) |
| `--json` | Prints the `crystal` object for a scenario instead of the summary |

The summary lists the cell parameters, volume, atoms per cell, atomic and
mass density (standard atomic weights), composition, the sites, the
neighbour shells of the first lattice site (8 at 2.48 A and 6 at 2.87 A
identify bcc iron) and the crystal directions of the laboratory axes.

## `stopping`

Prints electronic, nuclear and total stopping and the CSDA range on a
logarithmic energy grid (like SRIM's SR Module):

| Option | Meaning |
|---|---|
| `--ion X` | Projectile symbol or Z |
| `--ion-mass M` | Projectile mass (default: most abundant isotope) |
| `--target SPEC` | `Fe` or `"Fe:0.91,Cr:0.09"` (atom fractions, normalised) |
| `--density D` | Target density (default: the element's density; required for compounds) |
| `--gas` | Gas-phase target |
| `--from E --to E` | Energy range (default 10 keV .. 10 MeV) |
| `--per-decade N` | Points per decade (default 10) |
| `--model MODEL` | `srim`, `builtin`, `bethe-bloch`, `lss` (default `srim` with the bundled or the given tables; `builtin` if no tables are found) |
| `--units U` | `eV/A`, `keV/um`, `MeV/(mg/cm2)`, `eV/(1e15 atoms/cm2)` |
| `--csv` | Comma-separated output |

`--lammps` prints the table LAMMPS's `fix electron/stopping` reads (metal
units): column 1 the kinetic energy in eV, then the electronic stopping
N·S_e in eV/Å of each atom type, energies increasing, `#` comment lines on
top. `--ion Fe,Cr` gives one column per atom type in that order (each ion
in the target of `--target`/`--density`); `--ion Fe --ntypes 3` repeats one
column for three types. Electronic stopping only (LAMMPS's fix applies
electronic stopping; nuclear stopping is MD's own collisions), from the same
models as SRIM mode (`--model`, bundled SRIM-2013 tables by default).

```bash
amcset stopping --lammps --ion Fe,Cr --target "Fe:0.9,Cr:0.1" --density "7.8 g/cm3" \
               --from "10 eV" --to "10 MeV" --per-decade 20 > elstop.txt
```

Nuclear stopping uses the ZBL universal fit, or SRIM's SNUC03 pair fits
with `--srim-data`. Above the SRIM tables (and for `builtin` above ~1 MeV/u)
the electronic stopping is relativistic Bethe-Bloch with effective charge
and density-effect correction, so the energy range extends to 100 GeV.
