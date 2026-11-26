# Output files and how to shape them

A SRIM-mode run writes its results in up to four formats, chosen per run:

| Format | Files | Contents |
| --- | --- | --- |
| `trim` | `RANGE.txt`, `VACANCY.txt`, `NOVAC.txt`, `IONIZ.txt`, `PHONON.txt`, `E2RECOIL.txt`, `LATERAL.txt`, `TDATA.txt`, `BACKSCAT.txt`, `TRANSMIT.txt`, `SPUTTER.txt`, `TRIMOUT.txt`, `RANGE_3D.txt`, `COLLISON.txt`, `EXYZ.txt`, `PKA_SPECTRUM.txt` | TRIM's own text layout, so scripts written for TRIM keep working |
| `csv` | `csv/<table>.csv` | one file per table; headings `Label [unit]`; `# ...` comment lines with the description and the run data |
| `json` | `results.json` | everything in one file: run metadata, the scenario, the scalar results (value, unit, dimension) and every table column-wise |
| `parquet` | `parquet/<table>.parquet` | one Apache Parquet file per table (doubles, 64-bit integers, UTF-8 text); units, labels and run data in the file's key-value metadata |

`summary.json` (AMCSET's short summary) is always written. The default
formats are `trim`, `csv` and `json`. The files are written concurrently
(one task per file); their contents do not depend on it. A run with the
fast approximations (`physics.fast_math`) says so in `TDATA.txt`
("Fast approximations = on") and in `summary.json` (`"fast_math": true`);
without them neither appears, so the files are the same as before the
option existed.

## Tables

Every result is a table of typed columns (each numeric column has a physical
dimension and a unit). The tables, in TRIM's units:

| Table | Columns | TRIM file / plot |
| --- | --- | --- |
| `range` | `depth`, `ions`, `recoils_<element>` (1/cm = (atoms/cm³)/(atoms/cm²)) | RANGE.txt |
| `lateral` | `depth`, `lateral_range`, `lateral_straggling`, `radial_range`, `radial_straggling` (Å) | LATERAL.txt |
| `vacancies` | `depth`, `by_ions`, `by_recoils` (quick) or `knock_ons`, `vacancies_<element>` (full), `total` (1/Å per ion) | VACANCY.txt |
| `replacements` | `depth`, `replacements` | NOVAC.txt |
| `damage` | `depth`, `displacements`, `vacancies`, `replacements`, `interstitials` | collision-events plot |
| `dpa` | `depth`, `trim_vacancies`, `nrt`, `kinchin_pease`, `arc_dpa`, `rpa` (dpa at `output.dpa_fluence`), `damage_energy` (eV/Å per ion) | (AMCSET; trim2013 only, see "dpa" below) |
| `ionization` | `depth`, `by_ions`, `by_recoils` (eV/Å per ion) | IONIZ.txt |
| `phonons` | `depth`, `by_ions`, `by_recoils` | PHONON.txt |
| `energy_to_recoils` | `depth`, `from_ions`, `absorbed_<element>` | E2RECOIL.txt |
| `sputtering` | `element`, `z`, `yield`, `mean_energy`, `transmitted_recoils` | calculation window |
| `surface_energy` | `energy` (0.1–20 eV), `arrivals_<element>`, `yield_above_<element>` | integral / differential sputtering plots |
| `recoil_spectrum` | `energy_low`, `energy_high`, `primary_knock_ons`, `all_recoils` | (AMCSET) |
| `exits` | `kind` (B/T/S), `ion`, `z`, `element`, `energy`, `x`, `y`, `z_position`, `cos_x`, `cos_y`, `cos_z` | BACKSCAT/TRANSMIT/SPUTTER/TRIMOUT |
| `ion_positions` | `ion`, `x`, `y`, `z` | RANGE_3D.txt |
| `exyz` | `ion`, `energy`, `x`, `y`, `z`, `electronic_stopping`, `last_recoil_energy` | EXYZ.txt |
| `collisions` | `ion`, `ion_energy`, `x`, `y`, `z`, `electronic_stopping`, `atom`, `atom_z`, `recoil_energy`, `displacements`, `vacancies`, `replacements`, `interstitials` | COLLISON.txt |
| `ion_damage` | `ion`, `displacements`, `replacements`, `vacancies`, `interstitials`, `sputtered`, `transmitted` | COLLISON.txt (summary after each ion) |
| `cascade_atoms` | `ion`, `cascade`, `number`, `atom`, `atom_z`, `energy`, `x`, `y`, `z`, `vacancy`, `replacement` | COLLISON.txt with `output.collision_recoils` (TRIM's "Ion + Recoils") |

`depth` is the end of each depth bin, as in TRIM's files. Scalars (in
`results.json` and the GUI's statistics) include the range moments,
backscattered/transmitted fractions, sputtering yield, vacancies,
displacements, replacements and interstitials per ion and the energy
balance. The quick (Kinchin-Pease) calculation follows no recoils: its
displacements equal its vacancies, it has no replacements or
interstitials, and `knock_ons_per_ion` counts the atoms the ion itself
displaces.

`surface_energy` is the data behind TRIM's sputtering plots: the number of
target atoms per ion that reach the surface with a given normal energy
(E cos²θ, 0.1 eV bins as TRIM). `yield_above_<element>` at energy E is the
sputtering yield the target would have if its surface binding energy were E.

## Choosing what is written

The scenario's `output.files` object (or a separate file given with
`--output-spec`, or the GUI's *Output > Tables, units, numbers...* dialog)
controls the files:

```json
{
  "formats": ["csv", "json", "parquet"],
  "units": { "length": "nm", "energy": "keV", "energy_per_length": "keV/um" },
  "number_format": { "notation": "scientific", "digits": 4 },
  "csv": { "delimiter": ",", "comments": true },
  "json": { "pretty": true },
  "parquet": { "units_in_names": false },
  "tables": {
    "range": {
      "columns": ["depth", "ions"],
      "units": { "depth": "um" },
      "labels": { "ions": "Fe concentration per fluence" },
      "formats": { "ions": { "notation": "engineering", "digits": 3 } },
      "filters": [ { "column": "depth", "min": "100 nm", "max": 2 } ],
      "drop_zero_rows": true,
      "every": 2,
      "max_rows": 50,
      "file_name": "fe_range"
    },
    "exits": { "exclude": ["cos_y", "cos_z"] },
    "collisions": false,
    "*": { "number_format": "general" }
  }
}
```

- `formats`: any of `trim`, `csv`, `json`, `parquet`.
- `units`: the unit of every column of a dimension. Dimensions and units:
  - `length`: A (Å), nm, um, mm, cm, m
  - `inverse_length`: 1/A, 1/nm, 1/um, 1/mm, 1/cm, 1/m
  - `energy`: eV, meV, keV, MeV, GeV, J
  - `energy_per_length`: eV/A, eV/nm, keV/nm, keV/um, MeV/mm, MeV/cm, eV/cm, J/m
  - `inverse_energy`: 1/eV, 1/keV, 1/MeV, 1/J
  - `dimensionless`: 1, %, ppm
  - `angle`: deg, rad, mrad; `time`: s, ms, us, ns, ps, fs

  Conversions come from the Boost.Units definitions of the units, and a unit
  of the wrong dimension is rejected.
- `number_format`: `notation` is `general` (printf `%g`, the default, with
  12 significant digits), `scientific` (`1.2346E+04`), `engineering`
  (`12.346E+03`), `fixed` (decimals) or `shortest` (the shortest text that
  reads back exactly); `digits` are significant digits for `general` and
  decimals otherwise; `uppercase` selects `E` or `e`. It applies to CSV and
  JSON (TRIM text files keep TRIM's layout).
- `tables.<key>` (or `"*"` for all others; `false` leaves a table out):
  `include`; `columns` (which, in this order); `exclude`; `units`, `labels`
  and `formats` per column; `filters` (keep rows whose value lies in
  `[min, max]`, bounds as numbers in the column's output unit or as text
  with a unit); `drop_zero_rows` (rows whose columns after the first are all
  zero); `every` (every n-th row); `max_rows`; `number_format`;
  `file_name`.

On the command line:

```bash
amcset run scenario.json --formats csv,parquet --output-units nm,keV \
       --number-format scientific:4 --output-spec my_tables.json
```

## Reading the files

```python
import pandas as pd, pyarrow.parquet as pq, json
r = pd.read_csv("out/csv/range.csv", comment="#")          # 'Depth [A]', ...
t = pq.read_table("out/parquet/range.parquet")              # columns: depth, ions, ...
units = json.loads(t.schema.metadata[b"amcset.units"])      # {'depth': 'A', ...}
res = json.load(open("out/results.json"))                   # metadata, scalars, tables
```

The GUI opens `results.json` (File > Open Results) to show the plots and
tables of an earlier run.

## Crystals: orientation (SRIM mode, layers with a crystal)

Every run with a crystal layer (`layers.N.crystal`) reports how each crystal
sits in the target and meets the beam, whichever transport runs (`followed`
says whether the crystal transport follows the lattice). Not TRIM.

**summary.json `crystals`**, one object per crystal layer: `layer`,
`layer_name`, `cell`, `followed`, `top_A`; `lattice_vectors_A` (rows a, b,
c in the crystal frame); `rotation_crystal_to_trim_axes` (TRIM's axes: x =
depth, the beam tilting in x-y) and `rotation_crystal_to_output_axes` (the
axes of the trajectory and event files: x = TRIM's y, y = TRIM's z, z =
depth), both mapping a crystal-frame vector v to R v; `euler_bunge_deg`
(phi1, Phi, phi2 of the crystal relative to the sample frame: x, y in the
surface, z the outward normal, as `orientation.euler_deg`, so that they
give the same rotation back); the low-index crystal directions nearest the
depth, TRIM's y and TRIM's z (`{"uvw", "off_deg"}`) and the beam
(`beam_direction`, and `beam_unit_vector_crystal`); `tilt_deg`, `twist_deg`,
`temperature_K`; `vibrations` (per element the rms displacement along one
axis the transport uses, 0 for a static atom, and whether it came from
Gemmell's table); `axes_near_beam`, the six low-index axes (indices up to
2) nearest the beam with the angle to it, the mean spacing of the atoms
along the row (the row's period over its atoms: a/√2 along ⟨110⟩ of
diamond, a√3/2 along ⟨111⟩) and Lindhard's critical angle ψ₁ =
(2 Z₁ Z₂ e²/(E d))^½ for the ion at its energy (Z₂: the cell's mean atomic
number). An ion within about ψ₁ of an axis channels along it.

**Table `crystals`** (CSV, JSON, Parquet): one row per crystal layer with
`layer`, `crystal`, `followed`, `phi1`, `Phi`, `phi2`, `depth_direction`,
`beam_direction`, `beam_off`, `tilt`, `twist`, `temperature_K`, the nearest
axis (`nearest_axis`, `axis_angle`, `axis_spacing`, `psi1`) and the
rotation to the output axes `r11` … `r33`.

## dpa: displacements per atom by damage model (SRIM mode)

TRIM gives vacancies per ion and depth; dpa is usually computed from them
afterwards. AMCSET does it during the run, by several damage models, for
the `trim2013` transport (both damage modes). It is **not a TRIM output**
and changes nothing TRIM writes.

**What is booked.** Every *primary knock-on atom* (PKA: a target atom
struck by the ion, or the atom of a TRIM.DAT recoil source) with an energy
transfer T of at least 1 eV inside the depth window gets a damage energy
T_d, TRIM's Robinson partition with the layer's mean Z and M (the same T_d
the quick calculation uses, docs/SRIM_MODE.md §2.5). The models are then
applied to T_d with the PKA element's Ed and booked **at the depth where
the PKA starts** (not spread over its cascade):

| Column | Model | Per PKA |
| --- | --- | --- |
| `trim_vacancies` | TRIM's own vacancies (VACANCY.txt, summed over elements, exact sums): the quick calculation's NRT values, or the vacancies of the followed full cascades, booked where they are left | — |
| `nrt` | Norgett-Robinson-Torrens (ASTM E521) | 0 below Ed, 1 up to 2Ed/0.8, 0.8 T_d/(2Ed) above |
| `kinchin_pease` | Kinchin-Pease with T_d | 0 below Ed, 1 up to 2Ed, T_d/(2Ed) above (= 1.25 × NRT above 2.5 Ed) |
| `arc_dpa` | athermal-recombination-corrected dpa, Nordlund et al., Nat. Commun. 9, 1084 (2018), Eqs. 6-7 | NRT × ξ_arc, ξ_arc = (1 − c)/(2Ed/0.8)^b · T_d^b + c above 2Ed/0.8 |
| `rpa` | replacements per atom (atomic mixing), same paper, Eq. 9 | NRT × ξ_rpa, ξ_rpa = (b^c/(2Ed/0.8)^c + 1) · T_d^c/(b^c + T_d^c) |
| `damage_energy` | the PKAs' T_d per ion and Å | — |

dpa in a bin = (events per ion and Å) × fluence / (atomic density of the
layer at the bin centre), with the fluence `output.dpa_fluence` (ions/cm²,
default 1e15; dpa scale linearly with it, so the table at 1e15 can be
rescaled by hand). In the quick calculation `nrt` equals `trim_vacancies`
exactly (TRIM's quick vacancies are NRT).

**Constants.** arc-dpa and rpa need four material constants per element.
AMCSET has the published fits for Fe, Cu, Ni, Pd, Pt and W (Table 1 of the
paper, fitted to MD with Ed = 40, 33, 39, 41, 42 and 70 eV); other elements
take them from `layers.N.elements.M.damage_efficiency` (`{"arc_b": ,
"arc_c": , "rpa_b": "eV", "rpa_c": }`, docs/CLI.md), which also overrides
the published values. If any PKA belongs to an element without constants,
the `arc_dpa` and `rpa` columns are left out (and the table description
says why) rather than mixing models.

**summary.json** gets `damage_models`: `damage_energy_eV_per_ion`,
`nrt_per_ion`, `kinchin_pease_per_ion`, `arc_dpa_per_ion` and
`rpa_per_ion` (when available), `trim_vacancies_per_ion` (TRIM's total,
which in full cascades also counts vacancies outside the depth window) and
`dpa_fluence_per_cm2`.

**What it is and is not.**

- arc-dpa and rpa estimate the *primary* damage: Frenkel pairs and
  replaced atoms surviving the cascade's first picoseconds, from MD in
  pure metals. They do not include thermally activated recombination or
  migration, clustering, dose overlap or the defect structure; the paper
  calls arc-dpa possibly not applicable to non-metals that amorphize.
- The constants were fitted with specific Ed values. With another Ed the
  efficiency curve moves with 2Ed/0.8; use the fitted Ed when comparing
  with MD.
- T_d is TRIM's Robinson partition, applied once per PKA; MD uses
  electronic friction along the cascade (ZBL96 above 10 eV in the paper).
  They are not the same quantity; even within SRIM, T_d of the quick and
  the full-cascade calculation are reported to differ by up to ±25 %
  (Agarwal et al., as quoted by Mitsi et al. 2023).
- Booking at the PKA's starting depth puts high-energy PKAs' damage
  slightly too shallow (by about their range); `trim_vacancies` in full
  cascades is booked where the vacancies are.
- TRIM's full-cascade vacancies are a count of every atom moved above Ed,
  not of surviving defects: about 2× NRT and 4-20× arc-dpa for PKAs of
  1-100 keV in metals. Which model to use, measured against MD and
  experiments: docs/VALIDATION.md §5.
- The GUI plots `dpa` (all but rpa) and `rpa` and has the fluence on its
  Output page; it does not edit the constants.

## Trajectory file and MD event file (SRIM mode, opt-in)

Two optional files for analysing collision cascades and for coupling them
to molecular dynamics (the hybrid BCA-MD method). They are **not TRIM
outputs**: TRIM.exe has no equivalent. Every transport writes them
(`trim2013`, `refined` and `crystal`; "Other transports" below) and they are
written while the run goes on, block by block, so they
can be larger than the memory. Both are written in ion order for any number
of threads, and are identical for a given seed whatever the thread count.

| | Trajectory file | MD event file |
|---|---|---|
| Scenario key | `output.trajectory_file` | `output.md_events` |
| CLI | `--trajectory-file F`, `--trajectory-points births` | `--md-events F`, `--md-ec LIST`, `--md-emin E`, `--md-handoff` |
| GUI | Output tab, *Trajectory and MD event files* | the same group |
| Needs | `output.trajectories` > 0 or `"all"` (the CLI option sets `"all"` if it is 0) | full cascades (`damage` full or monolayer), 1 to 32 values of E_c |
| Size | ~89 000 rows (14.5 MB of Parquet) per 5 MeV Fe ion with `points: all` | ~3 100 rows per 5 MeV Fe ion and E_c = 1 keV |
| Changes the results | no | record mode: no; handoff mode: yes (see below) |

**Format.** The file name decides it: `*.csv` or `*.parquet`. A relative
name goes into the output directory (`output.directory`, CLI `-o`); an
absolute one is used as it is. CSV files have exactly one header line, the
column keys (no comment lines, so any CSV reader takes them as they are),
and numbers in the shortest form that reads back to the same double.
Parquet files (format 1, PLAIN, uncompressed, one row group per ~256 000
rows; the same subset of the format as AMCSET's other Parquet files) carry the units, labels
and descriptions of the columns and the run's facts (program, title, seed,
first ion, number of ions, axes, E_c, E_min, mode) in their key-value
metadata (`amcset.units`, `amcset.descriptions`, `amcset.E_c_eV`, ...).
Parquet is the better choice for trajectory files (they are 30 % smaller
and much faster to read).

**Axes.** Both files use the axes of AMCSET's trajectory output: **z is the
depth below the surface** (TRIM's x), x and y are lateral (TRIM's y and z).
Positions are in Å; directions are unit vectors in the same axes.

**Time.** `t_ps` is the time since the ion entered the target, in
picoseconds: the sum of flight length / speed along the atom's path, with
the speed (relativistic) at the energy at the start of each flight, and a
recoil starting at its parent's time at the collision that set it in
motion. It is an estimate for ordering and for handing events to MD: TRIM
has no time, the electronic loss along a flight is not resolved, and the
binary-collision picture ignores the time spent in a collision.

### Trajectory file columns

| Column | Meaning |
|---|---|
| `ion` | 1-based number of the incident particle (TRIM.DAT line for particle sources) |
| `point` | Index of the point in its ion's trajectory, from 0 (the ion's entry); `parent_point` refers to it |
| `track` | Path within the ion's trajectory: 0 is the ion, every recoil that moves has its own number |
| `parent_track`, `parent_point` | The track and point of the collision that set this track's atom in motion; -1 for the ion and for a recoil source's atoms. A track's first point has the parent point's position (exactly, except a recoil reflected at the surface, whose parent point is moved to where it is reflected) |
| `kind` | `birth` (first point of a track: the ion's entry or the site where an atom starts to move), `collision` (end of a flight, where the next collision happens), `stop` (the atom came to rest: last point of its track), `exit` (it left the target: last point) |
| `gen` | 0 = the ion, 1 = atoms the ion set in motion, 2 = atoms those set in motion, ... |
| `Z`, `mass_amu` | The moving atom |
| `x_A`, `y_A`, `z_A` | Position (z = depth) |
| `ux`, `uy`, `uz` | Direction after the collision at the point (0 0 0 for an atom at rest) |
| `E_eV` | Kinetic energy after the collision at the point |
| `E_before_eV`, `T_eV`, `dE_e_eV` | At collision points: the energy just before the collision, the energy given to the struck atom, and the electronic loss over the flight that led to the point (E_eV = E_before_eV − T_eV − dE_e_eV for recoils; for the ion the electronic loss comes first and E_before_eV is the energy after it). 0 at other points |
| `t_ps` | Time (see above) |

Which points exist follows the transport. In `trim2013` the ion has a
point at its entry and after every flight; a recoil has its birth point and
a point after every step of one atomic layer (TRIM's recoils move in
monolayer steps), with a collision at every step except the last. In
`refined` and `crystal` every moving atom (ion or recoil) has its birth
point, a point at every collision (the deflection point; a collision in
which several atoms are struck together is one point, with the summed
transfer in `T_eV`), and a last point where it stops or leaves; a recoil's
birth point is the struck atom's site, within an impact parameter of its
parent's collision point, and E_eV = E_before_eV − T_eV − dE_e_eV holds for
every collision point, the ion's included. Atoms TRIM sets in
motion but never moves (a struck atom below Ed deeper than a few atomic
layers, or the primary of a sub-threshold ion collision that TRIM restarts
with energy 0, see below) have no points. `output.trajectory_file.points:
"births"` keeps only the birth points (one row per moving atom: where and
with which energy and direction every atom starts); `min_energy` keeps
only points with at least that energy.

### Other transports and crystals

The `refined` and `crystal` transports write both files with the same
columns and rules (tracks, parents, times, events, record and handoff
modes); recording changes no tally in any transport (tested), and the
handoff mode conserves the energy as in `trim2013` (tested for each). In
the handoff mode the struck atom MD takes over is not booked by the BCA at
all (no displacement, vacancy or binding energy for it), as in `trim2013`.
The crystal transport's ion starts a few Å in front of the surface (outside
the target; its time starts there).

With the crystal transport and a crystal layer the files have two more
columns, and the Parquet files carry the crystals' orientation in their
metadata:

| Column | Meaning |
|---|---|
| `lattice_dx_A`, `lattice_dy_A` | The lateral offset of the crystal lattice for this ion: its entry point on the beam spot (`ion.spot`). Positions in the files (as in every result) are relative to the ion's entry point; a point (x, y, z) is at (x + `lattice_dx_A`, y + `lattice_dy_A`, z − top of the layer) in the frame of the layer's lattice, whose origin is a lattice point. Constant within an ion. |

| Metadata key (Parquet) | Value |
|---|---|
| `crystal.layerN.name`, `.followed`, `.top_A` | The cell, whether the crystal transport follows it, the depth of the layer's top |
| `crystal.layerN.lattice_vectors_A` | a, b, c in the crystal frame (9 numbers, row by row) |
| `crystal.layerN.rotation_crystal_to_file_axes` | R (9 numbers, row-major): a vector of the crystal frame v is R v in the files' axes (x, y lateral, z depth) |
| `crystal.layerN.euler_bunge_deg`, `.depth_direction`, `.beam_direction`, `.tilt_twist_deg`, `.temperature_K` | As in the "crystals" table |

CSV files have no metadata: the same information is in summary.json
(`crystals`) and in the `crystals` table (below), for every run that has a
crystal.

### MD event file

The events are what an MD code takes over from the binary-collision
calculation. For each transition energy E_c (`E_c`, a list in increasing
order: `"1 keV"` or `["0.5 keV", "1 keV", "2 keV"]`):

1. A **birth** event: an atom set in motion with an energy transfer
   E_min ≤ T < E_c. Its energy is T, its position the collision site, its
   direction the one it starts with.
2. A **crossing** event: a moving atom (the ion or a recoil) whose energy
   was at least E_c and falls below E_c at a collision (or when it loses
   its lattice binding energy Eb as it is displaced). Its energy is the
   energy it has after that collision, its direction the one it leaves
   with.
3. An ion that enters the target below E_c is a birth event at its entry.
4. Nothing that descends from an event's atom is an event again for the
   same E_c: the MD calculation is expected to follow that atom and all it
   sets in motion.
5. Atoms set in motion with T < E_min are not events; they are counted
   (`below_E_min_per_ion` in summary.json) and stay with the BCA, which
   treats them as TRIM does. `E_min: "Ed"` (the default) uses the struck
   element's displacement energy.
6. Transfers that TRIM treats as sub-threshold in the ion's collisions
   (T ≤ Ed) are never events: TRIM books them as the ion's phonons. With
   TRIM's `stale_primary_energy` bug on (the default) TRIM still runs a
   cascade from such a collision with a stale energy; that cascade
   produces no events (it is counted below E_min). Run with
   `--fix-trim-bugs` for MD coupling (see "Which TRIM bugs" below).
7. Atoms that leave the target (backscattered, transmitted, sputtered)
   before falling below E_c are not events.

| Column | Meaning |
|---|---|
| `ion_id` | 1-based number of the incident particle |
| `event_id` | 1-based within the ion and E_c, in time order (the file is sorted by ion, then E_c, then time) |
| `E_c_eV` | The transition energy the event belongs to (one file holds every E_c) |
| `kind` | `birth` or `crossing` (rules 1-3 above) |
| `t_ps` | Time since the ion entered |
| `x_A`, `y_A`, `z_A` | Position (z = depth) |
| `E_eV` | Kinetic energy handed to MD |
| `ux`, `uy`, `uz` | Unit direction |
| `gen`, `Z`, `mass_amu` | Generation (0 = the ion), atomic number, mass |
| `track`, `parent_track` | The atom's track and its parent's, as in the trajectory file (-1: no parent) |
| `parent_event` | Always -1: by rule 4 no event descends from another event of the same E_c. (Kept because coupling drivers expect the column.) |

**Record and handoff.** `mode: "record"` (the default) only writes the
events: the transport, the random numbers and every TRIM table stay exactly
as without them (tested at 1 and 4 threads), so one run gives the events of
several E_c and the TRIM results together. `mode: "handoff"` (one E_c)
also stops following every event's atom, so that nothing MD will simulate
is counted twice and the run is faster (5 MeV Fe → Fe at E_c = 1 keV, 20
ions on one thread: 0.6 s recorded, 0.1 s handed off). The TRIM tables of a handoff run cover the BCA part
only (fewer vacancies, phonons, ...), and are not comparable with TRIM.

**Totals.** `summary.json` gets an `md_events` block (mode, file, E_min and
per E_c: events, birth and crossing events and their energy per ion, the
energy as a percentage of the ion energy, atoms below E_min and their
energy, ions that became events), the CLI prints one line per E_c, and the
dataset has a table `md_events_summary` (GUI Tables tab; CSV, JSON and
Parquet outputs). In handoff mode the energy deposited by the BCA
(ionization, phonons, binding energy), the energy leaving the target and
the energy handed to MD add up to the ion energy: 80.355 % + 19.645 % =
100.000 % for 20 ions of 5 MeV Fe → Fe at E_c = 1 keV (the command below,
handoff mode) (`TrimMdEvents.HandoffConservesTheEnergy` checks it to 10⁻³).

**Which TRIM bugs.** For MD coupling, fix TRIM's bugs: `--fix-trim-bugs`
(`physics.trim_bugs: "fixed"`). **One caveat:** in a scenario an energy of 0
means "TRIM's default", and only TRIM's `zero_energy_defaults_skip_first_atom`
bug lets the first target atom keep Ed, Eb or Es = 0. With every bug fixed,
an Eb of 0 becomes 2 eV (the run says so: "a zero energy is replaced by
TRIM's default"). For Eb = 0 keep that one bug:
`--fix-trim-bugs --set physics.trim_bugs.zero_energy_defaults_skip_first_atom=true`
(the `--set` applies after `--fix-trim-bugs`), or in a scenario an object
of flags with that one `true` and the others `false`. With TRIM's behaviour,
`recoil_energy_created` creates energy in slow recoil collisions (Fe → Fe at
1 keV ends with about 105 % of the incident energy), `stale_primary_energy`
starts cascades with energies left over from earlier ones, and
`struck_atom_direction` gives atoms struck by recoils a random transverse
direction (docs/SRIM_MODE.md §7). The events are then what TRIM would
produce, which is rarely what MD should start from. Even with every bug
fixed, the trim2013 bookkeeping of a full run adds up to slightly more than
the ion energy (100.06 % for 60 ions of 20 keV Fe → Fe); the same run in
handoff mode adds up to 99.995 %.

**Recoil sources** (TRIM.DAT, cascade code 6): every atom of the list is a
root (`gen` 1, `parent_track` -1); an atom below E_c is one birth event, an
atom above it is followed and its pieces are events.

**Example** (the reference case of the hybrid method: 5 MeV Fe → Fe, Ed 40
eV, Eb 0, full cascades):

```bash
amcset run fe.json -n 20 --damage full --fix-trim-bugs \
       --set physics.trim_bugs.zero_energy_defaults_skip_first_atom=true \
       --set "ion.energy=5 MeV" --set "layers.0.thickness=4 um" --set "layers.0.elements.0.Eb=0 eV" \
       --md-events events.parquet --md-ec "0.5,1,2,5,10 keV" \
       --trajectory-file births.parquet --trajectory-points births -o out_md
```

At E_c = 1 keV this gives about 3 090 events and 0.98 MeV handed to MD per
ion (record mode, 20 ions; 3 075 events in handoff mode, where the
transport after an event differs).

```python
import pyarrow.parquet as pq, json
ev = pq.read_table("out_md/events.parquet").to_pandas()
one = ev[ev.E_c_eV == 1000]                                   # one transition energy
meta = pq.read_table("out_md/events.parquet").schema.metadata
print(json.loads(meta[b"amcset.units"]), meta[b"amcset.mode"])
```

**Limits.** The files list the ions of one `run()`: a GUI *Continue* writes
the added ions to a new file of the same name (the earlier one is replaced);
a run stopped early lists the ions that finished (in order, except blocks
that finished after the stop). A trajectory file of every point is large:
about 11 MB of CSV or 8 MB of Parquet per 1 MeV Fe ion, 14.5 MB of Parquet
per 5 MeV Fe ion; use `points: births` or `min_energy` when the full paths
are not needed. The refined transport
(`physics.transport: refined`) and the native engine do not support these
files.
