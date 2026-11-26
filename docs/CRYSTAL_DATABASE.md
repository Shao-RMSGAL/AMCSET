# The crystal database

AMCSET 1.0.1 comes with a database of crystals for the crystal transport:
elements and compounds in AMCSET's built-in structures (bcc, fcc, hcp,
diamond, zinc blende, wurtzite, rock salt, CsCl, fluorite, perovskite, sc),
with the lattice constants at room temperature, the Debye temperature where
a source gives one, and the source of every value. It is compiled into the
program (`common/src/crystal/crystal_database.cpp`, `crystal_database.h`).

## Using it

| Where | How |
|---|---|
| GUI | Crystal tab, *Database...*: search, see the sources of the selected row, *Give the layer this crystal*. An amorphous layer whose composition is a database crystal shows *Give it this crystal*. |
| GUI | Crystal tab, crystal transport settings: *Amorphous layers of a database crystal's composition become that crystal* (`physics.crystal.auto_crystal`). |
| CLI | `amcset crystals [WORDS ...]` lists the rows that match every word (a name, element symbols or names, a structure, a word of the note), with density and sources: `amcset crystals nitride`, `amcset crystals Ga As`, `amcset crystals hcp`. |
| Scenario | `layers.N.crystal: {"material": "GaAs", "orientation": {...}, "temperature": "300 K"}`: the database crystal (its cell and, unless `vibrations` is given, its Debye temperature). A name with several forms gives the preferred one (`"ZnS"`: sphalerite); words choose another (`"ZnS wurtzite"`). |
| Scenario | `physics.crystal.auto_crystal: true` (default false): with the crystal transport, every layer without a crystal whose elements and atomic fractions (within 0.01) are those of a database crystal becomes that crystal when the calculation starts ((0 0 1) surface, 293 K); the engine's warnings say which, and the layer's density becomes the crystal's. Gas layers and alloys are left amorphous. |

When a composition has several forms, the one most common at room
temperature is chosen (`preferred`): ZnS → sphalerite (zinc blende), CdS →
greenockite (wurtzite), ZnO → zincite (wurtzite; the rock-salt form is the
high-pressure phase).

## Sources

Gemmell's Tables IV and V (the 43 crystals that were built in before,
`crystal/materials.h`) come first and are kept for the materials they list,
so a crystal chosen from *Common crystals*, the database or a scenario is
the same. The other rows:

| Key | Source | Rows |
|---|---|---|
| `gemmell1974` | D. S. Gemmell, Rev. Mod. Phys. 46, 129 (1974), Appendix B, Tables IV and V | 19 elements with Debye temperatures, 24 compounds |
| `wolfram` | Wolfram Research, Mathematica ElementData, as listed at periodictable.com, "Lattice Constants" (retrieved 2026-11-26). Structures as Wikipedia, "Periodic table (crystal structure)", gives them at room temperature | 35 elements (bcc, fcc, hcp, sc) |
| `kittel2004` | C. Kittel, *Introduction to Solid State Physics*, 8th ed. (Wiley, 2004), Debye temperatures as tabulated in Wikipedia, "Debye model" | Debye temperatures of Be, Ti, Zn, Cd, Cs |
| `ioffe` | Ioffe Institute, "Semiconductor NSM" | AlAs, AlP, InP, AlN, GaN, InN |
| `aps` | Argonne National Laboratory, Advanced Photon Source, "Lattice Constants" calculator | CdS, CdSe, CdTe, ZnS, BN, BP, ZnO (rock salt), PbTe |
| `morkoc2009` | H. Morkoç, Ü. Özgür, *Zinc Oxide* (Wiley-VCH, 2009) | ZnO (wurtzite) |
| `toth1967` | L. E. Toth, *Transition Metal Carbides and Nitrides* (Academic Press, 1967) | TiC, ZrC, HfC, VC, NbC, TaC (measured on slightly substoichiometric carbides; the note says which) |
| `saha2010` | B. Saha et al., J. Appl. Phys. 107, 033715 (2010) | ScN |
| `landolt` | J. B. Goodenough, J. M. Longo, Landolt-Börnstein III/4a (perovskites) | SrTiO3, KTaO3, SrVO3, CaVO3 |
| `wikipedia` | Wikipedia, "Lattice constant", rows of its table without a reference | LiCl, LiBr, LiI, NaBr, KF, RbF, RbCl, RbI, CsF, CsCl, CsBr, TiN, ZrN, VN, CrN, NbN |

The compounds' values (all but Gemmell's) were transcribed from the table
"Lattice constants for various materials at 300 K" of Wikipedia, "Lattice
constant" (revision of 13 October 2026), with the reference each row
gives. `CrystalDatabase.EveryRowBuildsWithItsSources` checks every row
against its structure, and every element's density against AMCSET's element
table (within 4 %), which catches a mistyped lattice constant or structure.

## What is not in it

Only the built-in structures are possible, so graphite, h-BN, alpha-Mn,
beta-Sn, alpha-U, Ga, In, As, Sb, Bi, the lanthanides of the double-hexagonal
structure (La, Pr, Nd, Am, ...), Rb and Ac (the only source's lattice
constants disagree with their room-temperature densities: Rb's is a
low-temperature value, Ac's gives 8.3 g/cm3 instead of about 10), HfN
(the source repeats NbN's lattice constant for it, 9 % too dense), corundum (Al2O3), quartz, rutile and the
orthorhombic or tetragonal perovskites are missing; SiC is missing because
the source gives a polytype (4H) the wurtzite structure cannot hold. Use a
cell file (CIF) for those (docs/CLI.md "Cell files"). Compounds have no
Debye temperature: their atoms stand still unless `vibrations` gives one
(the plausibility warnings say so). Wurtzites use the ideal internal
parameter u = 3/8.

## Adding a crystal

Add a row to `crystal_database()` in `common/src/crystal/crystal_database.cpp`
with its source (a `database_sources()` key; add the source if it is new,
here and in the table above), and run `common-test
--gtest_filter=CrystalDatabase*`. Do not add values without a source.
