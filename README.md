<p align="center"><img src="resources/icons/amcset-128.png" width="96" alt="AMCSET logo"></p>

# AMCSET

AMCSET simulates ions (and electrons) moving through layered materials and
the radiation damage they cause: ranges, vacancies, replacements,
interstitials, displacements per atom, energy deposition, backscattering,
transmission and sputtering. Its SRIM mode reproduces SRIM-2013's TRIM
statistically and runs on every processor core; beyond TRIM, layers can be
crystals in which ions channel. It has a desktop program (`amcset-gui`) and a
command-line program (`amcset`) for scripts and computing clusters.

Developed at Texas A&M University. License: GPL-3.0-or-later.

## Installing

Packages for Windows 10/11 (MSI, ZIP, MSIX) and Linux (AppImage, `.deb`,
`.rpm`, a command-line tarball) are attached to each
[release on GitHub](https://github.com/Shao-RMSGAL/AMCSET/releases):
[docs/INSTALL.md](docs/INSTALL.md). The Windows files of 1.0 are not signed
yet, so Windows asks before it runs them and Smart App Control blocks them:
INSTALL.md, "Running the unsigned Windows release", says what to do. To
build from source, see "Building" below.

## What it does

- **SRIM mode** (default) follows TRIM.exe step by step (reverse engineered, [docs/TRIM2013_ALGORITHM.md](docs/TRIM2013_ALGORITHM.md)) with TRIM's inputs (per-element displacement, lattice and surface binding energies; quick Kinchin-Pease or full-cascade damage), tallies and output files (RANGE.txt, VACANCY.txt, ...), and reads TRIM.IN directly. It uses the SRIM-2013 stopping tables of every ion in every element, which come with it (`data/srim2013`). Against 161 SRIM-2013 TRIM runs of 75 distinct inputs, 1031 of 1032 comparisons of ranges, damage, sputtering and energy partition agree within the statistics; the exception, replacements in Xe gas (+4.6 %), is still open ([docs/VALIDATION.md](docs/VALIDATION.md) §2).
- **Speed.** The validation cases took AMCSET 0.1.0 53 times less machine time than TRIM.exe; the transport has since become 1.35 times faster with identical results, and optional fast approximations (`--fast`) add a factor of about 2 with statistically the same results (SRIM_MODE.md §8).
- **TRIM's bugs** are reproduced by default and can be corrected one by one (SRIM_MODE.md §7). Relativistic Bethe-Bloch stopping extends the energy range to 100 GeV.
- **Damage.** dpa by depth at a fluence, by four damage models (NRT, Kinchin-Pease, arc-dpa, rpa) besides TRIM's vacancies; for surviving defects in metals arc-dpa agrees best with MD and experiments, while TRIM's full-cascade vacancies are about twice NRT (VALIDATION.md §5).
- **Crystals.** Layers can be single crystals (built-in structures, or CIF, POSCAR, LAMMPS and XYZ cell files; any orientation; thermal vibrations), and a crystal transport follows their lattice, so that ions channel. It reproduces a measured channeling dip and measured boron profiles in silicon in peak and width; ranges in a crystal exceed the amorphous ones through channeling, also for slow heavy ions at the usual 7° tilt, as measured profiles show ([docs/CLI.md](docs/CLI.md) "Crystals"). The depth of channeling tails depends on a division of the electronic stopping that the data bracket but do not fix.
- **Validated** against TRIM.exe, measured stopping powers (the IAEA database, NIST), measured ranges (SIMS, RBS, nuclear reactions), sputtering yields, ionization yields, channeling data and MD damage functions, with the deviations found and the corrections that can be switched on (an optional fix for SRIM-2013's too-low stopping of Ni ions, which brings measured Ni ranges at 3.5-9 MeV within 9 %): [docs/VALIDATION.md](docs/VALIDATION.md).
- **Coupling to molecular dynamics.** Every trajectory point and hybrid BCA-MD recoil events (CSV or Parquet), LAMMPS electronic-stopping tables ([docs/OUTPUT.md](docs/OUTPUT.md)).
- **Outputs.** TRIM's text files, CSV, JSON and Parquet, with the tables, columns, rows, units and number format chosen by the user.
- **Safety nets.** Inputs are checked against the program's physical limits, implausible set-ups are warned about, the target depth can be sized automatically, and a crash leaves a report naming the ion that caused it.
- The **native engine** (AMCSET's original code) also transports electrons. Several of its physics bugs were fixed in 0.1.0, but its results are not yet validated ([docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md)).

## Quick start

The GUI ([docs/GUI.md](docs/GUI.md)) sets up a calculation (ion, target
layers, crystals, damage model, outputs), runs it, and shows TRIM's plots,
tables and the particle paths in 3D while it runs. It opens and saves the
same scenario files as the CLI (JSON or TRIM.IN).

```bash
amcset template > fe.json            # a complete example scenario (Fe into Fe)
amcset run fe.json -n 5000 -o out    # TRIM-format results in out/
amcset trim TRIM.IN -o out           # run an existing TRIM input
amcset stopping --ion He --target Fe --from "10 keV" --to "10 GeV"
```

Every scenario key and option: [docs/CLI.md](docs/CLI.md); `man amcset` on
Linux.

## How it is built

- `common/`: the simulation core (both engines, stopping and scattering physics, the SRIM-mode engine in `trim/`, crystals, element and isotope data, scenario files, result tables and their writers). Physical quantities are Boost.Units types.
- `cli/`: the `amcset` program. `gui/`: the Qt 6 Widgets program (live plots and statistics, GPU path views through QRhi).
- `packaging/`: installation and packages; `tools/`: the TRIM comparison, SRIM data extraction, crystal studies, build and release scripts.

C++20, CMake, Conan 2 for the C++ dependencies (Boost, abseil, GoogleTest,
Qt 6.7 or later). Repository map and rules: [AGENTS.md](AGENTS.md).

## Building

Details for every platform: [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md).

**Windows** (Visual Studio with MSVC, CMake, Python, Conan 2; `conan profile detect` once):

```bat
tools\windows\build_and_test.bat           :: configure, build, run all tests
tools\windows\build_and_test.bat --conan   :: (re)install the Conan dependencies first
```

It installs the dependencies with Conan if needed (building Qt from source
takes hours the first time), builds the GUI, the CLI and the tests in
`build\Release`, runs `ctest` and writes `build\build_and_test.log`. Set
`AMCSET_NO_PAUSE=1` to run it from scripts without the final pause.

**Linux**, the core, CLI and tests without Qt (GCC 13 or later, CMake,
Ninja, Boost 1.83 or later, abseil, GoogleTest):

```bash
cmake -S . -B build-noqt -G Ninja -DCMAKE_BUILD_TYPE=Release -DAMCSET_BUILD_GUI=OFF
cmake --build build-noqt -j && ctest --test-dir build-noqt
```

With the GUI, either through Conan (`pip install -r requirements.txt`,
`conan profile detect --force`, then `./qbuild.sh -c -b -t`; `-r` starts the
GUI; the first run builds Qt) or against a distribution's Qt 6.7 or later
(Debian 13, Ubuntu 25.04, Fedora 42: install its Qt development packages
with the private headers, then configure as above without
`-DAMCSET_BUILD_GUI=OFF`).

**Docker**: `./run-dev.sh` builds an Ubuntu 24.04 development image
(pinned GCC, CMake and Conan), mounts the repository, keeps the Conan cache
in a volume and forwards the display (Wayland or X11, and `/dev/dri`); inside
it, `./qbuild.sh -b -t -r` as above. `./run-dev.sh --rebuild` rebuilds the
image.

**Packages**: `cpack` in a build folder, and the release workflow;
[docs/RELEASING.md](docs/RELEASING.md).

## Documentation

Start with [docs/README.md](docs/README.md) (index of all documents).

- [docs/INSTALL.md](docs/INSTALL.md): installing the packages
- [docs/CLI.md](docs/CLI.md): command line and scenario files
- [docs/GUI.md](docs/GUI.md): the desktop program
- [docs/OUTPUT.md](docs/OUTPUT.md): output tables, units and formats
- [docs/SRIM_MODE.md](docs/SRIM_MODE.md): the TRIM-compatible engine and performance; [docs/VALIDATION.md](docs/VALIDATION.md): every validation
- [docs/TRIM2013_ALGORITHM.md](docs/TRIM2013_ALGORITHM.md): TRIM-2013's algorithm as recovered from TRIM.exe
- [docs/GLOSSARY.md](docs/GLOSSARY.md): terms and abbreviations
- [docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md): open bugs and limitations; [docs/ROADMAP.md](docs/ROADMAP.md): improvement backlog
- [CHANGELOG.md](CHANGELOG.md): changes by version
- For contributors: [CONTRIBUTING.md](CONTRIBUTING.md), [AGENTS.md](AGENTS.md) (rules, repository map, conventions; also read by coding agents, with [CLAUDE.md](CLAUDE.md) for the maintainers' workstation), [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md), [docs/HOWTO.md](docs/HOWTO.md), [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md), [docs/TESTING.md](docs/TESTING.md)
- For maintainers: [docs/RELEASING.md](docs/RELEASING.md) (packages, signing), [docs/PLATFORM.md](docs/PLATFORM.md) (Windows and Linux requirements)

## License

AMCSET is distributed under the GNU General Public License, version 3 or (at
your option) any later version; see [LICENSE](LICENSE). The bundled SRIM-2013
stopping tables (`data/srim2013`) are derived from SRIM (J. F. Ziegler,
www.srim.org) and are free for non-commercial use only
([data/srim2013/README.md](data/srim2013/README.md)).