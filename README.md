# AMCSET

AMCSET is a C++ application for simulating ion and electron bombardment of
layered materials: ranges, damage (vacancies, replacements, interstitials),
energy deposition, backscattering, transmission and sputtering. It has a Qt
desktop interface and a full command-line program.

Two engines are available:

- **SRIM mode** (default) reproduces SRIM-2013 TRIM statistically, with
  TRIM's input parameters (per-element displacement, lattice and surface
  binding energies; quick Kinchin-Pease or full-cascade damage), TRIM's
  tallies and TRIM-format output files (RANGE.txt, VACANCY.txt, ...). It
  reads TRIM.IN files directly. It follows TRIM.exe step by step (reverse
  engineered, `docs/TRIM2013_ALGORITHM.md`) and uses the SRIM-2013 stopping
  tables of every ion in every element that come with it (`data/srim2013`).
  Against 161 SRIM-2013 TRIM runs of 75 distinct inputs, 1031 of 1032
  comparisons of ranges, damage, sputtering and energy partition agree
  within the statistics; the exception, replacements in Xe gas (+4.6 %), is
  still open ([docs/SRIM_MODE.md](docs/SRIM_MODE.md) §6). The validation
  cases took AMCSET 0.1.0 53 times less machine time than TRIM.exe (§8);
  the transport has since become 1.35 times faster with identical results,
  and optional fast approximations (`--fast`) add a factor of about 2 with
  statistically the same results (§8.1-8.2). The bugs
  of TRIM.exe are reproduced by default and can be corrected one by one (§7).
  Relativistic Bethe-Bloch stopping extends the energy range to 100 GeV.
  Inputs are checked against the program's physical limits, implausible
  set-ups are warned about, the target depth can be sized automatically,
  and crashes leave a report with the ion that caused them
  ([docs/CLI.md](docs/CLI.md)). For coupling to molecular dynamics it can
  write every trajectory point and hybrid BCA-MD recoil events (CSV or
  Parquet) and LAMMPS electronic-stopping tables
  ([docs/OUTPUT.md](docs/OUTPUT.md)). dpa by depth at a fluence comes out
  directly, by four damage models (NRT, Kinchin-Pease, arc-dpa, rpa) besides
  TRIM's vacancies; for surviving defects in metals arc-dpa agrees best with
  MD and experiments, while TRIM's full-cascade vacancies are about twice
  NRT ([docs/SRIM_MODE.md](docs/SRIM_MODE.md) §6.3). Beyond TRIM, layers can
  be single crystals (built-in structures or CIF, POSCAR, LAMMPS and XYZ
  cell files, any orientation, thermal vibrations) and a separate crystal
  transport follows their lattice, so that ions channel; the GUI's Crystal
  tab edits the unit cell in a 3D view. It reproduces a measured channeling
  dip and measured boron profiles in silicon in peak and width; the depth
  of channeling tails depends on a division of the electronic stopping that
  those data bracket but do not fix ([docs/CLI.md](docs/CLI.md) "Crystals").
- The **native** engine (the original AMCSET code) also transports
  electrons. Several of its physics bugs were fixed in 0.1.0, but its
  results are not yet validated ([docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md)).

```bash
amcset template > fe.json            # a complete example scenario (Fe into Fe)
amcset run fe.json -n 5000 -o out    # TRIM-format results in out/
amcset trim TRIM.IN -o out           # run an existing TRIM input
amcset stopping --ion He --target Fe --from "10 keV" --to "10 GeV"
```

See [docs/CLI.md](docs/CLI.md). The GUI ([docs/GUI.md](docs/GUI.md)) opens
and saves the same scenario files (JSON or TRIM.IN), shows TRIM's plots live
and draws the particle paths in 2D and 3D on the GPU. Results can be written
in TRIM's text format, CSV, JSON and Parquet, with the tables, columns, rows,
units and number format chosen by the user ([docs/OUTPUT.md](docs/OUTPUT.md)).

## How it works

- `common/` — the simulation core: both engines, stopping and scattering
  physics (`physics/`), the SRIM-mode engine (`trim/`), element and isotope
  data, scenario files. All physical quantities use Boost.Units.
- `cli/` — the `amcset` command-line program.
- `gui/` — the Qt 6 Widgets desktop application (setup, live plots and
  statistics, GPU path views through QRhi); see docs/GUI.md.
- `common/src/output/` — result tables with units, and the CSV, JSON and
  Parquet writers.
- `tools/` — TRIM comparison (`trim_compare`), SRIM data extraction and TRIM
  batch runs (`srim`), the Windows build script and a determinism harness.

The project uses C++20, CMake, and Conan to manage its C++ dependencies.
Relative paths (the electron parameter tables in `common/assets/`) are
resolved against the working directory, so launch the GUI from the
repository root (as `./qbuild.sh -r` does) when using the native engine
with electrons.

## Build and run on Linux

Ubuntu 24.04 is the environment used by the project's Docker image. Building
Qt and other dependencies from source can take a long time, especially on the
first build.

Install the build tools and Conan 2.x. You will need a C++20-capable compiler
(GCC 14 is used in the Docker image), CMake 3.23 or newer, Ninja, Python 3,
pipx, and the system development libraries required by Qt 6 (including the
X11/Wayland and OpenGL/Vulkan development libraries). On Ubuntu, the base
tools can be installed with:

```bash
sudo apt update
sudo apt install build-essential cmake ninja-build python3 python3-pip python3-venv pipx
pipx install "conan>=2,<3"
export PATH="$HOME/.local/bin:$PATH"
```

If your distribution does not provide a recent enough CMake, Ninja, compiler,
or Qt system dependency, use the Docker workflow below. Then, from the
repository root:

```bash
conan profile detect --force
./qbuild.sh -c -b
```

`-c` installs the dependencies through Conan and `-b` configures and builds
the Release application. The first dependency install may compile packages
from source and take substantial time. The detected Conan profile must match
your installed compiler; the build requests the GNU C++20 and C17 standards.

This builds the GUI, the `amcset` CLI and the tests (CMake options
`AMCSET_BUILD_GUI`, `AMCSET_BUILD_CLI`, `AMCSET_BUILD_TESTS`; the CLI and
tests need no Qt). Once built, use the helper script to run tests or launch
the GUI:

```bash
./qbuild.sh -t  # Run the tests
./qbuild.sh -r  # Launch the GUI
```

You can also perform all four steps in sequence with
`./qbuild.sh -c -b -t -r`. The build and test files are placed under
`build/Release/` (tests in `build/Release/test`).

## Build and run in Docker

The Docker image provides an Ubuntu 24.04 development environment with pinned
compiler and build-tool versions. The included `run-dev.sh` helper is designed
for Linux hosts: it builds the image if needed, mounts this repository into
the container, persists the Conan cache in a Docker volume, and forwards an
available Wayland or X11 display (and `/dev/dri` when available).

Requirements:

- Docker Engine and the Docker CLI
- A Linux desktop session for displaying the GUI
- For VS Code integration, the Dev Containers extension

From the repository root, start the development container:

```bash
./run-dev.sh
```

The first image build and Conan dependency setup can take a while. When the
script finishes, attach VS Code to `amcset-dev-container` using **Dev
Containers: Attach to Running Container**. In a terminal in the container,
build, test, and run as usual:

```bash
./qbuild.sh -b
./qbuild.sh -t
./qbuild.sh -r
```

`run-dev.sh --rebuild` forces the Docker image to be rebuilt. The workspace is
mounted from the host, so source changes and generated build files remain in
the repository. Docker keeps downloaded Conan packages in a named volume so
they can be reused on later runs. If the host has no usable display, the
container can still be used to build and test, but the GUI will not appear.

## Windows

Install Conan 2, CMake, Python and Visual Studio (MSVC, C++20) and run
`conan profile detect`. Then double-click `tools\windows\build_and_test.bat`
or run it from a command prompt in the repository root:

```bat
tools\windows\build_and_test.bat           :: configure, build, run all tests
tools\windows\build_and_test.bat --conan   :: (re)install the Conan dependencies first
```

It installs the dependencies with Conan only if they are not installed yet
(building Qt from source takes hours the first time), builds the GUI, the
CLI and the tests in `build\Release`, runs `ctest` and writes
`build\build_and_test.log`. Set `AMCSET_NO_PAUSE=1` to run it from scripts
without the final pause. The PowerShell helper `qbuild.ps1` (`-c`, `-b`,
`-t`, `-r`) also works, but the batch file is the tested path.

## Documentation

Start with [docs/README.md](docs/README.md) (index of all documents).

- [AGENTS.md](AGENTS.md) — rules, repository map, build, conventions (for
  contributors and coding agents; [CLAUDE.md](CLAUDE.md) adds the
  workstation specifics)
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — modules, data flow, threads
- [docs/HOWTO.md](docs/HOWTO.md) — recipes for common changes
- [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md) — building, debugging,
  validation runs, benchmarks
- [docs/SRIM_MODE.md](docs/SRIM_MODE.md) — the TRIM-compatible engine,
  validation and performance
- [docs/TRIM2013_ALGORITHM.md](docs/TRIM2013_ALGORITHM.md) — TRIM-2013's
  algorithm as recovered from TRIM.exe
- [docs/CLI.md](docs/CLI.md) — command line and scenario files
- [docs/GUI.md](docs/GUI.md) — the desktop application
- [docs/OUTPUT.md](docs/OUTPUT.md) — output tables, units and formats
- [docs/TESTING.md](docs/TESTING.md) — tests, TRIM comparison, determinism
  harness
- [docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md) — open bugs and physics issues
- [docs/ROADMAP.md](docs/ROADMAP.md) — improvement backlog
- [docs/GLOSSARY.md](docs/GLOSSARY.md) — terms and abbreviations
- [CHANGELOG.md](CHANGELOG.md) — changes by version
- [docs/CLEANUP_2026-11.md](docs/CLEANUP_2026-11.md) — log of the November
  2026 cleanup (historical)
- [CONTRIBUTING.md](CONTRIBUTING.md)

## License

AMCSET is distributed under the GNU General Public License, version 3 or (at
your option) any later version. See [LICENSE](LICENSE).
