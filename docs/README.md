# AMCSET documentation

Start with `AGENTS.md` in the repository root (rules, map, build, workflow).
Then read what your task needs:

| If you are going to… | Read |
|---|---|
| change anything | `AGENTS.md`, then the `CLAUDE.md` of the directory you work in |
| understand how the pieces fit | [ARCHITECTURE.md](ARCHITECTURE.md) |
| install AMCSET from a package | [INSTALL.md](INSTALL.md) |
| make a release, sign it | [RELEASING.md](RELEASING.md), [PLATFORM.md](PLATFORM.md) |
| build, run, validate, benchmark | [DEVELOPMENT.md](DEVELOPMENT.md) |
| add a setting, bug switch, table, plot, option, format | [HOWTO.md](HOWTO.md) |
| touch SRIM-mode physics | [SRIM_MODE.md](SRIM_MODE.md) and [TRIM2013_ALGORITHM.md](TRIM2013_ALGORITHM.md) |
| write or fix tests | [TESTING.md](TESTING.md) |
| work on the command line or scenario files | [CLI.md](CLI.md) |
| work on result files and tables | [OUTPUT.md](OUTPUT.md) |
| work on the GUI | [GUI.md](GUI.md) |
| pick something to improve | [ROADMAP.md](ROADMAP.md), [KNOWN_ISSUES.md](KNOWN_ISSUES.md) |
| look up a term (Ed, KP, PKA, ε, monolayer, z, …) | [GLOSSARY.md](GLOSSARY.md) |

## All documents

| Document | Contents | Audience |
|---|---|---|
| [ARCHITECTURE.md](ARCHITECTURE.md) | modules, scenario pipeline, SRIM-mode engine internals, native engine, GUI threading, data files, performance, extension points | developers |
| [DEVELOPMENT.md](DEVELOPMENT.md) | build matrix, Linux/Windows/Docker builds, Windows hurdles, running, debugging, validation against TRIM, benchmark, CI, git hygiene | developers |
| [HOWTO.md](HOWTO.md) | recipes for the common changes | developers |
| [TESTING.md](TESTING.md) | suites, what they cover, optional SRIM-data tests, TRIM comparison, GUI smoke test, determinism harness | developers |
| [SRIM_MODE.md](SRIM_MODE.md) | SRIM-mode physics, parameters, stopping tables, tallies, TRIM's bugs, performance | users, developers |
| [GPU.md](GPU.md) | the GPU backends (OpenCL, HIP, CUDA, Kokkos): using them from the CLI and the GUI, how a run works (diagrams, animations), precision, the comparison with the CPU, speed, profile, building, changing the kernel; figures in `gpu/figures/` | users, developers |
| [VALIDATION.md](VALIDATION.md) | Every comparison with references: 161 TRIM runs and published SRIM calculations, the SRIM-2013 stopping tables against the IAEA database and NIST, measured ranges (Ni ions, 350 ranges in Si and light targets), sputtering yields, damage models against MD, the crystal transport against measurements, the ionization of slow recoils; the deviations and the corrections that can be switched on | users, developers |
| [validation/](validation/README.md) | The validation as a review with 35 figures, in nine chapters: TRIM.exe, stopping, the Ni-ion problem, ranges, sputtering, damage, crystals, ionization, conclusions and recommendations | users, developers |
| [TRIM2013_ALGORITHM.md](TRIM2013_ALGORITHM.md) | TRIM.exe's algorithm as recovered by disassembly, with routine addresses | developers |
| [INSTALL.md](INSTALL.md) | the packages for Windows and Linux, checking a download, where the data are found | users |
| [CLI.md](CLI.md) | the `amcset` program and the full scenario (JSON) schema | users |
| [OUTPUT.md](OUTPUT.md) | output formats, tables, columns, OutputSpec | users |
| [GUI.md](GUI.md) | the desktop application and its code layout | users, developers |
| [KNOWN_ISSUES.md](KNOWN_ISSUES.md) | open bugs and limitations (native engine, SRIM mode, GUI, build) | developers |
| [RELEASING.md](RELEASING.md) | the packages and how they are built, the release workflow, code signing (options, setting up Azure Artifact Signing), catalogues | maintainers |
| [PLATFORM.md](PLATFORM.md) | Windows' desktop-app requirements and Linux standards, and how AMCSET 1.0 meets each | maintainers |
| [ROADMAP.md](ROADMAP.md) | prioritised improvement backlog | developers |
| [CRYSTAL_DATABASE.md](CRYSTAL_DATABASE.md) | the crystal database: what it holds, where every value comes from, how to search it (GUI, `amcset crystals`, `"material"`) and the automatic crystals | users, developers |
| [CRYSTAL_DESIGN.md](CRYSTAL_DESIGN.md) | crystalline targets: the design, what is built (unit cells, orientation, cell files, the scenario's `crystal`, the crystal transport, common crystals, the GUI tab, a first comparison with measurements; user guide in CLI.md "Crystals" and GUI.md) and what is still a proposal (polycrystals, the dose-dependent lattice, a fitted stopping division) | developers |
| [GLOSSARY.md](GLOSSARY.md) | terms and their names in the code | everyone |

Elsewhere: `README.md` (overview), `CHANGELOG.md` (versions),
`CONTRIBUTING.md`, `data/srim2013/README.md` (bundled tables and their
terms), `tools/srim/README.md` (SRIM data and TRIM reference runs).
