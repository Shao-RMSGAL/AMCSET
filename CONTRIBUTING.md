# Contributing to AMCSET

Thanks for helping. Before changing code, read [AGENTS.md](AGENTS.md) (the
layout, conventions, units and the rules for changes) and
[docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md).

## Reporting a problem

Open an issue on GitHub with:

- the version (`amcset version`, or Help > About in the GUI) and the system
  (Windows or Linux, which package: MSI, ZIP, AppImage, ...);
- the scenario file or TRIM.IN, and the command or the steps in the GUI;
- what you expected and what happened. After a crash the program prints the
  path of a crash report: attach it. The CLI writes them to
  `%TEMP%\amcset-crashes` or `/tmp/amcset-crashes` (with `--log FILE`, next
  to the log), the GUI to its log folder (docs/GUI.md).

Questions about a comparison with TRIM: say which TRIM version and settings
were used (docs/SRIM_MODE.md lists what SRIM mode reproduces).

## Workflow

1. Work on a branch off `develop`; `main` holds releases.
2. Build and run the tests: Windows `tools\windows\build_and_test.bat`, Linux
   `cmake --build <dir> && ctest --test-dir <dir>`
   ([docs/DEVELOPMENT.md](docs/DEVELOPMENT.md)).
3. A change that must **not** change results (refactoring, speed):
   - SRIM mode (`trim2013`, `refined`): `tools/perf/bench.py run` with the
     old and the new executable, then `bench.py compare`: every output file
     identical. Compare builds of one compiler (MSVC and GCC differ in the
     last bits of `exp`, `log` and `pow`).
   - the native engine: the determinism harness before and after
     ([docs/TESTING.md](docs/TESTING.md)).
4. A change that **does** change results (a physics fix): one thing at a
   time, with the evidence (TRIM.exe's routine, a published value, or a
   statistically significant comparison) in the code comment and the pull
   request, a test, and a [CHANGELOG.md](CHANGELOG.md) entry.
5. SRIM mode (`common/src/trim/`, the stopping tables, `common/src/physics/`):
   run the TRIM comparison (`tools\windows\compare_with_trim.bat`,
   [docs/VALIDATION.md](docs/VALIDATION.md) §2) and confirm the agreement still
   holds or explain the change. Behaviour that TRIM.exe gets wrong is
   reproduced by default and corrected only through a `physics.trim_bugs`
   switch, never silently.
6. The crystal transport is AMCSET's own physics: keep the trends of
   `tools/crystal/channeling_checks.py` (a random direction behaves like the
   amorphous target) and leave the other transports untouched.
7. Add or update tests for what you changed, and the documents that
   describe it ([AGENTS.md](AGENTS.md) §7 lists which document owns what).

## Requirements for every change

- Every `.h`, `.cpp` and `.py` file carries the exact GPL header of the
  existing files (`python3 scripts/check_headers.py` runs in CI; a shebang
  goes on line 1 above it).
- Warning-free with the project's warning flags, on MSVC and GCC.
- Physical quantities are Boost.Units types; conversions go through units,
  not bare factors ([AGENTS.md](AGENTS.md) §5).
- A new scenario setting goes into `config::RunConfig` (JSON read and
  write, validation, `--set`), the GUI panel that edits it, and
  [docs/CLI.md](docs/CLI.md); see [docs/HOWTO.md](docs/HOWTO.md).
- Never edit `common/generated/*.h` by hand.
- No TRIM, MARLOWE or other third-party simulation code, disassembly or
  binaries in the repository: AMCSET re-implements documented behaviour.
  The only SRIM-derived data are `data/srim2013/srim2013_stopping.amcset`
  (non-commercial terms, see its README).
- Line endings follow `.gitattributes` (LF; CRLF for `.bat` and `.ps1`).
  Do not commit build outputs or the synchronisation artefacts of a
  shared folder (`*.conflict*`); `.gitignore` lists both.

## Releases

Maintainers: [docs/RELEASING.md](docs/RELEASING.md) (packages, the release
workflow, code signing) and [docs/PLATFORM.md](docs/PLATFORM.md) (what the
packages must meet on Windows and Linux).

## License

By contributing you agree that your contributions are licensed under the
GNU General Public License, version 3 or (at your option) any later version.
