# Contributing to AMCSET

Thanks for helping. Before changing code, read [AGENTS.md](AGENTS.md) — it
describes the layout, conventions, units and the rules for changes — and
[docs/KNOWN_ISSUES.md](docs/KNOWN_ISSUES.md).

## Workflow

1. Work on a branch off `develop`.
2. Build and run the tests (Linux: `cmake --build <dir> && ctest --test-dir
   <dir>`; Windows: `tools\windows\build_and_test.bat`, see
   [docs/DEVELOPMENT.md](docs/DEVELOPMENT.md)).
3. If you touched `common/` and did **not** intend to change results, run the
   determinism harness before and after and confirm the outputs are
   byte-identical ([docs/TESTING.md](docs/TESTING.md)).
4. If you **did** intend to change results (a physics fix), change one thing
   at a time, explain the evidence in the commit/PR, update
   `docs/KNOWN_ISSUES.md`, and regenerate the harness reference.
5. If you touched SRIM mode (`common/src/trim/`, the stopping tables, the
   physics headers), run the TRIM comparison
   (`tools\windows\compare_with_trim.bat`) and confirm 1031/1032 within 3σ
   still holds, or explain the change ([docs/SRIM_MODE.md](docs/SRIM_MODE.md)
   §6). Behaviour that TRIM.exe gets wrong is reproduced by default and
   fixable through a `physics.trim_bugs` switch, never silently corrected.
6. Add or update tests for the behaviour you changed, and update the docs
   that describe it (AGENTS.md §7 lists which doc owns what).

## Requirements for every change

- Every `.h`, `.cpp` and `.py` file carries the exact GPL header used by the
  existing files (`scripts/check_headers.py` runs in CI; a shebang, if any,
  goes on line 1 above it).
- Code compiles warning-free with the project's warning flags.
- Physical quantities use Boost.Units types; conversions go through units,
  not bare factors ([AGENTS.md](AGENTS.md) §5).
- A new scenario setting goes into `config::RunConfig` (JSON read/write,
  validation, `--set`), the GUI panel that edits it, and
  [docs/CLI.md](docs/CLI.md); see [docs/HOWTO.md](docs/HOWTO.md) §1.
- Never edit `common/generated/isotope_database.h` by hand.

## License

By contributing you agree that your contributions are licensed under the
GNU General Public License, version 3 or (at your option) any later version.
