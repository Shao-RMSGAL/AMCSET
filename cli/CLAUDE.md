# cli/ — the `amcset` command-line program

User documentation and the full scenario schema: `docs/CLI.md`.

## Structure

| File | Contents |
|---|---|
| `include/cli.h`, `src/cli.cpp` (library `amcset-cli-lib`) | `Options` (parsed command line), `parse_arguments` (throws `UsageError`, exit status 2), `apply_options` (options → `config::RunConfig`, options win over the file), `looks_like_trim_in`, `parse_composition`, `print_stopping_table` (the `stopping` command), `print_cell` (the `cell` command: unit cells from files or built-in structures, conversion, the scenario's `crystal` object, beam directions) |
| `include/interrupts.h`, `src/interrupts.cpp` (in `amcset-cli-lib`) | Signal and console-event handling during a run (`InterruptGuard`, `exit_after_interrupt`) |
| `src/main.cpp` (executable `amcset`) | Commands `run`, `trim`, `convert`, `stopping`, `elements`, `cell`, `crystals`, `grains`, `devices` (GPU backends), `template`, `validate`, `version`, `help`, the hidden `diag-crash` (crashes on purpose to check crash reports); the help text; progress bar; summary printing (`--json`); `diag::initialize` (`--log`, `-v`) |

The library is separate so that `test/cli/cli_test.cpp` can test parsing and
option mapping without spawning processes; the end-to-end tests
(`ctest -R cli`) run the executable (declared in `test/CMakeLists.txt`).

## Rules

- **The CLI is a thin front end.** It maps options onto a `RunConfig` and
  calls `config::make_engine`. Anything that changes a calculation must be a
  scenario key first (`common/CLAUDE.md`); a CLI option is a shortcut for
  it (`--set KEY=VALUE` already reaches every key).
- **New option:** add the field to `Options`, parse it in
  `parse_arguments` (unknown options are usage errors), map it in
  `apply_options`, document it in the help text in `main.cpp` and in
  `docs/CLI.md`, test it in `cli_test.cpp`.
- **TRIM bug flags** (`--fix-trim-bug NAME`, `--fix-trim-bugs`) iterate
  `trim::TrimBugs::members()`; a new bug switch needs no CLI change.
- **Exit codes:** 0 success, 1 run/config error, 2 usage error, 3 after a
  crash (report written by `diagnostics.h`). An interrupted run (SIGINT,
  SIGTERM, SIGHUP; Windows console events: `interrupts.h`) writes its
  partial results, then re-raises the signal (POSIX: status 128 + signal
  number) or exits with 130 (Windows). Signal handlers only set lock-free
  atomics. Messages to stderr; results/JSON to stdout.
- Keep behaviour stable for scripts: `tools/trim_compare/batch.py` calls
  `amcset trim TRIM.IN -n N -s SEED -q --formats trim -o DIR` and
  `compare.py` reads the TRIM-format files and `summary.json` it writes.
