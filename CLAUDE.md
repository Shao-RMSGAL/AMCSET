@AGENTS.md

# Claude Code notes for AMCSET

Everything in AGENTS.md (imported above) applies. This part covers how to
work in this repository with Claude Code: the machine, long-running jobs,
verification, and the things that have gone wrong before.

## The workstation

- Windows 11, AMD Ryzen 7 7800X3D (8 cores / 16 threads). Visual Studio 18
  (MSVC 14.5, "x64 Native Tools Command Prompt for VS"), Conan 2, CMake,
  Ninja, Python 3. Qt 6.11.1 and the other dependencies are built by Conan
  into `build\Release\generators` (do not delete `build\Release`: rebuilding
  Qt takes hours).
- Repository: `C:\Users\natha\Nextcloud\Code\C++\AMCSET` — inside a
  **Nextcloud-synchronised folder**. Sync can briefly lock or replace files
  (the scripts retry writes for that reason). Prefer writing large scratch
  output under `build\` and never commit sync artefacts (`*.conflict*`).
- SRIM-2013: `C:\Users\natha\Applications\SRIM-2013-Pro` (TRIM.exe, SR
  Module, Data). Its `_amcset\` subfolder is the TRIM reference workspace:
  `results\trim\` and `results\trim3\` hold every TRIM reference run
  (one folder per case with TRIM.IN and all TRIM outputs, `_done.txt` with
  TRIM's run time), `trim_cases\`, `trim_cases3\` the inputs. Do not modify
  these folders; they took days of TRIM time.
- Latest comparison and its statistics: `build\trim_validation\trim3\`
  (`analysis.txt`, `analysis.json`, `progress.txt` with AMCSET run times);
  earlier rounds in `trim3_v1`…`trim3_v3`.

## Building and testing from Claude Code

- Use the batch file non-interactively: `cmd /c "set AMCSET_NO_PAUSE=1 &&
  tools\windows\build_and_test.bat"`, then read the last line of
  `build\build_and_test.status` and, on failure, search
  `build\build_and_test.log` for `error` / `FAILED`. A no-change rebuild
  plus tests takes ~15 s; touching `trim_engine.cpp` recompiles in
  a few minutes (LTO).
- Running a single test binary: `call build\Release\generators\conanrun.bat`
  first (DLLs), then `build\Release\common-test.exe --gtest_filter=Trim*`.
  CTest runs tests from the repository root; do the same.
- `amcset-gui.exe` cannot be relinked while it is running (LNK1104); the
  rest of the build still succeeds. Ask the user to close the GUI.
- Windows may block a freshly built unsigned executable (Smart App
  Control / Device Guard); `build_and_test.bat` reports it in the log.
  `windeployqt` can be blocked the same way; the GUI then needs
  `conanrun.bat` on PATH.
- Linux alternative (no GUI) is in AGENTS.md §3; any Linux box with GCC ≥ 13,
  Boost and a recent abseil works for the core, CLI and tests.

## Long-running jobs

| Job | Command | Time | Needs |
|---|---|---|---|
| Full TRIM comparison | `tools\windows\compare_with_trim.bat` | ~1 h (all cores) | nothing else heavy running |
| One case vs TRIM | `amcset trim <case>\TRIM.IN -o out -s 101` then `python tools\trim_compare\compare.py <case> out` | seconds–minutes | |
| Benchmark vs TRIM.exe | `tools\windows\benchmark_vs_trim.bat` | 30–60 min | **idle machine**; opens TRIM windows |
| GUI smoke test | `tools\windows\gui_smoke_test.bat` | 1–2 min | a desktop session; writes `build\gui_smoke\<n>\` |
| New TRIM reference runs | `tools\srim\` scripts inside the SRIM folder | hours–days | the user's consent; TRIM opens many windows |

Run long jobs in the background and poll their progress files
(`progress.txt`, `*.status`, logs) instead of waiting on the process. Ask
before starting anything that takes the machine for more than a few minutes
or opens windows (TRIM, the GUI).

## Verifying claims

- Quote numbers only from files you produced or read in this session
  (analysis.txt, test output, benchmark CSV). The docs' numbers must match
  the latest analysis; if you rerun the comparison, update
  `docs/SRIM_MODE.md` §6.2 and the summary lines in `README.md` and
  `CHANGELOG.md` from `analysis.txt`.
- When you cite a TRIM.exe address in code or docs, it must already appear
  in `docs/TRIM2013_ALGORITHM.md` or be verified by disassembling TRIM.exe
  in this session. Do not guess addresses.
- A statistical "agreement" claim needs the comparison output (z values),
  not a single run.

## Things that went wrong before (and the fix)

- Comparing AMCSET with TRIM using a different stopping table than TRIM's:
  the bundled `data/srim2013` tables are the right default; the built-in
  stopping differs by up to 25 %.
- Treating a TRIM quirk as an AMCSET bug and "fixing" it in the trim2013
  transport: TRIM agreement broke. Bugs go behind `TrimBugs` switches.
- Single runs looked like systematic differences; they were run-to-run
  scatter. Use several seeds and `analyze.py`'s z statistics.
- Copying files to the workstation through a stale intermediate copy
  overwrote newer docs. After syncing files, compare checksums or grep for
  the new text.
- Long scripts with `pause`/`timeout` at the end hang non-interactive
  shells: use `AMCSET_NO_PAUSE=1` or call the `.ps1` directly.
