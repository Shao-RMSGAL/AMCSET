# How to make common changes

Step-by-step recipes. Each ends with the checks to run; the general rules
are in `AGENTS.md` (read it first). Paths are relative to the repository
root.

---

## 1. Add a scenario setting (SRIM mode or both engines)

1. **Field.** Add it to `trim::TrimConfig` (`common/include/trim/trim.h`)
   with a Doxygen comment saying what it does, its unit and its TRIM
   equivalent (or "AMCSET option, not in TRIM"). Use a quantity type for
   dimensional values; pick the default that reproduces TRIM.
   Native-only settings go to `config::NativeOptions` (`run_config.h`).
2. **Validate** in `TrimConfig::validate()` (`trim_engine.cpp`): throw
   `std::invalid_argument` naming the key.
3. **Parse** in `config::parse_config` (`common/src/run_config.cpp`): add
   the key to the section's `reject_unknown({...})` list and read it with
   the `Reader` helpers (`flag`, `count`, `num`, `str`, `quantity(...,
   parse_energy)`).
4. **Serialise** in `to_object` (same file) so that `to_json` round-trips.
   `apply_override` (`--set`) then works without further code.
5. **Use** it in the engine. In the trim2013 transport, an option that
   departs from TRIM must default to TRIM's behaviour.
6. **TRIM.IN** — only if TRIM.IN has the field: read it in
   `parse_trim_in` and write it in `write_trim_in` (`trim_in.cpp`).
7. **Front ends** (if users should set it directly): a CLI option
   (§6) and/or a GUI widget bound with `FormBinder` in the right setup
   panel (§4), with a tooltip.
8. **Docs:** the key in `docs/CLI.md` (scenario schema), the parameter table
   in `docs/SRIM_MODE.md` §3 if it is physics, `CHANGELOG.md`.
9. **Tests:** a round-trip/override test in `test/cli/cli_test.cpp`
   (`RunConfig` suite) and a behaviour test in `test/common/trim_test.cpp`.

Checks: `ctest`; `amcset template | amcset validate` style smoke test:
`amcset template -o t.json && amcset validate t.json`.

## 2. Add a TRIM bug switch

Use this when TRIM.exe does something wrong and you want AMCSET to be able
to do it right, without losing TRIM reproducibility.

1. **Evidence.** Identify the TRIM.exe routine (address or location label)
   and what it does; quantify the effect with AMCSET (fix on vs off, several
   seeds, the reference cases of `docs/SRIM_MODE.md` §7).
2. **Member.** Add `bool <name> = true;` to `trim::TrimBugs`
   (`common/include/trim/trim.h`) with a comment: what TRIM does, where,
   what the fix does. `true` always means "as TRIM".
3. **Register** it in `TrimBugs::members()` (`trim_engine.cpp`). That makes
   it known to the scenario parser (`physics.trim_bugs`), `to_json`, the
   CLI (`--fix-trim-bug NAME`, `--fix-trim-bugs`) and `set_all` /
   `all_reproduced`.
4. **Branch** in the engine on `cfg_.bugs.<name>` (or `config_.bugs` in the
   model construction). Keep the TRIM branch byte-for-byte as before so
   existing results do not change; keep the random-draw sequence of the
   TRIM branch identical.
5. **GUI:** add a `BugRow` in `CalculationPanel` (`gui/src/setup/
   calculation_panel.cpp`) with a label and a tooltip in the established
   form: "TRIM: …", "Fixed: …", "Effect of the fix: …" (measured numbers),
   followed by the shared caveat.
6. **Tests** (`test/common/trim_test.cpp`): update the expected count in
   `TrimBugs` (`members().size()`), and extend
   `TrimBugs.EachFixChangesWhatItShould` with the quantity the fix must
   change (and one it must not).
7. **Docs:** a row in the table of `docs/SRIM_MODE.md` §7, the bug list in
   `docs/TRIM2013_ALGORITHM.md`, `docs/CLI.md` if the names are listed,
   `CHANGELOG.md`.

Checks: `common-test --gtest_filter='TrimBugs*:Trim*'`; the default
(reproduce-all) results must be unchanged — exact tallies in `trim_test.cpp`
stay as they were.

## 3. Add a result (tally, table, column or scalar)

1. **Tally:** add the accumulator to `trim::TrimTallies` (`trim.h`); size it
   in `resize()`, merge it in `add()` (and `without_lists()` if it is a
   list).
2. **Accumulate** it in the `Transporter` (`trim_engine.cpp`). Depth
   tallies use `add(hist, x, value)` / `deposit(hist, a, b, value)` and the
   depth window; totals go in separate fields.
3. **Dataset:** in `trim::make_dataset` (`common/src/trim/trim_tables.cpp`)
   add a column to an existing table or a new `output::Table` with a unique
   `key`, `title`, `description` and typed columns
   (`Column::real(key, label, Dimension, unit, values)`, `Column::integer`,
   `Column::strings`), in TRIM's units. Scalars go to `data.scalars`. The
   GUI tables, CSV, JSON and Parquet get it automatically.
4. **TRIM file** — only if TRIM writes it: `trim_output.cpp` (`Writer`),
   matching TRIM's layout exactly (compare with a TRIM output file).
5. **GUI plot** (optional): §4.
6. **Docs:** the table list in `docs/OUTPUT.md`; `docs/SRIM_MODE.md` §5
   (tallies) if it is a new physical quantity.
7. **Tests:** `TrimTables` (table present, units, sizes) and an engine test
   for the value.

## 4. GUI: add a setup field or a plot

**Field.** Add the setting to `RunConfig` first (§1). In the panel
(`gui/src/setup/<x>_panel.cpp`) create the widget and bind it:

```cpp
binder_->text(edit, tr("Final energy"),
              [](const config::RunConfig &c) { return QString::fromStdString(config::format_energy(c.trim.final_energy)); },
              [](config::RunConfig &c, const QString &v) { c.trim.final_energy = config::parse_energy(v.toStdString(), "final energy"); });
```

`FormBinder` offers `text`, `check`, `combo`, `spin` and `custom`; edits go
through `Scenario::modify` (undoable; an exception rejects the edit and its
message is shown). Combo boxes bound to an enum by index must list the
items in enum order. Add a tooltip.

**Plot.** In `gui/src/views/plots_panel.cpp` add a `PlotDef` with an id, a
title, the dataset table key, the x column and the y columns (`*`
wildcards such as `recoils_*` expand per element), plus the axis options.
The data must exist in `trim::make_dataset` (§3).

Checks: build the GUI; `tools\windows\gui_smoke_test.bat` and look at the
screenshots in `build\gui_smoke\<n>\`.

## 5. Add an electronic stopping model

1. Append a value to `physics::ElectronicStoppingModel`
   (`common/include/physics/stopping.h`) — append, because the GUI combo box
   binds by index; extend `to_string` / `parse_electronic_stopping_model`.
2. Implement it as a function returning a `stopping_cross_section_quantity`
   for one projectile in one element and add the case to the dispatch in
   `stopping.cpp` (the `switch` used by `layer_stopping`). Bragg's rule,
   compound correction and the `StoppingCurve` tabulation are shared.
3. GUI: add the item to the SRIM-mode stopping combo box in
   `calculation_panel.cpp` (same order); CLI: the `--stopping` help text in
   `main.cpp`.
4. Tests in `physics_test.cpp` against reference values (PSTAR, SRIM,
   published data); docs: `docs/SRIM_MODE.md` §2.2, `docs/CLI.md`.

## 6. Add a CLI option

1. Field in `cli::Options` (`cli/include/cli.h`).
2. Parse it in `parse_arguments` (`cli/src/cli.cpp`); bad values throw
   `UsageError`.
3. Map it to the scenario in `apply_options` (prefer setting a `RunConfig`
   field over special handling in `main.cpp`).
4. Help text in `cli/src/main.cpp`; `docs/CLI.md`.
5. Tests in `test/cli/cli_test.cpp` (`Cli` suite: parsing and mapping);
   an end-to-end `add_test` in `test/CMakeLists.txt` if it changes what a
   command writes.

## 7. Add an output format

1. Value in `output::Format` (`common/include/output/table.h`),
   `to_string` / `parse_format`.
2. Writer in `common/src/output/` taking a shaped `Table`/`Dataset` (units
   and row selection already applied by `output::apply`); call it from
   `output::write_dataset` (`writers.cpp`).
3. CLI `--formats` help, GUI output panel checkbox, `docs/OUTPUT.md`.
4. Tests in `output_test.cpp` (write, read back).

## 8. Regenerate the bundled SRIM stopping tables

Only needed for a new SRIM version or a changed energy grid.

1. On the workstation, in `<SRIM>\_amcset` (a copy of `tools/srim`): run
   steps 1–3 of `tools/srim/README.md` (gen_cases.py → SR Module extraction
   for every pair, ~17 000 tables, several hours → `pack_stopping.py`).
2. `gunzip` the result to `data/srim2013/srim2013_stopping.amcset`; update
   the counts and energies in `data/srim2013/README.md`.
3. Run the tests (`Stopping*`, `TrimRegression`), then the full TRIM
   comparison; update `docs/SRIM_MODE.md` §4 and `docs/VALIDATION.md` §2.

## 9. Add a TRIM reference case

1. Add it to `tools/srim/gen_cases3.py` (or a new batch generator) with a
   descriptive name (`<prefix>_<ion>_<target>_<energy>keV_<kp|fc|...>`;
   replicate seeds as `_s<k>`), regenerate into `<SRIM>\_amcset\trim_cases3`.
2. Ask before running TRIM: `run_trim_cases3.bat` in the SRIM folder runs
   the unfinished cases into `results\trim3\` (TRIM.DAT recoil cases need
   `run_held_cases.bat`, which shows TRIM's message box).
3. Run `tools\windows\compare_with_trim.bat` and update the docs (§10).

## 10. Refresh the validation numbers in the docs

After a full comparison (`build\trim_validation\<set>\analysis.json`):

```bat
python tools\trim_compare\doc_tables.py build\trim_validation\trim3\analysis.json --units
```

Paste the summary and per-unit tables into `docs/VALIDATION.md` §2.2,
update the headline numbers there, in `README.md`, `AGENTS.md` §1 and
`CHANGELOG.md`, and explain any unit beyond |z| = 3.

## 11. Change the native engine

Keep `Simulation::Settings`' three constructor overloads, `print_settings()`
and the conversion in `config::make_engine` (native branch) in sync. Run the
determinism harness before and after (`docs/TESTING.md`); a refactor must be
bit-identical; a physics fix is a separate change with a reference
comparison and an entry in `docs/KNOWN_ISSUES.md` §0.
