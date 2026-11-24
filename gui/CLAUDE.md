# gui/ — amcset-gui (Qt 6 Widgets + QRhi)

User guide and code layout: `docs/GUI.md`. Architecture and threading:
`docs/ARCHITECTURE.md` §5.

## Structure (`gui/src/`)

| Folder | Role | Key types |
|---|---|---|
| `app/` | Window, actions, menus, docks; wires everything; unattended mode (`--run`, `--screenshots`, `--quit`, `--timeout`) | `MainWindow` |
| `model/` | The edited scenario and its table models | `Scenario` (RunConfig + `QUndoStack`; `modify()` is the only way to change it), `LayerTableModel`/`ElementTableModel` (`target_models.h`), unit-aware text (`quantity_text.h`) |
| `run/` | Running a calculation | `RunController` (engine on a `std::thread`, GUI-thread `QTimer` poll, live `output::Dataset`, pause/stop/continue), `TrackStore` (path segments), `datasets.cpp` (native-engine and `results.json` datasets) |
| `setup/` | Setup panels and dialogs | `IonPanel`, `TargetPanel`, `CrystalPanel` (a layer's unit cell, orientation, vibrations; the crystal transport's settings), `CalculationPanel` (incl. the TRIM bug checkboxes and tooltips), `OutputPanel`, `OutputSpecDialog`, `PeriodicTableDialog`, `FormBinder` |
| `views/` | Results, and the unit cell | `TrackView` (a `QRhiWidget`)/`TrackPanel` (path views), `PlotWidget`/`PlotsPanel` (TRIM's plots + custom), `TablesPanel`, `StatsPanel`, `CellView` (the unit cell in 3D, QPainter; used by `CrystalPanel`) |

`gui/shaders/`: GLSL sources and pre-baked `.qsb` in `baked/` (used by
default; `-DAMCSET_COMPILE_SHADERS=ON` compiles them with qsb; regenerate
with `shaders/bake.sh` after editing).

## Rules

- **The scenario is the model.** Widgets never hold state of their own:
  read from `Scenario::config()`, write through `Scenario::modify(what,
  edit, merge_id)` (undoable, validated: if `edit` throws, the change is
  rejected and the message shown). Bind simple widgets with `FormBinder`
  (`text`, `check`, `combo`, `spin`, …); it refreshes them on
  `Scenario::changed` with signals blocked.
- **No physics or file formats here.** Use `config::` (scenario JSON,
  `make_engine`), `trim::` (TRIM.IN/TRIM.DAT, `make_dataset`) and
  `output::` (tables, units, writers). A value the GUI shows must come from
  an `output::Dataset` or the scenario.
- **Threading.** The engine runs on `RunController`'s `std::thread`. Only the
  GUI thread touches widgets. The controller polls `completed_count()`,
  `snapshot_tallies()`/TRIM snapshots and `try_pop_trajectory()` from a
  `QTimer`; keep each tick short (`drain_tracks(budget_ms)`). Never block the
  GUI thread on the engine.
- **Exceptions** must not escape slots: connect through `guarded(context,
  what, f)` / `guarded_quiet(what, f)` (`app/guard.h`), wrap the bodies of
  virtual overrides (`data()`, `paintEvent`, ...) in try/catch with a safe
  return; `Application::notify` is only the last net. Catch, then show (dialog or the log
  dock; in unattended mode the log on stderr).
- **Naming:** `snake_case` methods and members ending in `_` (as in
  `common/`), Qt signal/slot connections with function pointers.
- **New result plot:** add a `PlotDef` in `views/plots_panel.cpp` referring
  to a table key and column keys of the dataset (wildcards like
  `recoils_*`); new data must first exist as a table in
  `trim::make_dataset`.
- **New setup field:** add it to RunConfig first (`common/CLAUDE.md`), then
  bind it in the panel with `FormBinder`; add a tooltip that says what the
  value does and its TRIM equivalent.
- **Crystals.** `CrystalPanel` changes a layer's crystal only through its
  `edit()` / `set_cell()`, which end with `trim::apply_crystal(layer)`: the
  layer's elements and density must stay the unit cell's or validation
  refuses the change. Geometry (Miller indices, tilt and twist, nearest
  direction, cell files) comes from `crystal::`; the panel and `CellView`
  compute none of it themselves. `CellView` holds no scenario state except
  the direction of view.

## Checking GUI changes

There are no GUI unit tests. Build, then run
`tools\windows\gui_smoke_test.bat [scenario]`: it runs the scenario, saves a
screenshot of every view and every plot (`plot_<id>.png`) to
`build\gui_smoke\<n>\` and writes
`exit_code.txt` (0 = finished, 1 = run failed, 2 = could not start,
3 = timeout, 10+n = shutdown hung). Look at the screenshots. The graphics API
can be forced with `AMCSET_RHI=opengl|vulkan|d3d11|d3d12`.

Without opening a window (everything except the path views, which need a
graphics device): after `conanrun.bat`, `set QT_QPA_PLATFORM=offscreen`,
`set QT_QPA_FONTDIR=C:\Windows\Fonts` (the offscreen platform has no fonts
of its own), then `amcset-gui.exe --screenshots DIR --quit --timeout 240
scenario.json`. A scenario with a crystal layer gives `crystal_tab.png` and
`crystal_cell.png` their content; give it its own `output.directory`.

The GUI cannot be built against Ubuntu 24.04's Qt 6.4 (needs ≥ 6.7); on
Linux use the Conan build.
