// Copyright 2025, Texas A&M University
//
// This file is part of AMCSET.
//
// AMCSET is free software: you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// AMCSET. If not, see <https://www.gnu.org/licenses/>.

#include "app/guard.h"
#include "app/main_window.h"

#include "config.h"
#include "diagnostics.h"
#include "model/output_location.h"
#include "run/datasets.h"
#include "run_config.h"
#include "setup/calculation_panel.h"
#include "setup/crystal_panel.h"
#include "setup/ion_panel.h"
#include "setup/output_panel.h"
#include "setup/target_panel.h"
#include "trim/trim_tables.h"
#include "views/cell_view.h"
#include "views/plots_panel.h"
#include "views/stats_panel.h"
#include "views/tables_panel.h"
#include "views/track_panel.h"
#include "views/track_view.h"

#include <QAction>
#include <QApplication>
#include <QCheckBox>
#include <QCloseEvent>
#include <QDateTime>
#include <QDir>
#include <QPainter>
#include <QTimer>
#include <QDesktopServices>
#include <QDockWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QInputDialog>
#include <QLabel>
#include <QMenuBar>
#include <QMessageBox>
#include <QPlainTextEdit>
#include <QProgressBar>
#include <QPushButton>
#include <QScrollArea>
#include <QSessionManager>
#include <QSettings>
#include <QSpinBox>
#include <QStandardPaths>
#include <QStatusBar>
#include <QStyleHints>
#include <QTabWidget>
#include <QToolBar>
#include <QUrl>

#include <filesystem>

namespace amcset {
namespace gui {

namespace {

QWidget *scrolled(QWidget *w) {
  auto *s = new QScrollArea();
  s->setWidget(w);
  s->setWidgetResizable(true);
  s->setFrameShape(QFrame::NoFrame);
  return s;
}

bool is_results_json(const QString &path) {
  return QFileInfo(path).fileName().compare(QStringLiteral("results.json"), Qt::CaseInsensitive) == 0;
}

} // namespace

MainWindow::MainWindow() : run_(&tracks_) {
  // Views.
  paths_ = new TrackPanel(&tracks_, this);
  plots_ = new PlotsPanel(this);
  tables_ = new TablesPanel(this);
  log_ = new QPlainTextEdit(this);
  log_->setReadOnly(true);
  log_->setMaximumBlockCount(5000);
  central_ = new QTabWidget(this);
  central_->addTab(paths_, tr("Paths"));
  central_->addTab(plots_, tr("Plots"));
  central_->addTab(tables_, tr("Tables"));
  central_->addTab(log_, tr("Messages"));
  central_->setDocumentMode(true);
  setCentralWidget(central_);

  build_actions();
  build_menus();
  build_docks();

  progress_ = new QProgressBar(this);
  progress_->setMaximumWidth(260);
  progress_->setTextVisible(true);
  status_ = new QLabel(tr("Ready"), this);
  hover_ = new QLabel(this);
  checks_ = new QLabel(this);
  checks_->setVisible(false);
  statusBar()->addWidget(status_, 1);
  statusBar()->addPermanentWidget(checks_);
  statusBar()->addPermanentWidget(hover_);
  // Plausibility checks run once the edits settle (they evaluate stopping
  // powers), not on every keystroke.
  checks_timer_ = new QTimer(this);
  checks_timer_->setSingleShot(true);
  checks_timer_->setInterval(300);
  connect(checks_timer_, &QTimer::timeout, this, guarded(this, tr("Checking the setup"), [this] { update_checks(); }));
  statusBar()->addPermanentWidget(progress_);

  connect(&scenario_, &Scenario::changed, this, guarded(this, tr("Main window (changed)"), [this] { scenario_changed(); }));
  connect(&scenario_, &Scenario::file_path_changed, this, guarded(this, tr("Main window (file path changed)"), [this] { update_title(); }));
  connect(scenario_.undo_stack(), &QUndoStack::cleanChanged, this, guarded(this, tr("Main window (cleanChanged)"), [this] { update_title(); }));
  connect(paths_, &TrackPanel::hover_text, hover_, &QLabel::setText);
  connect(paths_->view(), &TrackView::graphics_ready, this,
          guarded(this, tr("Main window (graphics ready)"), [this](const QString &d) { log(tr("Path views: %1").arg(d)); }));
  connect(plots_, &PlotsPanel::hover_text, hover_, &QLabel::setText);
  connect(paths_, &TrackPanel::clear_requested, this, guarded(this, tr("Main window (clear requested)"), [this] {
    tracks_.clear();
    paths_->data_changed();
  }));

  connect(&run_, &RunController::state_changed, this, guarded(this, tr("Main window (state changed)"), [this] { update_actions(); }));
  connect(&run_, &RunController::progress, this, guarded(this, tr("Main window (progress)"), [this](std::size_t done, std::size_t total, double rate) {
    progress_->setMaximum(static_cast<int>(std::min<std::size_t>(total, 2'000'000'000)));
    progress_->setValue(static_cast<int>(std::min<std::size_t>(done, 2'000'000'000)));
    stats_->set_progress(done, total, run_.busy() ? rate : 0.0, run_.elapsed_seconds());
  }));
  connect(&run_, &RunController::dataset_changed, this, guarded(this, tr("Main window (dataset changed)"), [this] {
    const auto d = run_.dataset();
    plots_->set_dataset(d);
    tables_->set_dataset(d);
    stats_->set_dataset(d);
  }));
  connect(&run_, &RunController::tracks_changed, paths_,
          guarded_quiet(tr("Drawing new paths"), [p = paths_] { p->data_changed(); }));
  connect(&run_, &RunController::finished, this, guarded(this, tr("Main window (finished)"), [this](bool ok, const QString &message) {
    status_->setText(message);
    log(message, !ok);
    if (escaped_ > 0) {
      log(tr("%L1 of %L2 ions left through the back of the target (transmitted).").arg(escaped_).arg(run_.completed()),
          true);
    }
    update_actions();
    if (retry_after_stop_) {
      retry_after_stop_ = false;
      if (ok) {
        QTimer::singleShot(0, this, guarded(this, tr("Retrying with the automatic target depth"), [this] { run(); }));
        return;
      }
    }
    if (!ok && !unattended_) QMessageBox::warning(this, tr("Calculation"), message);
  }));
  connect(&run_, &RunController::started, this, guarded(this, tr("Main window (started)"), [this](bool continued) {
    for (const auto &w : run_.warnings()) log(w, true);
    if (!continued) log(tr("Started: %1").arg(QString::fromStdString(run_.config().trim.title)));
    if (!continued && run_.config().mode == config::Mode::srim && run_.config().trim.random_seed) {
      log(tr("Random seed drawn for this run: %1 (enter it as the seed to run the same calculation again)")
              .arg(run_.config().trim.seed));
    }
    update_actions();
  }));
  connect(&run_, &RunController::target_sized, this,
          guarded(this, tr("Main window (target sized)"),
                  [this](const QStringList &notes, bool changed) { apply_sized_target(notes, changed); }));
  connect(&run_, &RunController::ions_escaped, this,
          guarded(this, tr("Main window (ions escaped)"),
                  [this](std::size_t transmitted, std::size_t done) { ions_escaped(transmitted, done); }));
#if QT_CONFIG(sessionmanager)
  // The desktop session ends (log-off, shutdown, an installer's Restart
  // Manager): Qt asks for the data to be saved first (Windows:
  // WM_QUERYENDSESSION).
  connect(qApp, &QGuiApplication::commitDataRequest, this, [this](QSessionManager &m) {
    try {
      session_ending(m);
    } catch (...) {
      report_exception_quietly(QStringLiteral("Ending the session"), std::current_exception());
    }
  });
#endif

  QSettings s;
  restoreGeometry(s.value(QStringLiteral("window/geometry")).toByteArray());
  restoreState(s.value(QStringLiteral("window/state")).toByteArray());
  paths_->restore_settings();
  if (!s.contains(QStringLiteral("window/geometry"))) resize(1500, 950);
  scenario_changed();
  update_actions();
  log(tr("AMCSET %1. Set up a calculation on the left and press Run.").arg(QStringLiteral(AMCSET_VERSION)));
  // Scenarios kept when an earlier session ended with unsaved changes.
  const QStringList kept = QDir(recovery_dir()).entryList({QStringLiteral("scenario-*.json")}, QDir::Files, QDir::Time);
  if (!kept.isEmpty()) {
    log(tr("When the desktop session ended, %n scenario(s) with unsaved changes were kept in %1 (newest: %2). "
           "Open them with File > Open; delete them when no longer needed.",
           nullptr, static_cast<int>(kept.size()))
            .arg(QDir::toNativeSeparators(recovery_dir()), kept.front()),
        true);
  }
}

MainWindow::~MainWindow() = default;

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

void MainWindow::build_actions() {
  const auto act = [this](const QString &text, const QKeySequence &key, QIcon::ThemeIcon icon) {
    auto *a = new QAction(QIcon::fromTheme(icon), text, this);
    if (!key.isEmpty()) a->setShortcut(key);
    return a;
  };
  new_act_ = act(tr("&New"), QKeySequence::New, QIcon::ThemeIcon::DocumentNew);
  open_act_ = act(tr("&Open..."), QKeySequence::Open, QIcon::ThemeIcon::DocumentOpen);
  save_act_ = act(tr("&Save"), QKeySequence::Save, QIcon::ThemeIcon::DocumentSave);
  save_as_act_ = act(tr("Save &As..."), QKeySequence::SaveAs, QIcon::ThemeIcon::DocumentSaveAs);
  import_act_ = act(tr("&Import TRIM.IN..."), {}, QIcon::ThemeIcon::DocumentOpen);
  export_act_ = act(tr("&Export TRIM.IN (and TRIM.DAT)..."), {}, QIcon::ThemeIcon::DocumentSend);
  results_act_ = act(tr("Open &Results (results.json)..."), {}, QIcon::ThemeIcon::DocumentOpen);
  quit_act_ = act(tr("&Quit"), QKeySequence::Quit, QIcon::ThemeIcon::ApplicationExit);
  undo_act_ = scenario_.undo_stack()->createUndoAction(this, tr("&Undo"));
  undo_act_->setShortcut(QKeySequence::Undo);
  undo_act_->setIcon(QIcon::fromTheme(QIcon::ThemeIcon::EditUndo));
  redo_act_ = scenario_.undo_stack()->createRedoAction(this, tr("&Redo"));
  redo_act_->setShortcut(QKeySequence::Redo);
  redo_act_->setIcon(QIcon::fromTheme(QIcon::ThemeIcon::EditRedo));
  run_act_ = act(tr("&Run"), QKeySequence(Qt::Key_F5), QIcon::ThemeIcon::MediaPlaybackStart);
  run_act_->setToolTip(tr("Start a new calculation (F5)"));
  pause_act_ = act(tr("&Pause"), QKeySequence(Qt::Key_F6), QIcon::ThemeIcon::MediaPlaybackPause);
  pause_act_->setCheckable(true);
  stop_act_ = act(tr("S&top"), QKeySequence(Qt::SHIFT | Qt::Key_F5), QIcon::ThemeIcon::MediaPlaybackStop);
  extend_act_ = act(tr("&Continue"), QKeySequence(Qt::Key_F7), QIcon::ThemeIcon::MediaSeekForward);
  extend_act_->setToolTip(tr("Add more ions to the finished calculation, with the current ion settings "
                             "(TRIM's continue, also to superimpose implants of other ions or energies)"));

  connect(new_act_, &QAction::triggered, this, guarded(this, tr("Main window (triggered)"), [this] { new_scenario(); }));
  connect(open_act_, &QAction::triggered, this, guarded(this, tr("Main window (triggered)"), [this] { open_dialog(); }));
  connect(save_act_, &QAction::triggered, this, guarded(this, tr("Main window (triggered)"), [this] { save(); }));
  connect(save_as_act_, &QAction::triggered, this, guarded(this, tr("Main window (triggered)"), [this] { save_as(); }));
  connect(import_act_, &QAction::triggered, this, guarded(this, tr("Main window (triggered)"), [this] { import_trim_in(); }));
  connect(export_act_, &QAction::triggered, this, guarded(this, tr("Main window (triggered)"), [this] { export_trim_in(); }));
  connect(results_act_, &QAction::triggered, this, guarded(this, tr("Main window (triggered)"), [this] { open_results(); }));
  connect(quit_act_, &QAction::triggered, this, &QWidget::close);
  connect(run_act_, &QAction::triggered, this, guarded(this, tr("Main window (triggered)"), [this] { run(); }));
  connect(pause_act_, &QAction::toggled, this, guarded(this, tr("Pause"), [this](bool on) { run_.set_paused(on); }));
  connect(stop_act_, &QAction::triggered, this, guarded(this, tr("Stop"), [this] { run_.stop(); }));
  connect(extend_act_, &QAction::triggered, this, guarded(this, tr("Main window (triggered)"), [this] { extend(); }));

  auto *tb = addToolBar(tr("Main"));
  tb->setObjectName(QStringLiteral("main_toolbar"));
  tb->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
  tb->addAction(new_act_);
  tb->addAction(open_act_);
  tb->addAction(save_act_);
  tb->addSeparator();
  tb->addAction(run_act_);
  tb->addAction(pause_act_);
  tb->addAction(stop_act_);
  tb->addSeparator();
  tb->addAction(extend_act_);
  more_ions_ = new QSpinBox(this);
  more_ions_->setRange(1, 2'000'000'000);
  more_ions_->setValue(1000);
  more_ions_->setPrefix(tr("+"));
  more_ions_->setSuffix(tr(" ions"));
  more_ions_->setGroupSeparatorShown(true);
  tb->addWidget(more_ions_);
}

void MainWindow::build_menus() {
  auto *file = menuBar()->addMenu(tr("&File"));
  file->addAction(new_act_);
  file->addAction(open_act_);
  recent_menu_ = file->addMenu(tr("Open Recent"));
  file->addAction(save_act_);
  file->addAction(save_as_act_);
  file->addSeparator();
  file->addAction(import_act_);
  file->addAction(export_act_);
  file->addSeparator();
  file->addAction(results_act_);
  file->addAction(tr("Output folder for new scenarios..."), this, guarded(this, tr("Menu"), [this] {
    const QString now = output_root();
    QMessageBox box(QMessageBox::Question, tr("Output folder for new scenarios"),
                    tr("New scenarios write their results to %1, and a relative output folder of an unsaved "
                       "scenario is taken relative to it. Choose another folder, or go back to the default (%2)?")
                        .arg(QDir::toNativeSeparators(QDir(now).filePath(QStringLiteral("output"))),
                             QDir::toNativeSeparators(default_output_root())),
                    QMessageBox::Cancel, this);
    auto *choose = box.addButton(tr("Choose..."), QMessageBox::AcceptRole);
    auto *reset = box.addButton(tr("Default"), QMessageBox::ResetRole);
    box.exec();
    if (box.clickedButton() == choose) {
      const QString d = QFileDialog::getExistingDirectory(this, tr("Output folder for new scenarios"), now);
      if (d.isEmpty()) return;
      set_output_root(d);
    } else if (box.clickedButton() == reset) {
      set_output_root({});
    } else {
      return;
    }
    log(tr("New scenarios write to %1.").arg(QDir::toNativeSeparators(QDir(output_root()).filePath(QStringLiteral("output")))));
  }));
  file->addSeparator();
  file->addAction(quit_act_);
  rebuild_recent_menu();

  auto *edit = menuBar()->addMenu(tr("&Edit"));
  edit->addAction(undo_act_);
  edit->addAction(redo_act_);

  auto *runm = menuBar()->addMenu(tr("&Calculation"));
  runm->addAction(run_act_);
  runm->addAction(pause_act_);
  runm->addAction(stop_act_);
  runm->addAction(extend_act_);
  runm->addSeparator();
  runm->addAction(tr("Open output folder"), this, guarded(this, tr("Menu"), [this] {
    const QString d = QString::fromStdString(
        absolute_output_directory(scenario_.config().trim.output_directory, scenario_.file_path()));
    if (!QFileInfo::exists(d)) {
      log(tr("The output folder %1 does not exist yet (it is made by the first run).").arg(d));
      return;
    }
    QDesktopServices::openUrl(QUrl::fromLocalFile(d));
  }));

  auto *view = menuBar()->addMenu(tr("&View"));
  view->setObjectName(QStringLiteral("view_menu"));
  auto *theme = view->addMenu(tr("Colour scheme"));
  const auto scheme = [this](Qt::ColorScheme s) {
    QGuiApplication::styleHints()->setColorScheme(s);
    QSettings().setValue(QStringLiteral("view/color_scheme"), static_cast<int>(s));
    paths_->view()->set_target({}, paths_->view()->incident_z());
  };
  theme->addAction(tr("System"), this, guarded(this, tr("Menu"), [scheme] { scheme(Qt::ColorScheme::Unknown); }));
  theme->addAction(tr("Light"), this, guarded(this, tr("Menu"), [scheme] { scheme(Qt::ColorScheme::Light); }));
  theme->addAction(tr("Dark"), this, guarded(this, tr("Menu"), [scheme] { scheme(Qt::ColorScheme::Dark); }));
  view->addAction(tr("Path memory limit..."), this, guarded(this, tr("Menu"), [this] {
    bool ok = false;
    const int m = QInputDialog::getInt(this, tr("Path memory limit"),
                                       tr("Most path segments kept for drawing (millions; 1 million uses about 72 MB):"),
                                       static_cast<int>(tracks_.max_segments() / 1'000'000), 1, 100, 1, &ok);
    if (!ok) return;
    tracks_.set_max_segments(static_cast<std::size_t>(m) * 1'000'000);
    QSettings().setValue(QStringLiteral("view/max_segments"), m);
  }));
  view->addSeparator();

  auto *help = menuBar()->addMenu(tr("&Help"));
  help->addAction(tr("&About AMCSET"), this, guarded(this, tr("Menu"), [this] {
    QMessageBox::about(this, tr("About AMCSET"),
                       tr("<h3>AMCSET %1</h3><p>Monte Carlo transport of ions and electrons in matter.</p>"
                          "<p>SRIM mode reproduces the algorithms of SRIM-2013's TRIM (Ziegler, Biersack and "
                          "Ziegler). The native engine follows ions and electrons with AMCSET's own "
                          "models.</p><p>Texas A&amp;M University. GNU General Public License v3.</p>")
                           .arg(QStringLiteral(AMCSET_VERSION)));
  }));
  help->addAction(tr("About &Qt"), qApp, &QApplication::aboutQt);
  help->addSeparator();
  auto *detailed = help->addAction(tr("&Detailed diagnostic log"));
  detailed->setCheckable(true);
  detailed->setChecked(QSettings().value(QStringLiteral("diagnostics/detailed"), false).toBool());
  detailed->setToolTip(tr("Also record debug messages (per-thread progress) in the diagnostic log. The log "
                          "only records run-level events, never single collisions, so it costs no speed "
                          "either way."));
  connect(detailed, &QAction::toggled, this, guarded(this, tr("Main window (toggled)"), [](bool on) {
    QSettings().setValue(QStringLiteral("diagnostics/detailed"), on);
    diag::set_level(on ? diag::Level::debug : diag::Level::info);
    AMCSET_LOG(info) << "detailed diagnostic log " << (on ? "on" : "off");
  }));
  help->addAction(tr("Open &log folder"), this, guarded(this, tr("Menu"), [] {
    const auto file = diag::log_file();
    const QString dir = file.empty() ? QString::fromStdString(diag::crash_directory().string())
                                     : QFileInfo(QString::fromStdString(file.string())).absolutePath();
    QDesktopServices::openUrl(QUrl::fromLocalFile(dir));
  }));

  QSettings s;
  const int max_m = s.value(QStringLiteral("view/max_segments"), 4).toInt();
  tracks_.set_max_segments(static_cast<std::size_t>(std::clamp(max_m, 1, 100)) * 1'000'000);
  if (s.contains(QStringLiteral("view/color_scheme"))) {
    QGuiApplication::styleHints()->setColorScheme(static_cast<Qt::ColorScheme>(s.value(QStringLiteral("view/color_scheme")).toInt()));
  }
}

void MainWindow::build_docks() {
  auto *setup = new QTabWidget(this);
  auto *ion = new IonPanel(&scenario_, setup);
  auto *target = new TargetPanel(&scenario_, setup);
  auto *crystal = new CrystalPanel(&scenario_, setup);
  crystal->setObjectName(QStringLiteral("crystal_panel"));
  auto *calc = new CalculationPanel(&scenario_, setup);
  auto *out = new OutputPanel(&scenario_, [this] { return sample_dataset(); }, setup);
  setup->addTab(scrolled(ion), tr("Ion"));
  setup->addTab(target, tr("Target"));
  setup->addTab(scrolled(crystal), tr("Crystal"));
  setup->addTab(scrolled(calc), tr("Calculation"));
  setup->addTab(scrolled(out), tr("Output"));
  const auto msg = [this](const QString &m) {
    status_->setText(m);
    log(m, true);
  };
  connect(ion, &IonPanel::message, this, guarded_quiet(tr("Showing a message"), msg));
  connect(target, &TargetPanel::message, this, guarded_quiet(tr("Showing a message"), msg));
  connect(crystal, &CrystalPanel::message, this, guarded_quiet(tr("Showing a message"), msg));
  connect(calc, &CalculationPanel::message, this, guarded_quiet(tr("Showing a message"), msg));
  connect(out, &OutputPanel::message, this, guarded_quiet(tr("Showing a message"), msg));

  auto *setup_dock = new QDockWidget(tr("Setup"), this);
  setup_dock->setObjectName(QStringLiteral("setup_dock"));
  setup_dock->setWidget(setup);
  setup_dock->setMinimumWidth(380);
  addDockWidget(Qt::LeftDockWidgetArea, setup_dock);

  stats_ = new StatsPanel(this);
  auto *stats_dock = new QDockWidget(tr("Statistics"), this);
  stats_dock->setObjectName(QStringLiteral("stats_dock"));
  stats_dock->setWidget(stats_);
  addDockWidget(Qt::RightDockWidgetArea, stats_dock);

  if (auto *view = menuBar()->findChild<QMenu *>(QStringLiteral("view_menu"))) {
    view->addAction(setup_dock->toggleViewAction());
    view->addAction(stats_dock->toggleViewAction());
  }
  stats_dock->setMinimumWidth(260);
  if (!QSettings().contains(QStringLiteral("window/state"))) {
    QTimer::singleShot(0, this, guarded_quiet(tr("Arranging the docks"), [this, setup_dock, stats_dock] {
      resizeDocks({setup_dock, stats_dock}, {440, 360}, Qt::Horizontal);
    }));
  }
}

// ---------------------------------------------------------------------------
// State
// ---------------------------------------------------------------------------

void MainWindow::scenario_changed() {
  const auto &c = scenario_.config();
  std::vector<double> b{0.0};
  for (const auto &l : c.trim.layers) b.push_back(b.back() + common::in_units(l.width, common::angstrom));
  paths_->set_target(b, c.mode == config::Mode::native && c.native.electron ? 0u : c.trim.ion.z);
  // Recoil atoms start inside the target: the views show their cascades.
  paths_->view()->set_fit_to_paths(c.mode == config::Mode::srim &&
                                   c.trim.source == trim::TrimConfig::Source::particle_recoils);
  plots_->set_output_units(c.trim.outputs);
  tables_->set_output_spec(c.trim.outputs);
  stats_->set_output_units(c.trim.outputs);
  if (!c.trim.layers.empty() && !c.trim.layers.front().elements.empty()) {
    plots_->set_surface_binding_energy(
        common::in_units(c.trim.layers.front().elements.front().surface_binding_energy, common::electron_volt));
  }
  update_title();
  update_actions();
  checks_timer_->start();
}

void MainWindow::update_checks() {
  const auto &c = scenario_.config();
  QStringList w;
  if (c.mode == config::Mode::srim) {
    for (const auto &s : c.trim.plausibility_warnings()) w << QString::fromStdString(s);
  }
  if (w == shown_checks_) return;
  for (const auto &s : w) {
    if (!shown_checks_.contains(s)) log(tr("Setup check: %1").arg(s), true);
  }
  shown_checks_ = w;
  checks_->setVisible(!w.isEmpty());
  checks_->setText(tr("%n setup warning(s)", nullptr, static_cast<int>(w.size())));
  checks_->setToolTip(w.join(QStringLiteral("\n\n")));
}

void MainWindow::update_title() {
  const QString path = scenario_.file_path();
  const QString name = path.isEmpty() ? tr("Untitled") : QFileInfo(path).fileName();
  setWindowTitle(QStringLiteral("%1%2 - AMCSET").arg(name, scenario_.is_modified() ? QStringLiteral(" *") : QString()));
}

void MainWindow::update_actions() {
  const auto st = run_.state();
  const bool busy = run_.busy();
  run_act_->setEnabled(!busy);
  pause_act_->setEnabled(busy && st != RunController::State::stopping);
  {
    QSignalBlocker b(pause_act_);
    pause_act_->setChecked(st == RunController::State::paused);
  }
  stop_act_->setEnabled(busy && st != RunController::State::stopping);
  extend_act_->setEnabled(run_.can_extend());
  more_ions_->setEnabled(!busy);
  for (auto *a : {new_act_, open_act_, import_act_, results_act_}) a->setEnabled(!busy);
  if (run_.sizing()) pause_act_->setEnabled(false); // nothing to pause yet
  switch (st) {
  case RunController::State::running:
    status_->setText(run_.sizing() ? tr("Sizing the target (automatic depth)...") : tr("Running..."));
    break;
  case RunController::State::paused: status_->setText(tr("Paused")); break;
  case RunController::State::stopping: status_->setText(tr("Stopping...")); break;
  case RunController::State::idle: break;
  }
}

void MainWindow::log(const QString &text, bool warning) {
  const QString t = QDateTime::currentDateTime().toString(QStringLiteral("HH:mm:ss"));
  log_->appendPlainText(QStringLiteral("%1  %2%3").arg(t, warning ? tr("! ") : QString(), text));
  if (unattended_) {
    if (warning) qWarning().noquote() << text;
    else qInfo().noquote() << text;
  }
}

std::shared_ptr<const output::Dataset> MainWindow::sample_dataset() const {
  if (auto d = run_.dataset()) return d;
  // No results yet: empty tables of the current target, for the names.
  try {
    auto c = scenario_.config().trim;
    c.stopping_model = physics::ElectronicStoppingModel::builtin;
    c.srim_data_directory.clear();
    c.ion_count = std::max<std::size_t>(c.ion_count, 1);
    if (c.source != trim::TrimConfig::Source::beam && c.particles.empty()) c.source = trim::TrimConfig::Source::beam;
    trim::TrimSimulation sim(c);
    return std::make_shared<output::Dataset>(trim::make_dataset(sim, sim.snapshot(true)));
  } catch (const std::exception &) {
    return nullptr;
  }
}

// ---------------------------------------------------------------------------
// Files
// ---------------------------------------------------------------------------

bool MainWindow::maybe_save() {
  if (!scenario_.is_modified()) return true;
  const auto r = QMessageBox::question(this, tr("Unsaved changes"), tr("Save the changes to the scenario?"),
                                       QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
  if (r == QMessageBox::Cancel) return false;
  if (r == QMessageBox::Save) return save();
  return true;
}

void MainWindow::new_scenario() {
  if (!maybe_save()) return;
  auto c = config::example_config(config::Mode::srim);
  c.trim.trajectory_ions = 200;
  c.trim.output_directory = new_scenario_output_directory();
  scenario_.reset(c);
  log(tr("New scenario."));
}

void MainWindow::open_path(const QString &path) {
  if (is_results_json(path)) {
    output::Dataset d;
    QString err;
    if (!read_results_json(path, d, &err)) {
      QMessageBox::warning(this, tr("Open results"), err);
      return;
    }
    run_.show_dataset(std::make_shared<output::Dataset>(std::move(d)));
    central_->setCurrentWidget(plots_);
    log(tr("Opened results %1").arg(path));
    return;
  }
  const QString err = scenario_.load(path);
  if (!err.isEmpty()) {
    QMessageBox::warning(this, tr("Open"), tr("Cannot open %1:\n%2").arg(path, err));
    return;
  }
  add_recent(path);
  log(tr("Opened %1").arg(path));
  const QString problem = scenario_.problem();
  if (!problem.isEmpty()) {
    log(tr("The scenario cannot run as it is: %1").arg(problem), true);
    QMessageBox::warning(this, tr("Open"),
                         tr("%1 was opened, but it cannot run as it is:\n%2\n\nCorrect it in the setup "
                            "panels before pressing Run.")
                             .arg(QFileInfo(path).fileName(), problem));
  }
}

void MainWindow::open_dialog() {
  if (!maybe_save()) return;
  const QString path = QFileDialog::getOpenFileName(
      this, tr("Open scenario"), QSettings().value(QStringLiteral("dirs/scenario")).toString(),
      tr("AMCSET scenarios and TRIM inputs (*.json TRIM.IN *.IN *.in);;All files (*)"));
  if (path.isEmpty()) return;
  QSettings().setValue(QStringLiteral("dirs/scenario"), QFileInfo(path).absolutePath());
  open_path(path);
}

bool MainWindow::save() {
  if (scenario_.file_path().isEmpty()) return save_as();
  const QString err = scenario_.save(scenario_.file_path());
  if (!err.isEmpty()) {
    QMessageBox::warning(this, tr("Save"), err);
    return false;
  }
  status_->setText(tr("Saved %1").arg(scenario_.file_path()));
  return true;
}

bool MainWindow::save_as() {
  const QString path = QFileDialog::getSaveFileName(this, tr("Save scenario"),
                                                    QSettings().value(QStringLiteral("dirs/scenario")).toString(),
                                                    tr("AMCSET scenario (*.json)"));
  if (path.isEmpty()) return false;
  const QString err = scenario_.save(path);
  if (!err.isEmpty()) {
    QMessageBox::warning(this, tr("Save"), err);
    return false;
  }
  add_recent(path);
  return true;
}

void MainWindow::import_trim_in() {
  if (!maybe_save()) return;
  const QString path = QFileDialog::getOpenFileName(this, tr("Import TRIM.IN"), {}, tr("TRIM input (TRIM.IN *.IN *.in)"));
  if (!path.isEmpty()) open_path(path);
}

void MainWindow::export_trim_in() {
  const QString dir = QFileDialog::getExistingDirectory(this, tr("Folder for TRIM.IN"));
  if (dir.isEmpty()) return;
  const QString err = scenario_.export_trim(dir);
  if (!err.isEmpty()) QMessageBox::warning(this, tr("Export"), err);
  else log(tr("Wrote TRIM.IN to %1").arg(dir));
}

void MainWindow::open_results() {
  const QString path = QFileDialog::getOpenFileName(this, tr("Open results"), {}, tr("AMCSET results (results.json)"));
  if (!path.isEmpty()) open_path(path);
}

void MainWindow::add_recent(const QString &path) {
  QSettings s;
  QStringList r = s.value(QStringLiteral("recent")).toStringList();
  r.removeAll(path);
  r.prepend(path);
  while (r.size() > 10) r.removeLast();
  s.setValue(QStringLiteral("recent"), r);
  rebuild_recent_menu();
}

void MainWindow::rebuild_recent_menu() {
  recent_menu_->clear();
  const QStringList r = QSettings().value(QStringLiteral("recent")).toStringList();
  for (const auto &p : r) {
    recent_menu_->addAction(QFileInfo(p).fileName(), this, guarded(this, tr("Menu"), [this, p] {
      if (maybe_save()) open_path(p);
    }))->setToolTip(p);
  }
  recent_menu_->setEnabled(!r.isEmpty());
}

// ---------------------------------------------------------------------------
// Running
// ---------------------------------------------------------------------------

void MainWindow::run() {
  reset_run_notices();
  QString err;
  // The output folder made absolute (relative ones: to the scenario's
  // folder or the output root, never to the working directory).
  config::RunConfig to_run = scenario_.config();
  to_run.trim.output_directory = absolute_output_directory(to_run.trim.output_directory, scenario_.file_path());
  if (!run_.start(to_run, &err)) {
    log(err, true);
    if (unattended_) QCoreApplication::exit(2);
    else QMessageBox::warning(this, tr("Cannot start"), err);
    return;
  }
  // (RunController::started logs the warnings and "Started".)
  if (run_.sizing()) {
    log(tr("Sizing the target: quick test runs of %L1 ions (automatic target depth)...")
            .arg(scenario_.config().trim.auto_depth.test_ions));
  }
  if (central_->currentWidget() == log_) central_->setCurrentWidget(paths_);
}

void MainWindow::extend() {
  reset_run_notices();
  QString err;
  config::RunConfig to_run = scenario_.config();
  to_run.trim.output_directory = absolute_output_directory(to_run.trim.output_directory, scenario_.file_path());
  if (!run_.extend(to_run, static_cast<std::size_t>(more_ions_->value()), &err)) {
    QMessageBox::warning(this, tr("Cannot continue"), err);
    return;
  }
  log(tr("Continuing with %L1 more ions.").arg(more_ions_->value()));
}

void MainWindow::reset_run_notices() {
  escaped_ = 0;
  escape_asked_ = false;
  retry_after_stop_ = false;
  if (escape_box_) escape_box_->close();
}

void MainWindow::ions_escaped(std::size_t transmitted, std::size_t done) {
  escaped_ = transmitted;
  if (escape_asked_) return;
  escape_asked_ = true;
  const auto &ad = run_.config().trim.auto_depth;
  const bool was_enabled = ad.enabled;
  const std::size_t test_ions = ad.test_ions;
  const std::size_t more_test_ions = std::min<std::size_t>(test_ions * 10, 10'000'000);
  const QString what = tr("Ions are leaving through the back of the target: %L1 of the first %L2 ions were "
                          "transmitted. The calculation continues.")
                           .arg(transmitted)
                           .arg(done);
  log(what, true);
  if (unattended_) return;
  const bool can_retry = !was_enabled || more_test_ions > test_ions;
  QString question;
  if (!was_enabled) {
    question = tr("Stop and run again with the automatic target depth? Quick test runs then make the last layer "
                  "just thick enough that the deepest ion stops in its last %1 %.")
                   .arg(ad.margin * 100.0, 0, 'g', 3);
  } else if (can_retry) {
    question = tr("The target was sized with test runs of %L1 ions; a longer run reaches deeper. Stop and run again "
                  "with test runs of %L2 ions?")
                   .arg(test_ions)
                   .arg(more_test_ions);
  } else {
    question = tr("The target was sized with the largest test runs (%L1 ions); a thicker last layer can be set "
                  "on the Target tab.")
                   .arg(test_ions);
  }
  question += QStringLiteral("\n\n") +
              tr("For a transmission study, continue. The warning can be turned off here or under Calculation, "
                 "Target depth.");
  auto *box = new QMessageBox(QMessageBox::Warning, tr("Ions escaped"), what, QMessageBox::NoButton, this);
  box->setInformativeText(question);
  box->setAttribute(Qt::WA_DeleteOnClose);
  box->setWindowModality(Qt::NonModal);
  QAbstractButton *retry = nullptr;
  if (can_retry) {
    retry = box->addButton(was_enabled ? tr("Retry with more test ions") : tr("Retry with automatic depth"),
                           QMessageBox::AcceptRole);
  }
  box->addButton(tr("Continue"), QMessageBox::RejectRole);
  auto *quiet = new QCheckBox(tr("Do not warn about transmitted ions in this scenario"), box);
  box->setCheckBox(quiet);
  connect(box, &QMessageBox::finished, this,
          guarded(this, tr("Ions escaped (answer)"), [this, box, retry, quiet, was_enabled, test_ions](int) {
            if (quiet->isChecked()) {
              const QString err = scenario_.modify(tr("Warn when ions escape"), [](config::RunConfig &c) {
                c.trim.auto_depth.warn_transmitted = false;
              }, -1);
              if (!err.isEmpty()) log(err, true);
            }
            if (retry != nullptr && box->clickedButton() == retry) retry_with_auto_depth(was_enabled, test_ions);
          }));
  escape_box_ = box;
  box->show();
}

void MainWindow::retry_with_auto_depth(bool was_enabled, std::size_t test_ions) {
  const QString err = scenario_.modify(tr("Automatic target depth"), [&](config::RunConfig &c) {
    auto &ad = c.trim.auto_depth;
    if (was_enabled) ad.test_ions = std::min<std::size_t>(test_ions * 10, 10'000'000);
    ad.enabled = true;
  }, -1);
  if (!err.isEmpty()) {
    log(err, true);
    QMessageBox::warning(this, tr("Automatic target depth"), err);
    return;
  }
  log(tr("Running again with the automatic target depth."));
  if (run_.busy()) {
    retry_after_stop_ = true; // finished() starts the new run
    run_.stop();
  } else {
    run();
  }
}

void MainWindow::apply_sized_target(const QStringList &notes, bool changed) {
  for (const auto &n : notes) log(n);
  const auto &sized = run_.config().trim;
  const auto &now = scenario_.config().trim;
  if (!changed && sized.window_end == now.window_end) return;
  // Show the thickness that is simulated in the setup (undoable).
  const QString err = scenario_.modify(tr("Automatic target depth"), [&sized](config::RunConfig &c) {
    if (c.trim.layers.size() != sized.layers.size() || sized.layers.empty()) {
      throw std::invalid_argument("the target was changed while it was sized");
    }
    c.trim.layers.back().width = sized.layers.back().width;
    c.trim.window_end = sized.window_end;
  }, -1);
  if (!err.isEmpty()) log(tr("The sized thickness is not shown in the setup: %1").arg(err), true);
}

void MainWindow::save_screenshots(const QString &dir) {
  QDir().mkpath(dir);
  grab().save(QDir(dir).filePath(QStringLiteral("window.png")));
  for (int i = 0; i < central_->count(); ++i) {
    central_->setCurrentIndex(i);
    QCoreApplication::processEvents();
    QWidget *w = central_->widget(i);
    QImage img;
    if (w == paths_) {
      img = paths_->view()->grabFramebuffer();
      if (!img.isNull()) {
        img = img.convertToFormat(QImage::Format_ARGB32_Premultiplied);
        img.setDevicePixelRatio(paths_->view()->devicePixelRatioF());
        if (auto *overlay = paths_->view()->findChild<QWidget *>()) {
          QPainter p(&img);
          overlay->render(&p, QPoint(), QRegion(), QWidget::DrawChildren);
        }
      }
    } else {
      img = w->grab().toImage();
    }
    img.save(QDir(dir).filePath(QStringLiteral("tab%1.png").arg(i)));
  }
  // Every plot, not only the selected one.
  central_->setCurrentWidget(plots_);
  QCoreApplication::processEvents();
  plots_->save_plot_images(dir);
  central_->setCurrentIndex(0);
  if (auto *dock = findChild<QDockWidget *>(QStringLiteral("setup_dock"))) {
    if (auto *tabs = qobject_cast<QTabWidget *>(dock->widget())) {
      for (int i = 0; i < tabs->count(); ++i) {
        tabs->setCurrentIndex(i);
        QCoreApplication::processEvents();
        dock->grab().save(QDir(dir).filePath(QStringLiteral("setup%1.png").arg(i)));
        // A tab that scrolls in the dock: also the whole of it.
        if (auto *area = qobject_cast<QScrollArea *>(tabs->widget(i)); area != nullptr && area->widget() != nullptr) {
          area->widget()->grab().save(QDir(dir).filePath(QStringLiteral("setup%1_full.png").arg(i)));
        }
      }
      tabs->setCurrentIndex(0);
    }
  }
  // The whole Crystal tab (it scrolls in the dock) and its view of the cell.
  if (auto *crystal = findChild<CrystalPanel *>(QStringLiteral("crystal_panel"))) {
    crystal->grab().save(QDir(dir).filePath(QStringLiteral("crystal_tab.png")));
    crystal->cell_view()->grab().save(QDir(dir).filePath(QStringLiteral("crystal_cell.png")));
  }
}

void MainWindow::run_unattended(const QString &dir, bool quit, int timeout_s) {
  unattended_ = true;
  // What was logged before (e.g. the graphics device of the path views).
  if (!log_->toPlainText().isEmpty()) qInfo().noquote() << log_->toPlainText();
  if (timeout_s > 0) {
    QTimer::singleShot(std::chrono::seconds(timeout_s), this, guarded(this, tr("Time-out"), [this] {
      log(tr("Timed out."), true);
      run_.stop();
      QCoreApplication::exit(3);
    }));
  }
  connect(&run_, &RunController::finished, this, guarded(this, tr("Main window (finished)"), [this, dir, quit](bool ok, const QString &message) {
    QTimer::singleShot(500, this, guarded(this, tr("Saving the views"), [this, dir, quit, ok, message] {
      if (!ok) log(message, true);
      if (!dir.isEmpty()) {
        save_screenshots(dir);
        log(tr("Saved images of the views to %1").arg(QDir::toNativeSeparators(dir)));
      }
      if (quit) QCoreApplication::exit(ok ? 0 : 1);
    }));
  }));
  // Start once the window is shown, so the path views have a graphics device.
  QTimer::singleShot(0, this, guarded(this, tr("Starting the calculation"), [this] { run(); }));
}

void MainWindow::closeEvent(QCloseEvent *e) {
  // Called by Qt: no exception may leave. If saving the settings fails, the
  // window still closes.
  try {
    close_event(e);
  } catch (...) {
    report_exception_quietly(QStringLiteral("Closing the window"), std::current_exception());
    e->accept();
  }
}

void MainWindow::close_event(QCloseEvent *e) {
  if (quiet_quit_) {
    // A signal or the session's end: no questions (the session's end has
    // asked already, where it could).
    if (run_.busy()) run_.stop();
    save_window_settings();
    e->accept();
    return;
  }
  if (run_.busy()) {
    const auto r = QMessageBox::question(this, tr("Calculation running"),
                                         tr("A calculation is running. Stop it and quit?"));
    if (r != QMessageBox::Yes) {
      e->ignore();
      return;
    }
    run_.stop();
  }
  if (!maybe_save()) {
    e->ignore();
    return;
  }
  save_window_settings();
  e->accept();
}

void MainWindow::save_window_settings() {
  QSettings s;
  s.setValue(QStringLiteral("window/geometry"), saveGeometry());
  s.setValue(QStringLiteral("window/state"), saveState());
  paths_->save_settings();
}

void MainWindow::quit_quietly(int exit_code) {
  // SIGTERM, SIGHUP or SIGINT (POSIX): stop the calculation and quit without
  // dialogs; the window's layout is kept, unsaved scenario edits are not.
  AMCSET_LOG(info) << "GUI: quitting (exit code " << exit_code << ")";
  quiet_quit_ = true;
  if (run_.busy()) run_.stop();
  close();
  QCoreApplication::exit(exit_code);
}

QString MainWindow::recovery_dir() {
  return QStandardPaths::writableLocation(QStandardPaths::AppLocalDataLocation) + QStringLiteral("/recovery");
}

#if QT_CONFIG(sessionmanager)
void MainWindow::session_ending(QSessionManager &) {
  // Log-off, shutdown, an installer's Restart Manager: never block or veto
  // it (Windows' rules for desktop programs: answer at once, do not cancel).
  // A scenario with unsaved changes is kept in the recovery folder instead
  // of asking; the next start says where it is.
  AMCSET_LOG(info) << "GUI: the desktop session ends";
  if (scenario_.is_modified()) {
    const QString dir = recovery_dir();
    QDir().mkpath(dir);
    const QString file = dir + QStringLiteral("/scenario-%1.json")
                                   .arg(QDateTime::currentDateTime().toString(QStringLiteral("yyyyMMdd-HHmmss")));
    try {
      config::save_config(scenario_.config(), std::filesystem::path(file.toStdWString()));
      AMCSET_LOG(info) << "GUI: unsaved scenario kept as " << file.toStdString();
    } catch (...) {
      report_exception_quietly(QStringLiteral("Keeping the unsaved scenario"), std::current_exception());
    }
  }
  if (run_.busy()) run_.stop();
  save_window_settings();
  quiet_quit_ = true; // the window then closes without asking
}
#endif

} // namespace gui
} // namespace amcset
