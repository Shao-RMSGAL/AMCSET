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

/*! \file main_window.h
 *  \brief The application window: it owns the scenario, the run controller
 *  and the path store and connects them to the setup panels and the views.
 */

#pragma once

#include "model/scenario.h"
#include "run/run_controller.h"
#include "run/track_store.h"

#include <QMainWindow>
#include <QPointer>
#include <QStringList>

class QAction;
class QMessageBox;
class QTimer;
class QLabel;
class QMenu;
class QPlainTextEdit;
class QProgressBar;
class QSpinBox;
class QTabWidget;

namespace amcset {
namespace gui {

class PlotsPanel;
class StatsPanel;
class TablesPanel;
class TrackPanel;

class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  MainWindow();
  ~MainWindow() override;
  //! Open a scenario, TRIM.IN or results.json given on the command line.
  void open_path(const QString &path);
  /*!
   * \brief Start the calculation; when it finishes, save an image of every
   * view to \p screenshot_dir (if not empty) and quit if \p quit (for
   * automated checks of the GUI). Messages go to the log on stderr instead
   * of dialogs; after \p timeout_s seconds (0 = none) the program exits with
   * code 3.
   */
  void run_unattended(const QString &screenshot_dir, bool quit, int timeout_s = 0);

protected:
  void closeEvent(QCloseEvent *e) override;

private:
  void close_event(QCloseEvent *e);
  void build_actions();
  void build_menus();
  void build_docks();
  void scenario_changed();
  void update_title();
  void update_actions();
  void log(const QString &text, bool warning = false);
  bool maybe_save();
  void new_scenario();
  void open_dialog();
  bool save();
  bool save_as();
  void import_trim_in();
  void export_trim_in();
  void open_results();
  void run();
  void extend();
  void reset_run_notices();
  //! RunController::ions_escaped: log, and once per run offer a retry.
  void ions_escaped(std::size_t transmitted, std::size_t done);
  //! Enable the automatic target depth (or more test ions) and run again.
  void retry_with_auto_depth(bool was_enabled, std::size_t test_ions);
  //! RunController::target_sized: show the sized thickness in the setup.
  void apply_sized_target(const QStringList &notes, bool changed);
  //! The scenario's plausibility warnings in the status bar and the log.
  void update_checks();
  void add_recent(const QString &path);
  void rebuild_recent_menu();
  std::shared_ptr<const output::Dataset> sample_dataset() const;
  void save_screenshots(const QString &dir);

  Scenario scenario_;
  TrackStore tracks_;
  RunController run_;

  TrackPanel *paths_ = nullptr;
  PlotsPanel *plots_ = nullptr;
  TablesPanel *tables_ = nullptr;
  StatsPanel *stats_ = nullptr;
  QPlainTextEdit *log_ = nullptr;
  QTabWidget *central_ = nullptr;
  QProgressBar *progress_ = nullptr;
  QLabel *status_ = nullptr;
  QLabel *hover_ = nullptr;
  QLabel *checks_ = nullptr;        //!< "n setup warnings" (tooltip: the list)
  QTimer *checks_timer_ = nullptr;  //!< update_checks() after edits settle
  QStringList shown_checks_;        //!< warnings already in the log
  QSpinBox *more_ions_ = nullptr;
  QMenu *recent_menu_ = nullptr;

  QAction *new_act_ = nullptr;
  QAction *open_act_ = nullptr;
  QAction *save_act_ = nullptr;
  QAction *save_as_act_ = nullptr;
  QAction *import_act_ = nullptr;
  QAction *export_act_ = nullptr;
  QAction *results_act_ = nullptr;
  QAction *quit_act_ = nullptr;
  QAction *undo_act_ = nullptr;
  QAction *redo_act_ = nullptr;
  QAction *run_act_ = nullptr;
  QAction *pause_act_ = nullptr;
  QAction *stop_act_ = nullptr;
  QAction *extend_act_ = nullptr;

  bool unattended_ = false; //!< no dialogs; log to stderr; quit on errors
  // Ions transmitted during the current run (RunController::ions_escaped).
  std::size_t escaped_ = 0;
  bool escape_asked_ = false;          //!< the retry question was shown this run
  bool retry_after_stop_ = false;      //!< run() again when the stopped run finishes
  QPointer<QMessageBox> escape_box_;   //!< the open retry question
};

} // namespace gui
} // namespace amcset
