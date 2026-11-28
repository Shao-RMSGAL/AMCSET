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

/*! \file run_controller.h
 *  \brief Runs a calculation on a worker thread and feeds the views.
 *
 * The controller owns the engine and its thread. A GUI-thread timer polls
 * the engine for progress, moves finished trajectories into the TrackStore
 * and rebuilds the result dataset from a snapshot of the tallies, so the
 * plots and statistics are live while the calculation runs. A finished
 * SRIM-mode calculation can be extended with more ions (TRIM's "continue"),
 * also with a different ion, energy or angle on the same target (TRIM's
 * superposition of several implants).
 */

#pragma once

#include "output/table.h"
#include "run_config.h"
#include "trim/trim.h"

#include <QElapsedTimer>
#include <QObject>
#include <QStringList>
#include <QTimer>

#include <atomic>
#include <exception>
#include <memory>
#include <thread>

namespace amcset {
namespace gui {

class TrackStore;

class RunController : public QObject {
  Q_OBJECT

public:
  enum class State { idle, running, paused, stopping };

  explicit RunController(TrackStore *tracks, QObject *parent = nullptr);
  ~RunController() override;

  /*!
   * \brief Start a new calculation (clears the paths and the previous
   * results). With the automatic target depth (SRIM mode, beam source,
   * auto_depth.enabled) quick test runs first size the last layer on the
   * worker thread (sizing() is true meanwhile; target_sized() reports the
   * result), then the calculation runs on the sized target.
   */
  bool start(const config::RunConfig &config, QString *error);
  /*!
   * \brief Add \p more_ions ions to the finished SRIM-mode calculation, with
   * the ion and settings of \p config (the target must be the same).
   */
  bool extend(const config::RunConfig &config, std::size_t more_ions, QString *error);
  void set_paused(bool paused);
  void stop();

  State state() const noexcept { return state_; }
  bool busy() const noexcept { return state_ != State::idle; }
  //! The test runs of the automatic target depth are running (state running).
  bool sizing() const noexcept { return sizing_; }
  //! Whether extend() can be called (a SRIM-mode calculation is finished).
  bool can_extend() const noexcept { return state_ == State::idle && has_previous_; }
  //! Results of the running or last calculation (null before the first one).
  std::shared_ptr<const output::Dataset> dataset() const { return dataset_; }
  //! Replace the shown results (e.g. a results.json opened from disk).
  void show_dataset(std::shared_ptr<const output::Dataset> d);
  //! The configuration of the running or last calculation.
  const config::RunConfig &config() const noexcept { return config_; }
  std::size_t completed() const noexcept { return completed_; }
  std::size_t total() const noexcept { return total_; }
  double elapsed_seconds() const;
  QStringList warnings() const { return warnings_; }
  //! The last run's GPU report (trim::TrimSimulation::ComputeReport::lines(),
  //! the CLI's lines); empty after a run on the CPU. Valid in finished().
  QStringList compute_lines() const { return compute_lines_; }

signals:
  void state_changed(amcset::gui::RunController::State state);
  void progress(std::size_t done, std::size_t total, double per_second);
  void dataset_changed();
  void tracks_changed();
  void finished(bool ok, const QString &message);
  //! The engine started (\p continued: extend()); warnings() are its warnings.
  void started(bool continued);
  /*!
   * \brief The automatic target depth sized the target: config() now has
   * the new last-layer thickness (and depth window); \p notes say what was
   * done. Emitted before the calculation starts on that target.
   */
  void target_sized(const QStringList &notes, bool changed);
  /*!
   * \brief Ions have left through the back of the target (SRIM mode, beam
   * source, auto_depth.warn_transmitted): \p transmitted of \p done so far.
   * Emitted when the count grows (at most every poll); the run continues.
   */
  void ions_escaped(std::size_t transmitted, std::size_t done);

private:
  bool launch(std::unique_ptr<common::SimulationEngine> engine, QString *error);
  void finish_sizing();
  void fail_before_start(const QString &message, bool ok);
  void poll();
  void finish();
  void set_state(State s);
  void refresh_dataset(bool final_results);
  bool drain_tracks(int budget_ms);

  TrackStore *tracks_;
  config::RunConfig config_;
  std::unique_ptr<common::SimulationEngine> engine_;
  trim::TrimSimulation *trim_ = nullptr; // engine_ if SRIM mode
  std::thread thread_;
  std::atomic<bool> done_{false};
  std::exception_ptr error_;
  State state_ = State::idle;
  QTimer timer_;
  QElapsedTimer clock_;
  qint64 last_snapshot_ms_ = -1000000;
  int refresh_failures_ = 0; // live-dataset refreshes that failed this run
  std::size_t escaped_reported_ = 0; // transmitted ions last reported (ions_escaped)
  // Automatic target depth. While sizing_, the worker thread owns
  // sizing_config_ and sizing_result_; the GUI thread reads them after the
  // join in finish_sizing(). engine_ is null meanwhile.
  bool sizing_ = false;
  std::atomic<bool> sizing_cancel_{false};
  config::RunConfig sizing_config_;
  trim::AutoDepthResult sizing_result_;
  std::size_t completed_ = 0;
  std::size_t total_ = 0;
  double elapsed_offset_ = 0.0;
  std::shared_ptr<const output::Dataset> dataset_;
  QStringList warnings_;
  QStringList compute_lines_;
  // For extend().
  bool has_previous_ = false;
  trim::TrimTallies previous_;
  std::size_t previous_ions_ = 0;
  double previous_seconds_ = 0.0;
};

} // namespace gui
} // namespace amcset
