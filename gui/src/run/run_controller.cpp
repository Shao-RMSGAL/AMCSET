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
#include "run/run_controller.h"

#include "amcset_common.h"
#include "diagnostics.h"
#include "run/datasets.h"
#include "run/track_store.h"
#include "trim/trim_tables.h"

#include <QFileInfo>

#include <system_error>

namespace amcset {
namespace gui {

namespace {

constexpr int poll_interval_ms = 100;
constexpr qint64 snapshot_interval_ms = 400;

// Everything that must not change when a calculation is extended.
std::string target_signature(const config::RunConfig &c) {
  config::RunConfig s;
  s.mode = c.mode;
  s.trim.ion = c.trim.ion;
  s.trim.layers = c.trim.layers;
  s.trim.window_start = c.trim.window_start;
  s.trim.window_end = c.trim.window_end;
  s.trim.depth_bins = c.trim.depth_bins;
  s.trim.damage = c.trim.damage;
  s.trim.source = c.trim.source;
  auto j = config::to_json(s);
  // The ion may change (TRIM superimposes implants); everything else may not.
  const auto a = j.find("\"ion\"");
  const auto b = j.find("\"layers\"");
  if (a != std::string::npos && b != std::string::npos && b > a) j.erase(a, b - a);
  return j;
}

QString what_of(const std::exception_ptr &p) {
  try {
    if (p) std::rethrow_exception(p);
  } catch (const std::exception &e) {
    return QString::fromUtf8(e.what());
  } catch (...) {
    return QStringLiteral("unknown error");
  }
  return {};
}

} // namespace

RunController::RunController(TrackStore *tracks, QObject *parent)
    : QObject(parent), tracks_(tracks) {
  timer_.setInterval(poll_interval_ms);
  connect(&timer_, &QTimer::timeout, this, guarded(nullptr, tr("Calculation (timeout)"), [this] { poll(); }));
}

RunController::~RunController() {
  sizing_cancel_.store(true);
  if (engine_) engine_->request_cancel();
  if (thread_.joinable()) thread_.join();
}

void RunController::set_state(State s) {
  if (s == state_) return;
  state_ = s;
  emit state_changed(s);
}

double RunController::elapsed_seconds() const {
  return elapsed_offset_ + (busy() ? static_cast<double>(clock_.elapsed()) / 1000.0 : 0.0);
}

bool RunController::start(const config::RunConfig &config, QString *error) {
  if (busy()) {
    if (error) *error = tr("A calculation is already running.");
    return false;
  }
  const bool size_first = config.mode == config::Mode::srim && config.trim.auto_depth.enabled &&
                          config.trim.source == trim::TrimConfig::Source::beam;
  std::unique_ptr<common::SimulationEngine> engine;
  try {
    if (size_first) {
      config.trim.validate(); // errors now, not after the test runs
    } else {
      engine = config::make_engine(config);
    }
  } catch (const std::exception &e) {
    if (error) *error = QString::fromUtf8(e.what());
    return false;
  }
  config_ = config;
  has_previous_ = false;
  previous_ions_ = 0;
  previous_seconds_ = 0.0;
  elapsed_offset_ = 0.0;
  tracks_->clear();
  emit tracks_changed();
  if (!size_first) return launch(std::move(engine), error);

  // Automatic target depth: the test runs go on the worker thread.
  if (thread_.joinable()) thread_.join();
  engine_.reset();
  trim_ = nullptr;
  sizing_config_ = config;
  sizing_result_ = trim::AutoDepthResult{};
  sizing_cancel_.store(false);
  done_.store(false);
  error_ = nullptr;
  warnings_.clear();
  completed_ = 0;
  total_ = config.trim.ion_count;
  AMCSET_LOG(info) << "GUI: sizing the target (automatic depth, " << config.trim.auto_depth.test_ions
                   << " test ions)";
  try {
    thread_ = std::thread([this] {
      try {
        sizing_result_ = trim::size_target(sizing_config_.trim, &sizing_cancel_);
      } catch (const std::exception &e) {
        AMCSET_LOG(error) << "GUI: automatic target depth failed: " << e.what();
        error_ = std::current_exception();
      } catch (...) {
        AMCSET_LOG(error) << "GUI: automatic target depth failed with an unknown exception";
        error_ = std::current_exception();
      }
      done_.store(true);
    });
  } catch (const std::system_error &e) {
    AMCSET_LOG(error) << "GUI: cannot start the sizing thread: " << e.what();
    if (error) *error = tr("Cannot start the calculation thread: %1").arg(QString::fromUtf8(e.what()));
    return false;
  }
  sizing_ = true;
  clock_.start();
  set_state(State::running);
  timer_.start();
  return true;
}

void RunController::finish_sizing() {
  if (thread_.joinable()) thread_.join();
  sizing_ = false;
  if (error_) {
    fail_before_start(tr("The automatic target depth failed: %1").arg(what_of(error_)), false);
    return;
  }
  QStringList notes;
  for (const auto &n : sizing_result_.notes) notes << QString::fromStdString(n);
  if (sizing_cancel_.load()) {
    fail_before_start(tr("Stopped while sizing the target; the calculation did not start."), true);
    return;
  }
  config_ = sizing_config_;
  emit target_sized(notes, sizing_result_.changed);
  std::unique_ptr<common::SimulationEngine> engine;
  QString err;
  try {
    engine = config::make_engine(config_);
  } catch (const std::exception &e) {
    err = QString::fromUtf8(e.what());
  }
  if (!engine || !launch(std::move(engine), &err)) {
    fail_before_start(err.isEmpty() ? tr("The calculation could not start.") : err, false);
    return;
  }
  emit state_changed(state_); // still running, but no longer sizing
}

void RunController::fail_before_start(const QString &message, bool ok) {
  timer_.stop();
  sizing_ = false;
  engine_.reset();
  trim_ = nullptr;
  has_previous_ = false;
  AMCSET_LOG(info) << "GUI: " << message.toStdString();
  set_state(State::idle);
  emit progress(0, total_, 0.0);
  emit finished(ok, message);
}

bool RunController::extend(const config::RunConfig &config, std::size_t more_ions, QString *error) {
  if (!can_extend()) {
    if (error) *error = tr("There is no finished SRIM-mode calculation to continue.");
    return false;
  }
  if (config.mode != config::Mode::srim || target_signature(config) != target_signature(config_)) {
    if (error) {
      *error = tr("The target, depth window, damage calculation or particle source changed. "
                  "A calculation can only be continued on the same target; start a new one.");
    }
    return false;
  }
  config::RunConfig c = config;
  c.trim.first_ion = previous_ions_;
  c.trim.ion_count = more_ions;
  std::unique_ptr<common::SimulationEngine> engine;
  try {
    auto sim = std::make_unique<trim::TrimSimulation>(c.trim);
    sim->continue_from(previous_);
    engine = std::move(sim);
  } catch (const std::exception &e) {
    if (error) *error = QString::fromUtf8(e.what());
    return false;
  }
  config_ = c;
  elapsed_offset_ = previous_seconds_;
  return launch(std::move(engine), error);
}

bool RunController::launch(std::unique_ptr<common::SimulationEngine> engine, QString *error) {
  // A previous worker must be gone before its engine is replaced (assigning
  // to a joinable std::thread would terminate the program).
  if (thread_.joinable()) thread_.join();
  engine_ = std::move(engine);
  trim_ = dynamic_cast<trim::TrimSimulation *>(engine_.get());
  refresh_failures_ = 0;
  // (A continued run starts with the earlier part's transmitted ions.)
  escaped_reported_ = trim_ != nullptr ? trim_->transmitted_so_far() : 0;
  warnings_.clear();
  if (trim_ != nullptr) {
    for (const auto &w : trim_->warnings()) warnings_ << QString::fromStdString(w);
  }
  completed_ = 0;
  total_ = engine_->total_count();
  done_.store(false);
  error_ = nullptr;
  clock_.start();
  last_snapshot_ms_ = -1000000;
  AMCSET_LOG(info) << "GUI: calculation started (" << total_ << " particles, "
                   << (trim_ != nullptr ? "SRIM mode" : "native engine") << ")";
  for (const auto &w : warnings_) AMCSET_LOG(warning) << "GUI: " << w.toStdString();
  try {
    thread_ = std::thread([this] {
      try {
        engine_->run();
      } catch (const std::exception &e) {
        AMCSET_LOG(error) << "GUI: calculation failed: " << e.what();
        error_ = std::current_exception();
      } catch (...) {
        AMCSET_LOG(error) << "GUI: calculation failed with an unknown exception";
        error_ = std::current_exception();
      }
      done_.store(true);
    });
  } catch (const std::system_error &e) {
    AMCSET_LOG(error) << "GUI: cannot start the calculation thread: " << e.what();
    if (error) *error = tr("Cannot start the calculation thread: %1").arg(QString::fromUtf8(e.what()));
    engine_.reset();
    trim_ = nullptr;
    set_state(State::idle);
    return false;
  }
  set_state(State::running);
  timer_.start();
  refresh_dataset(false);
  emit started(has_previous_);
  return true;
}

void RunController::set_paused(bool paused) {
  if (!engine_ || !busy() || state_ == State::stopping) return;
  if (!engine_->supports_pause()) return;
  engine_->set_paused(paused);
  set_state(paused ? State::paused : State::running);
}

void RunController::stop() {
  if (sizing_ && state_ == State::running) {
    AMCSET_LOG(info) << "GUI: stop requested while sizing the target";
    sizing_cancel_.store(true);
    set_state(State::stopping);
    return;
  }
  if (!engine_ || !busy()) return;
  AMCSET_LOG(info) << "GUI: stop requested after " << engine_->completed_count() << " particles";
  engine_->request_cancel();
  set_state(State::stopping);
}

void RunController::show_dataset(std::shared_ptr<const output::Dataset> d) {
  if (busy()) return;
  dataset_ = std::move(d);
  has_previous_ = false;
  emit dataset_changed();
}

bool RunController::drain_tracks(int budget_ms) {
  if (!engine_) return false;
  QElapsedTimer t;
  t.start();
  bool any = false;
  std::vector<common::CollisionPoint> points;
  while (engine_->try_pop_trajectory(points)) {
    if (!tracks_->full()) {
      tracks_->add_history(points);
      any = true;
    }
    if (t.elapsed() > budget_ms) break;
  }
  return any;
}

void RunController::refresh_dataset(bool final_results) {
  if (!engine_) return;
  try {
    std::shared_ptr<output::Dataset> d;
    if (trim_ != nullptr) {
      d = std::make_shared<output::Dataset>(final_results ? trim::make_dataset(*trim_)
                                                          : trim::make_dataset(*trim_, trim_->snapshot(false)));
    } else {
      d = std::make_shared<output::Dataset>(
          native_dataset(engine_->snapshot_tallies(), QString::fromStdString(config_.trim.title)));
    }
    if (d->configuration_json.empty()) d->configuration_json = config::to_json(config_);
    dataset_ = std::move(d);
    emit dataset_changed();
  } catch (const std::exception &e) {
    // A snapshot can fail only transiently (e.g. while tallies are merged);
    // keep the previous one, and say so in the log (once per run).
    if (refresh_failures_++ == 0 || final_results) {
      AMCSET_LOG(warning) << "GUI: results view not refreshed: " << e.what();
    }
  } catch (...) {
    if (refresh_failures_++ == 0 || final_results) {
      AMCSET_LOG(warning) << "GUI: results view not refreshed (unknown exception)";
    }
  }
}

void RunController::poll() {
  if (sizing_) {
    if (done_.load()) finish_sizing();
    return;
  }
  if (!engine_) return;
  completed_ = engine_->completed_count();
  const double secs = static_cast<double>(clock_.elapsed()) / 1000.0;
  emit progress(completed_, total_, secs > 0 ? static_cast<double>(completed_) / secs : 0.0);
  if (trim_ != nullptr && config_.trim.auto_depth.warn_transmitted &&
      config_.trim.source == trim::TrimConfig::Source::beam) {
    const std::size_t n = trim_->transmitted_so_far();
    if (n > escaped_reported_) {
      escaped_reported_ = n;
      emit ions_escaped(n, completed_);
    }
  }
  if (drain_tracks(30)) emit tracks_changed();
  if (done_.load()) {
    finish();
    return;
  }
  if (clock_.elapsed() - last_snapshot_ms_ >= snapshot_interval_ms) {
    last_snapshot_ms_ = clock_.elapsed();
    refresh_dataset(false);
  }
}

void RunController::finish() {
  timer_.stop();
  if (thread_.joinable()) thread_.join();
  while (drain_tracks(1000)) {
  }
  emit tracks_changed();
  completed_ = engine_->completed_count();
  const bool cancelled = engine_->is_cancelled();
  const bool ok = !error_;
  QString message;
  if (!ok) {
    message = tr("The calculation failed: %1").arg(what_of(error_));
  } else {
    refresh_dataset(true);
    const QString dir = QString::fromStdString(config_.trim.output_directory);
    message = cancelled ? tr("Stopped after %1 of %2 particles.").arg(completed_).arg(total_)
                        : tr("Finished %1 particles in %2 s.")
                              .arg(completed_)
                              .arg(static_cast<double>(clock_.elapsed()) / 1000.0, 0, 'f', 1);
    if (!dir.isEmpty()) message += tr(" Output files in %1.").arg(QFileInfo(dir).absoluteFilePath());
  }
  AMCSET_LOG(info) << "GUI: " << message.toStdString();
  elapsed_offset_ += static_cast<double>(clock_.elapsed()) / 1000.0;
  if (ok && trim_ != nullptr) {
    previous_ = trim_->tallies();
    // Particles of the whole calculation so far (the next ion number).
    previous_ions_ = previous_.ions;
    previous_seconds_ = elapsed_offset_;
    has_previous_ = true;
  } else {
    has_previous_ = false;
  }
  set_state(State::idle);
  emit progress(completed_, total_, 0.0);
  emit finished(ok, message);
}

} // namespace gui
} // namespace amcset
