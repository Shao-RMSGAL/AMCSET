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

#include "model/scenario.h"

#include "app/guard.h"
#include "diagnostics.h"
#include "trim/trim_dat.h"
#include "trim/trim_in.h"

#include <QDir>
#include <QFileInfo>
#include <QUndoCommand>

#include <filesystem>
#include <fstream>

namespace amcset {
namespace gui {

namespace {

class EditCommand : public QUndoCommand {
public:
  EditCommand(Scenario *s, config::RunConfig before, config::RunConfig after, const QString &text,
              int id)
      : QUndoCommand(text), s_(s), before_(std::move(before)), after_(std::move(after)), id_(id) {}
  void undo() override { s_->set_config_silently(before_); }
  void redo() override { s_->set_config_silently(after_); }
  int id() const override { return id_; }
  bool mergeWith(const QUndoCommand *other) override {
    const auto *o = static_cast<const EditCommand *>(other);
    if (o->id_ != id_ || id_ < 0) return false;
    after_ = o->after_;
    return true;
  }

private:
  Scenario *s_;
  config::RunConfig before_;
  config::RunConfig after_;
  int id_;
};

bool is_trim_in(const QString &path) {
  const QString name = QFileInfo(path).fileName().toUpper();
  return name == QLatin1String("TRIM.IN") || name.endsWith(QLatin1String(".IN"));
}

} // namespace

Scenario::Scenario(QObject *parent) : QObject(parent) {
  config_ = config::example_config(config::Mode::srim);
  config_.trim.trajectory_ions = 200;
  remember_transport();
  undo_.setClean();
}

void Scenario::remember_transport() {
  if (config_.trim.transport != trim::TrimConfig::Transport::crystal) amorphous_transport_ = config_.trim.transport;
}

QString Scenario::set_crystal_transport(bool on) {
  using Transport = trim::TrimConfig::Transport;
  if (on == (config_.trim.transport == Transport::crystal)) return {};
  const Transport back = amorphous_transport_;
  return modify(on ? tr("Crystal transport on") : tr("Crystal transport off"),
                [on, back](config::RunConfig &c) { c.trim.transport = on ? Transport::crystal : back; });
}

QString Scenario::modify(const QString &what, const std::function<void(config::RunConfig &)> &edit,
                         int merge_id) {
  try {
    config::RunConfig next = config_;
    edit(next);
    // Every edit must leave a scenario that can run (or is a step towards
    // one: a particle source may still lack its list). An out-of-range or
    // nonsense value is refused here, with the message saying what is
    // allowed, instead of failing later at Run.
    next.trim.validate(false);
    undo_.push(new EditCommand(this, config_, next, what, merge_id));
  } catch (...) {
    const QString err = exception_text(std::current_exception());
    AMCSET_LOG(info) << "GUI: edit '" << what.toStdString() << "' refused: " << err.toStdString();
    return err;
  }
  return {};
}

QString Scenario::problem() const {
  try {
    config_.trim.validate(true);
  } catch (...) {
    return exception_text(std::current_exception());
  }
  return {};
}

void Scenario::set_config_silently(const config::RunConfig &c) {
  config_ = c;
  remember_transport();
  emit changed();
}

void Scenario::reset(const config::RunConfig &c, const QString &path) {
  undo_.clear();
  config_ = c;
  remember_transport();
  path_ = path;
  undo_.setClean();
  emit changed();
  emit file_path_changed(path_);
}

QString Scenario::load(const QString &path) {
  try {
    config::RunConfig c;
    if (is_trim_in(path)) {
      c.mode = config::Mode::srim;
      c.trim = trim::read_trim_in(std::filesystem::path(path.toStdWString()));
      if (c.trim.output_directory.empty()) {
        c.trim.output_directory =
            QDir(QFileInfo(path).absolutePath()).filePath(QStringLiteral("AMCSET Outputs")).toStdString();
      }
      if (c.trim.trajectory_ions == 0) c.trim.trajectory_ions = 200;
      reset(c, {});
    } else {
      c = config::load_config(std::filesystem::path(path.toStdWString()));
      reset(c, path);
    }
  } catch (...) {
    return exception_text(std::current_exception());
  }
  return {};
}

QString Scenario::save(const QString &path) {
  try {
    config::save_config(config_, std::filesystem::path(path.toStdWString()));
  } catch (...) {
    return exception_text(std::current_exception());
  }
  path_ = path;
  undo_.setClean();
  emit file_path_changed(path_);
  return {};
}

QString Scenario::export_trim(const QString &directory) const {
  try {
    const std::filesystem::path dir(directory.toStdWString());
    std::filesystem::create_directories(dir);
    {
      std::ofstream o(dir / "TRIM.IN", std::ios::binary);
      if (!o) return tr("Cannot write %1").arg(QDir(directory).filePath("TRIM.IN"));
      trim::write_trim_in(o, config_.trim, trim::extras_for(config_.trim));
    }
    if (config_.trim.source != trim::TrimConfig::Source::beam) {
      std::ofstream o(dir / "TRIM.DAT", std::ios::binary);
      if (!o) return tr("Cannot write %1").arg(QDir(directory).filePath("TRIM.DAT"));
      trim::ParticleList list;
      list.title = config_.trim.particles_title;
      list.particles = config_.trim.particles;
      trim::write_trim_dat(o, list);
      if (!o) return tr("Error writing %1").arg(QDir(directory).filePath("TRIM.DAT"));
    }
  } catch (...) {
    return exception_text(std::current_exception());
  }
  return {};
}

} // namespace gui
} // namespace amcset
