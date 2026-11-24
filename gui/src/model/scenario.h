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

/*! \file scenario.h
 *  \brief The scenario being edited: one config::RunConfig with undo.
 *
 * Every editor reads the scenario from here and changes it through
 * modify(), which records an undoable step and notifies the other editors.
 * Files: scenario JSON (native format) and TRIM.IN (+ TRIM.DAT) import and
 * export.
 */

#pragma once

#include "run_config.h"

#include <QObject>
#include <QString>
#include <QUndoStack>

#include <functional>

namespace amcset {
namespace gui {

class Scenario : public QObject {
  Q_OBJECT

public:
  explicit Scenario(QObject *parent = nullptr);

  const config::RunConfig &config() const noexcept { return config_; }

  /*!
   * \brief Change the scenario. \p edit receives a copy; if it throws, or
   * the result fails TrimConfig::validate (physical limits, nonsense
   * values), the change is dropped and the message returned. Similar
   * consecutive edits with the same \p merge_id (e.g. typing in one field)
   * merge into one undo step.
   * \return An empty string on success, otherwise the error message.
   */
  QString modify(const QString &what, const std::function<void(config::RunConfig &)> &edit,
                 int merge_id = -1);
  //! Replace everything (new/opened file); clears the undo history.
  void reset(const config::RunConfig &config, const QString &path = {});

  QUndoStack *undo_stack() noexcept { return &undo_; }
  QString file_path() const { return path_; }
  bool is_modified() const { return !undo_.isClean(); }
  void set_clean() { undo_.setClean(); }

  //! Open a scenario JSON or a TRIM.IN (with its TRIM.DAT). \return error or empty
  QString load(const QString &path);
  //! Why the current scenario cannot run (TrimConfig::validate), or empty.
  //! A file with such a problem still opens, so that it can be fixed.
  QString problem() const;
  //! Save as scenario JSON. \return error or empty
  QString save(const QString &path);
  //! Write TRIM.IN (and TRIM.DAT for particle sources) to \p directory.
  QString export_trim(const QString &directory) const;

  //! Apply a value without an undo step (used by undo/redo).
  void set_config_silently(const config::RunConfig &c);

  /*!
   * \brief Switch the crystal transport on or off (one undoable step). On
   * also switches the automatic target depth off, which that transport does
   * not offer (\p note then says so); off returns to the transport used
   * before it (TRIM-2013 at first).
   * \return An empty string on success, otherwise the error message.
   */
  QString set_crystal_transport(bool on, QString *note = nullptr);
  //! The transport in use when the crystal transport is off: the last other one chosen.
  trim::TrimConfig::Transport amorphous_transport() const noexcept { return amorphous_transport_; }

signals:
  void changed();
  void file_path_changed(const QString &path);

private:
  void remember_transport();

  config::RunConfig config_;
  QString path_;
  QUndoStack undo_;
  trim::TrimConfig::Transport amorphous_transport_ = trim::TrimConfig::Transport::trim2013;
};

} // namespace gui
} // namespace amcset
