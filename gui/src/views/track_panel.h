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

/*! \file track_panel.h
 *  \brief The particle-path view with its display options.
 */

#pragma once

#include <QWidget>

#include <vector>

class QComboBox;
class QCheckBox;
class QSlider;
class QSpinBox;

namespace amcset {
namespace gui {

class TrackStore;
class TrackView;

class TrackPanel : public QWidget {
  Q_OBJECT

public:
  TrackPanel(TrackStore *store, QWidget *parent = nullptr);
  TrackView *view() const noexcept { return view_; }
  void set_target(const std::vector<double> &boundaries, unsigned incident_z);
  void data_changed();
  //! Persisted display options.
  void save_settings() const;
  void restore_settings();

signals:
  void hover_text(const QString &text);
  void clear_requested();

private:
  TrackStore *store_;
  TrackView *view_;
  QComboBox *arrangement_;
  QComboBox *color_by_;
  QCheckBox *incident_;
  QCheckBox *recoils_;
  QSpinBox *generation_;
  QSlider *opacity_;
};

} // namespace gui
} // namespace amcset
