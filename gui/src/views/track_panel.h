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

#include "trim/trim.h"

#include <QWidget>

#include <atomic>
#include <memory>
#include <vector>

class QComboBox;
class QCheckBox;
class QSlider;
class QSpinBox;
class QThread;
class QTimer;

namespace amcset {
namespace gui {

class TrackStore;
class TrackView;
struct SampleScene;

class TrackPanel : public QWidget {
  Q_OBJECT

public:
  TrackPanel(TrackStore *store, QWidget *parent = nullptr);
  ~TrackPanel() override;
  TrackView *view() const noexcept { return view_; }
  void set_target(const std::vector<double> &boundaries, unsigned incident_z);
  /*!
   * The scenario, for the whole-beam-spot view: its spot and its target's
   * layers and grains (drawn off the GUI thread). \p crystal_transport: the
   * runs follow the crystal transport, the only one with entry points.
   */
  void set_scenario(const trim::TrimConfig &c, bool crystal_transport);
  void data_changed();
  //! Persisted display options.
  void save_settings() const;
  void restore_settings();
  //! For the screenshots: the whole-beam-spot view, the arrangement (index
  //! of the combo box: 0 four views, 1-4 one of them) and an image of the views.
  bool whole_spot_available() const;
  bool whole_spot() const;
  void set_whole_spot(bool on);
  int arrangement_index() const;
  void set_arrangement_index(int i);
  QImage grab_image() const;

signals:
  void hover_text(const QString &text);
  void clear_requested();

private:
  void build_scene();
  void scene_built(quint64 generation, std::shared_ptr<const SampleScene> scene);

  TrackStore *store_;
  TrackView *view_;
  QComboBox *arrangement_;
  QComboBox *color_by_;
  QCheckBox *incident_;
  QCheckBox *recoils_;
  QSpinBox *generation_;
  QSlider *opacity_;
  QCheckBox *spot_;
  // The sample of the whole-beam-spot view, built on a worker thread.
  trim::TrimConfig scenario_;
  bool have_scenario_ = false;
  bool scene_stale_ = true;
  QTimer *scene_timer_;
  QThread *worker_ = nullptr;
  std::shared_ptr<std::atomic<bool>> cancel_;
  quint64 scene_generation_ = 0;
};

} // namespace gui
} // namespace amcset
