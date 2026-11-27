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

/*! \file grains_box.h
 *  \brief The Crystal tab's grains: a layer's crystal as a single crystal or
 *  a polycrystal (generated grains, a grain file, a Neper tessellation, an
 *  EBSD map), with a map of its surface and a summary of its grains.
 *
 * Edits `layers.N.crystal.polycrystal` (crystal::Polycrystal) through the
 * FormBinder; the preview and the summary (crystal::summarize_grains) are
 * made off the GUI thread from the grains the transport would build
 * (crystal::make_grain_map).
 */

#pragma once

#include "crystal/crystal.h"
#include "run_config.h"

#include <QGroupBox>
#include <QImage>

#include <atomic>
#include <functional>
#include <memory>
#include <optional>

class QCheckBox;
class QComboBox;
class QFormLayout;
class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;
class QThread;
class QTimer;

namespace amcset {
namespace gui {

class FormBinder;
class Scenario;

class GrainsBox : public QGroupBox {
  Q_OBJECT

public:
  //! \p layer gives the selected layer (-1: none).
  GrainsBox(Scenario *scenario, FormBinder *binder, std::function<int()> layer, QWidget *parent = nullptr);
  ~GrainsBox() override;

signals:
  void message(const QString &text);

private:
  struct Preview {
    QImage map;
    QString summary;
  };
  void refresh(const config::RunConfig &c);
  void choose_source(int index);
  void choose_file();
  void start_preview();
  void preview_done(quint64 generation, const Preview &p);
  const crystal::Polycrystal *current(const config::RunConfig &c) const;
  crystal::Polycrystal &edited(config::RunConfig &c) const;

  Scenario *scenario_;
  FormBinder *binder_;
  std::function<int()> layer_;
  QFormLayout *form_;
  QComboBox *source_;
  QLineEdit *grain_size_, *spread_, *period_, *boundary_, *registry_, *scale_, *misorientation_;
  QCheckBox *columnar_;
  QComboBox *texture_;
  QSpinBox *seed_, *phase_;
  QLabel *file_;
  QPushButton *choose_;
  QLabel *map_;
  QLabel *summary_;
  // The preview, made on a worker thread when the grains change.
  std::optional<crystal::Crystal> previewed_;
  double previewed_width_ = -1.0, previewed_tilt_ = -1.0;
  QTimer *timer_;
  QThread *worker_ = nullptr;
  std::shared_ptr<std::atomic<bool>> cancel_;
  quint64 generation_ = 0;
};

} // namespace gui
} // namespace amcset
