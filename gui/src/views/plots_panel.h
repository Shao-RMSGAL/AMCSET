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

/*! \file plots_panel.h
 *  \brief TRIM's result plots (ion and recoil distributions, damage, energy
 *  loss, sputtering, ...) and custom plots of any result table, live while a
 *  calculation runs.
 */

#pragma once

#include "output/table.h"

#include <QWidget>

#include <memory>

class QComboBox;
class QListWidget;
class QCheckBox;
class QLabel;

namespace amcset {
namespace gui {

class PlotWidget;

class PlotsPanel : public QWidget {
  Q_OBJECT

public:
  explicit PlotsPanel(QWidget *parent = nullptr);

  //! New results (live snapshot or final).
  void set_dataset(std::shared_ptr<const output::Dataset> data);
  //! Units chosen for the output files (used unless overridden here).
  void set_output_units(const output::OutputSpec &spec);
  //! Surface binding energy of the first element of layer 1 (sputtering plots).
  void set_surface_binding_energy(double eV);
  //! Save an image of every plot as plot_<id>.png in \p dir (unattended
  //! screenshots); the selection is restored afterwards.
  void save_plot_images(const QString &dir);

signals:
  void hover_text(const QString &text);

private:
  struct PlotDef {
    QString id;
    QString title;
    QString table;
    QString x;
    QStringList y; //!< column keys; "prefix*" matches several
    bool log_x = false;
    bool log_y = false;
    int style = 1; // PlotSeries::Style
    bool fill = false;
  };
  void rebuild_list();
  void refresh(bool keep_zoom);
  void update_custom_columns();
  output::OutputSpec unit_spec() const;
  const PlotDef *current() const;

  std::shared_ptr<const output::Dataset> data_;
  output::OutputSpec spec_;
  double sbe_ = 0.0;
  std::vector<PlotDef> defs_;
  PlotDef custom_;
  QListWidget *list_;
  PlotWidget *plot_;
  QComboBox *length_unit_;
  QComboBox *energy_unit_;
  QCheckBox *log_y_;
  QWidget *custom_box_;
  QComboBox *custom_table_;
  QComboBox *custom_x_;
  QListWidget *custom_y_;
  QLabel *note_;
  QString selected_id_;
};

} // namespace gui
} // namespace amcset
