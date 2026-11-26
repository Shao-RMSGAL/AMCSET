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

/*! \file crystal_panel.h
 *  \brief The crystal of a target layer: unit cell, orientation, thermal
 *  vibrations, and the settings of the crystal transport.
 *
 * Not TRIM (TRIM targets are amorphous): a layer's crystal is followed by
 * the crystal transport only. The panel edits `layers.N.crystal` of the
 * scenario (trim::TargetLayer::crystal) through Scenario::modify; the
 * layer's composition and density then follow from the unit cell
 * (trim::apply_crystal). The unit cell is shown in three dimensions
 * (CellView): clicking an atom edits its site, and the direction of view
 * can be made the direction of the beam.
 */

#pragma once

#include "crystal/crystal.h"
#include "run_config.h"

#include <QString>
#include <QWidget>

#include <functional>
#include <optional>
#include <string>
#include <vector>

class QCheckBox;
class QComboBox;
class QGroupBox;
class QLabel;
class QLineEdit;
class QPushButton;
class QSpinBox;
class QTableWidget;

namespace amcset {
namespace gui {

class CellView;
class FormBinder;
class Scenario;

class CrystalPanel : public QWidget {
  Q_OBJECT

public:
  CrystalPanel(Scenario *scenario, QWidget *parent = nullptr);
  //! The three-dimensional view (for screenshots).
  CellView *cell_view() const noexcept { return view_; }

signals:
  void message(const QString &text);

private:
  using Edit = std::function<void(crystal::Crystal &, trim::TargetLayer &, config::RunConfig &)>;

  int layer() const;
  const trim::TargetLayer *current_layer() const;
  const crystal::Crystal *current_crystal() const;
  trim::TargetLayer &layer_in(config::RunConfig &c) const;
  crystal::Crystal &crystal_in(config::RunConfig &c) const;
  //! Change the selected layer's crystal (undoable); a rejected change is reported.
  bool edit(const QString &what, const Edit &f);
  //! Give the selected layer the unit cell \p cell (making it a crystal).
  bool set_cell(const QString &what, const crystal::UnitCell &cell,
                const std::vector<crystal::Vibration> *vibrations, const std::string &layer_name,
                std::optional<common::mass_density_quantity> density);

  void refresh(const config::RunConfig &c);
  void refresh_layers(const config::RunConfig &c);
  void refresh_sites(const crystal::Crystal *crystal);
  void refresh_vibrations(const config::RunConfig &c, const trim::TargetLayer *layer);
  void refresh_directions(const config::RunConfig &c, const crystal::Crystal *crystal);

  void set_kind(int index);
  void choose_material(int index);
  //! The crystal database: search it and give the layer the chosen crystal.
  void choose_from_database(const QString &query);
  void build_structure();
  void import_cell();
  void export_cell();
  void add_site();
  void remove_site();
  void site_edited(int row, int column);
  void vibration_edited(int row, int column);
  void site_menu(int site, const QPoint &global_position);
  void select_site(int site);
  void beam_from_text();
  void beam_along_view();
  void orient_to_view();
  void look_from_text();
  void set_crystal_transport(bool on);

  Scenario *scenario_;
  FormBinder *binder_;
  QComboBox *layer_ = nullptr;
  QComboBox *kind_ = nullptr;
  QLabel *transport_note_ = nullptr;
  QLabel *match_ = nullptr; // an amorphous layer whose composition is a database crystal
  QCheckBox *follow_ = nullptr; // the crystal transport on/off (as in the Calculation tab)
  QComboBox *material_ = nullptr;
  QCheckBox *keep_shape_ = nullptr;
  QLabel *summary_ = nullptr;
  CellView *view_ = nullptr;
  QSpinBox *repeat_ = nullptr;
  QLineEdit *look_ = nullptr;
  QTableWidget *sites_ = nullptr;
  QLineEdit *beam_ = nullptr;
  QString beam_text_; // what the scenario says (the field shows it unless being edited)
  QLabel *directions_ = nullptr;
  QTableWidget *vibrations_ = nullptr;
  std::vector<QWidget *> needs_crystal_; // enabled only when the layer has a crystal
  std::vector<QWidget *> needs_layer_;
};

} // namespace gui
} // namespace amcset
