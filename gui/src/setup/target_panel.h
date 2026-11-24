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

/*! \file target_panel.h
 *  \brief The target: layers and their elements.
 */

#pragma once

#include <QWidget>

class QLabel;
class QTableView;

namespace amcset {
namespace gui {

class Scenario;
class LayerTableModel;
class ElementTableModel;

class TargetPanel : public QWidget {
  Q_OBJECT

public:
  TargetPanel(Scenario *scenario, QWidget *parent = nullptr);

signals:
  void message(const QString &text);

private:
  int current_layer() const;
  void select_layer(int row);
  void add_layer();
  void add_compound();
  void remove_layer();
  void move_layer(int delta);
  void add_element();
  void remove_element();
  void update_summary();

  Scenario *scenario_;
  LayerTableModel *layers_;
  ElementTableModel *elements_;
  QTableView *layer_view_;
  QTableView *element_view_;
  QLabel *summary_;
  QLabel *element_title_;
};

} // namespace gui
} // namespace amcset
