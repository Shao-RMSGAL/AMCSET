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

/*! \file target_models.h
 *  \brief Table models of the target layers and of one layer's elements.
 *
 * Both models read the Scenario and write every edit back through
 * Scenario::modify(), so edits are undoable and the other editors follow.
 * Physical values are shown with units and accept any unit the scenario
 * files accept ("1.2 um", "7.87 g/cm3", "25 eV"); a number without a unit
 * takes the column's default unit.
 */

#pragma once

#include <QAbstractTableModel>

namespace amcset {
namespace gui {

class Scenario;

class LayerTableModel : public QAbstractTableModel {
  Q_OBJECT

public:
  enum Column { name, thickness, density, gas, compound_correction, column_count };
  explicit LayerTableModel(Scenario *scenario, QObject *parent = nullptr);

  int rowCount(const QModelIndex &parent = {}) const override;
  int columnCount(const QModelIndex &parent = {}) const override;
  QVariant data(const QModelIndex &index, int role) const override;
  bool setData(const QModelIndex &index, const QVariant &value, int role) override;
  QVariant headerData(int section, Qt::Orientation o, int role) const override;
  Qt::ItemFlags flags(const QModelIndex &index) const override;

signals:
  void edit_failed(const QString &message);

private:
  void sync();
  Scenario *scenario_;
  int rows_ = 0;
};

class ElementTableModel : public QAbstractTableModel {
  Q_OBJECT

public:
  enum Column { element, stoichiometry, mass, ed, eb, es, column_count };
  explicit ElementTableModel(Scenario *scenario, QObject *parent = nullptr);

  //! The layer whose elements are shown (-1 = none).
  void set_layer(int layer);
  int layer() const noexcept { return layer_; }

  int rowCount(const QModelIndex &parent = {}) const override;
  int columnCount(const QModelIndex &parent = {}) const override;
  QVariant data(const QModelIndex &index, int role) const override;
  bool setData(const QModelIndex &index, const QVariant &value, int role) override;
  QVariant headerData(int section, Qt::Orientation o, int role) const override;
  Qt::ItemFlags flags(const QModelIndex &index) const override;

signals:
  void edit_failed(const QString &message);

private:
  void sync();
  Scenario *scenario_;
  int layer_ = -1;
  int rows_ = 0;
};

} // namespace gui
} // namespace amcset
