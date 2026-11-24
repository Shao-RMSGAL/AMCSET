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

/*! \file tables_panel.h
 *  \brief Browse the result tables, as they are or as the output settings
 *  shape them (units, columns, rows, number format), and export them.
 */

#pragma once

#include "output/table.h"

#include <QAbstractTableModel>
#include <QWidget>

#include <memory>

class QCheckBox;
class QComboBox;
class QLabel;
class QTableView;

namespace amcset {
namespace gui {

//! Read-only model of one output table.
class TableModel : public QAbstractTableModel {
public:
  using QAbstractTableModel::QAbstractTableModel;
  void set_table(output::Table table, output::OutputSpec spec);
  int rowCount(const QModelIndex &parent = {}) const override;
  int columnCount(const QModelIndex &parent = {}) const override;
  QVariant data(const QModelIndex &index, int role) const override;
  QVariant headerData(int section, Qt::Orientation o, int role) const override;
  const output::Table &table() const noexcept { return table_; }

private:
  output::Table table_;
  output::OutputSpec spec_;
  std::vector<output::NumberFormat> formats_;
};

class TablesPanel : public QWidget {
  Q_OBJECT

public:
  explicit TablesPanel(QWidget *parent = nullptr);
  void set_dataset(std::shared_ptr<const output::Dataset> data);
  void set_output_spec(const output::OutputSpec &spec);

private:
  void refresh();
  void export_table(output::Format f);

  std::shared_ptr<const output::Dataset> data_;
  output::OutputSpec spec_;
  QComboBox *tables_;
  QCheckBox *shaped_;
  QLabel *info_;
  QTableView *view_;
  TableModel *model_;
};

} // namespace gui
} // namespace amcset
