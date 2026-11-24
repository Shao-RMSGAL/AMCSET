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

/*! \file output_spec_dialog.h
 *  \brief Choose what the output files contain and how: units, number
 *  format, file options, and per table the columns (order, units, labels,
 *  formats), the rows (value ranges, zero rows, every n-th, at most n) and
 *  the file name, with a live preview.
 */

#pragma once

#include "output/table.h"

#include <QDialog>

#include <map>
#include <memory>

class QCheckBox;
class QComboBox;
class QLabel;
class QLineEdit;
class QListWidget;
class QSpinBox;
class QTableView;
class QTableWidget;

namespace amcset {
namespace gui {

class TableModel;

class OutputSpecDialog : public QDialog {
  Q_OBJECT

public:
  /*!
   * \param sample Tables of the current target (for the column names and the
   * preview); may come from the last calculation or be empty.
   */
  OutputSpecDialog(const output::OutputSpec &spec, std::shared_ptr<const output::Dataset> sample,
                   QWidget *parent = nullptr);
  output::OutputSpec spec() const { return spec_; }

private:
  void load_general();
  void store_general();
  void select_table(int row);
  void load_table();
  void store_table();
  void select_column(int row);
  void store_column();
  void update_preview();
  output::TableSpec &table_spec();
  const output::Table *current_table() const;

  output::OutputSpec spec_;
  std::shared_ptr<const output::Dataset> sample_;
  bool loading_ = false;
  // General.
  std::map<output::Dimension, QComboBox *> units_;
  QComboBox *notation_;
  QSpinBox *digits_;
  QCheckBox *uppercase_;
  QComboBox *delimiter_;
  QCheckBox *comments_;
  QCheckBox *pretty_;
  QCheckBox *pq_units_;
  // Tables.
  QListWidget *tables_;
  QCheckBox *include_;
  QListWidget *columns_;
  QComboBox *col_unit_;
  QLineEdit *col_label_;
  QComboBox *col_notation_;
  QSpinBox *col_digits_;
  QTableWidget *filters_;
  QCheckBox *drop_zero_;
  QSpinBox *every_;
  QSpinBox *max_rows_;
  QLineEdit *file_name_;
  QComboBox *table_notation_;
  QSpinBox *table_digits_;
  QLabel *status_;
  QTableView *preview_;
  TableModel *preview_model_;
};

} // namespace gui
} // namespace amcset
