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

/*! \file stats_panel.h
 *  \brief The running totals TRIM shows in its calculation window: ranges,
 *  backscattering, transmission, sputtering, damage and the energy balance.
 */

#pragma once

#include "output/table.h"

#include <QWidget>

#include <memory>
#include <vector>

#include <QStringList>

class QTreeWidget;
class QTreeWidgetItem;

namespace amcset {
namespace gui {

class StatsPanel : public QWidget {
  Q_OBJECT

public:
  explicit StatsPanel(QWidget *parent = nullptr);
  void set_dataset(std::shared_ptr<const output::Dataset> data);
  void set_output_units(const output::OutputSpec &spec);
  //! Progress line: done / total, rate, remaining time.
  void set_progress(std::size_t done, std::size_t total, double per_second, double elapsed_s);

private:
  void refresh();
  QTreeWidgetItem *group(const QString &name);

  std::shared_ptr<const output::Dataset> data_;
  output::OutputSpec spec_;
  QTreeWidget *tree_;
  QTreeWidgetItem *progress_ = nullptr;
  QTreeWidgetItem *rate_ = nullptr;
  QTreeWidgetItem *remaining_ = nullptr;
  QStringList keys_;                       // rows shown, in order
  std::vector<QTreeWidgetItem *> values_;  // one per key
};

} // namespace gui
} // namespace amcset
