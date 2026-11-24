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

/*! \file output_panel.h
 *  \brief Output folder, file formats, depth window, the optional TRIM files
 *  and the shape of the output tables.
 */

#pragma once

#include "output/table.h"

#include <QWidget>

#include <functional>
#include <memory>

class QLabel;

namespace amcset {
namespace gui {

class FormBinder;
class Scenario;

class OutputPanel : public QWidget {
  Q_OBJECT

public:
  /*!
   * \param sample Returns tables for the current target (the last results,
   * or empty tables), for the table settings dialog.
   */
  OutputPanel(Scenario *scenario, std::function<std::shared_ptr<const output::Dataset>()> sample,
              QWidget *parent = nullptr);

signals:
  void message(const QString &text);

private:
  void customize();

  Scenario *scenario_;
  FormBinder *binder_;
  std::function<std::shared_ptr<const output::Dataset>()> sample_;
  QLabel *summary_;
};

} // namespace gui
} // namespace amcset
