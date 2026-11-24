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

/*! \file calculation_panel.h
 *  \brief Engine, damage calculation, number of ions, physics and the
 *  options of the native engine.
 */

#pragma once

#include <QWidget>

#include <vector>

class QGroupBox;
class QLabel;

namespace amcset {
namespace gui {

class FormBinder;
class Scenario;

class CalculationPanel : public QWidget {
  Q_OBJECT

public:
  CalculationPanel(Scenario *scenario, QWidget *parent = nullptr);

signals:
  void message(const QString &text);

private:
  Scenario *scenario_;
  FormBinder *binder_;
  QGroupBox *srim_;
  QGroupBox *bugs_;
  QLabel *bugs_note_ = nullptr;        // shown with the crystal transport
  std::vector<QWidget *> bug_widgets_; // the switches and buttons: usable with trim2013 only
  QGroupBox *native_;
};

} // namespace gui
} // namespace amcset
