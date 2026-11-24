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

/*! \file periodic_table.h
 *  \brief Pick an element from the periodic table (Z = 1-92).
 */

#pragma once

#include <QDialog>

namespace amcset {
namespace gui {

class PeriodicTableDialog : public QDialog {
  Q_OBJECT

public:
  explicit PeriodicTableDialog(QWidget *parent = nullptr, unsigned current = 0);
  unsigned selected() const noexcept { return selected_; }
  //! Show the dialog; 0 if cancelled.
  static unsigned pick(QWidget *parent, unsigned current = 0);

private:
  unsigned selected_ = 0;
};

} // namespace gui
} // namespace amcset
