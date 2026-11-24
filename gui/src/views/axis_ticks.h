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

/*! \file axis_ticks.h
 *  \brief Tick positions and labels shared by the plots and path views.
 */

#pragma once

#include <QString>

#include <vector>

namespace amcset {
namespace gui {

//! About \p max_count "nice" values (1, 2, 5 x 10^n steps) within [lo, hi].
std::vector<double> nice_ticks(double lo, double hi, int max_count);
//! Decades (and 2, 5 subdivisions when few) within [lo, hi] (both > 0).
std::vector<double> log_ticks(double lo, double hi, int max_count);
//! Short label of a tick value ("0", "250", "1.5", "2E-05").
QString format_tick(double v);

//! A length unit for axes spanning up to \p max_abs Angstrom.
struct AxisUnit {
  QString name;
  double factor; //!< Angstrom per unit
};
AxisUnit length_axis_unit(double max_abs_angstrom);
//! "123.4 nm" for a length in Angstrom.
QString format_length_value(double angstrom);

} // namespace gui
} // namespace amcset
