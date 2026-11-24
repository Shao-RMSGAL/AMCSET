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

/*! \file quantity_text.h
 *  \brief Physical values typed by the user, with units (Boost.Units).
 */

#pragma once

#include "amcset_units.h"

#include <QString>

namespace amcset {
namespace gui {

//! "1.2 um" or "12000" (= default unit). \throws std::invalid_argument
common::length_quantity parse_length(const QString &text, const char *default_unit = "A");
common::energy_quantity parse_energy(const QString &text, const char *default_unit = "eV");
common::mass_density_quantity parse_density(const QString &text, const char *default_unit = "g/cm3");
common::angle_quantity parse_angle(const QString &text, const char *default_unit = "deg");
common::mass_quantity parse_mass(const QString &text, const char *default_unit = "amu");

QString format_length(common::length_quantity v);
QString format_energy(common::energy_quantity v);
QString format_density(common::mass_density_quantity v);
QString format_number(double v, int significant = 6);

//! Element by symbol ("Fe", case-insensitive) or atomic number ("26"). \throws std::invalid_argument
unsigned parse_element(const QString &text);
//! "Fe".
QString element_symbol(unsigned z);

} // namespace gui
} // namespace amcset
