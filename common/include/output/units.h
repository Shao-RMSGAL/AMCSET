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

/*! \file units.h
 *  \brief Units of the values in output tables, and conversions between them.
 *
 * Every numeric column of an output table has a physical dimension and is
 * stored in one unit of it. Users choose the unit each column is written in;
 * the conversion factors come from the Boost.Units quantities of
 * amcset_units.h (each unit below is defined as a quantity of its dimension,
 * so a factor can only relate units of the same dimension), and a conversion
 * between different dimensions is rejected.
 */

#pragma once

#include <string>
#include <string_view>
#include <vector>

namespace amcset {
namespace output {

//! Physical dimension of a column.
enum class Dimension {
  dimensionless,     //!< Pure numbers (counts per ion, fractions, cosines, moments)
  length,            //!< Depths, positions, ranges
  inverse_length,    //!< Densities per ion and per unit depth, e.g. (atoms/cm3)/(atoms/cm2)
  energy,            //!< Energies
  energy_per_length, //!< Stopping powers, energy deposition per unit depth
  inverse_energy,    //!< Spectra per unit energy
  angle,             //!< Angles
  time,              //!< Times
};

//! Name of a dimension ("length", "energy_per_length", ...).
std::string to_string(Dimension d);
//! Inverse of to_string(). \throws std::invalid_argument
Dimension parse_dimension(std::string_view s);

//! One unit of a dimension.
struct Unit {
  std::string name;        //!< Canonical symbol, e.g. "A", "keV", "eV/A"
  Dimension dimension;
  double per_base = 1.0;   //!< Value of one base unit (SI) expressed in this unit
  std::vector<std::string> aliases;
};

/*!
 * \brief Look up a unit by name or alias (case-sensitive for SI prefixes;
 * "Angstrom", "angstrom", "Å" and "A" are the same unit).
 * \throws std::invalid_argument naming the known units of that dimension if
 * \p dimension is given, or all units otherwise.
 */
const Unit &find_unit(std::string_view name);
//! find_unit() and check the dimension. \throws std::invalid_argument
const Unit &find_unit(std::string_view name, Dimension dimension);
//! All units of a dimension, the default (TRIM's) one first.
std::vector<const Unit *> units_of(Dimension d);
//! The unit TRIM uses for a dimension (A, eV, eV/A, ...).
const Unit &trim_unit(Dimension d);
//! The SI base unit of a dimension (m, J, J/m, ...).
const Unit &base_unit(Dimension d);

/*!
 * \brief Factor that converts a value in \p from into \p to.
 * \throws std::invalid_argument if the dimensions differ.
 */
double conversion_factor(const Unit &from, const Unit &to);
//! conversion_factor() by unit names.
double conversion_factor(std::string_view from, std::string_view to);

} // namespace output
} // namespace amcset
