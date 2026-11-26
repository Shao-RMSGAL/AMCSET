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

/*! \file crystal_database.h
 *  \brief A searchable database of crystals in the built-in structures:
 *  elements and compounds with lattice constants (room temperature) and,
 *  where known, Debye temperatures, each value with its source.
 *
 * Gemmell's Tables IV and V (materials.h) come first and win for the
 * materials they list (so a crystal chosen there or here is the same); the
 * other rows come from the sources of database_sources() (docs/CRYSTAL_DATABASE.md
 * lists them). Only crystals whose structure is one of builtin_structures()
 * are included: graphite, alpha-Mn, beta-Sn, alpha-U, corundum, ... are not.
 */

#pragma once

#include "crystal/crystal.h"

#include <string_view>
#include <utility>
#include <vector>

namespace amcset {
namespace crystal {

//! One crystal of the database.
struct DatabaseCrystal {
  std::string_view name;       //!< formula or symbol: "GaAs", "Fe", "ZnS"
  std::string_view structure;  //!< a builtin_structures() key
  double a_A = 0.0;            //!< lattice constant a (Angstrom)
  double c_A = 0.0;            //!< c (Angstrom) for hcp and wurtzite; 0 otherwise
  //! The elements of the structure's sublattices, in its order (one each).
  std::vector<unsigned> sublattices;
  double debye_temperature_K = 0.0; //!< 0: not known here
  std::string_view lattice_source;  //!< a database_sources() key
  std::string_view debye_source;    //!< a database_sources() key, or empty
  std::string_view note;            //!< phase, conditions ("high-pressure phase", "ZrC0.97")
  //! The form chosen for its composition by match_composition() (false for
  //! a second, less common form of the same composition).
  bool preferred = true;
};

//! A source of database values.
struct DatabaseSource {
  std::string_view key;       //!< "ioffe", "gemmell1974", ...
  std::string_view citation;  //!< the full reference
};

//! Every crystal of the database (Gemmell's rows first).
const std::vector<DatabaseCrystal> &crystal_database();
//! The sources the rows name.
const std::vector<DatabaseSource> &database_sources();
//! The citation of source \p key (empty if unknown).
std::string_view database_citation(std::string_view key);

/*!
 * \brief The rows that match \p query: every word of it (case-insensitive)
 * must occur in the row's name, its elements' symbols or names, the
 * structure's key or title, or its note ("Fe", "iron", "zinc blende",
 * "nitride rocksalt", "Ga As"). An empty query gives every row.
 */
std::vector<const DatabaseCrystal *> search_crystal_database(std::string_view query);

/*!
 * \brief The row a scenario names ("material": "GaAs"): the preferred row
 * of that name (case-insensitive; "ZnS" is sphalerite), else the only row
 * search_crystal_database() finds ("ZnS wurtzite").
 * \throws std::invalid_argument naming the candidates if there are none or several.
 */
const DatabaseCrystal &find_database_crystal(std::string_view name);

//! Atoms of each element per unit cell of a row's structure (Z -> atoms).
std::vector<std::pair<unsigned, double>> database_composition(const DatabaseCrystal &row);

/*!
 * \brief The preferred row whose composition is \p composition (Z ->
 * amount, in any units): the same elements, every atomic fraction within
 * \p tolerance. Null if none (an alloy, an amorphous oxide, ...).
 */
const DatabaseCrystal *match_composition(const std::vector<std::pair<unsigned, double>> &composition,
                                         double tolerance = 0.01);

/*!
 * \brief The crystal of a row: its cell at the tabulated constants, a
 * (0 0 1) surface, 293 K, and the Debye temperature for every element if
 * the row has one. The cell's name is the row's name and structure.
 */
Crystal make_database_crystal(const DatabaseCrystal &row);

} // namespace crystal
} // namespace amcset
