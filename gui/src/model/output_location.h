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

/*! \file output_location.h
 *  \brief Where the GUI writes results: never relative to its working
 *  directory (an installed GUI starts in a folder the user cannot write).
 */

#pragma once

#include <QString>

#include <string>

namespace amcset {
namespace gui {

/*!
 * \brief The folder that holds the results of new scenarios: the user's
 * choice (File, Output folder for new scenarios; QSettings
 * "paths/output_root"), else AMCSET in the user's documents folder
 * (Windows: Documents\\AMCSET, macOS: ~/Documents/AMCSET, Linux: the XDG
 * documents folder, usually ~/Documents/AMCSET). Absolute.
 */
QString output_root();
//! Remember \p dir as output_root() (empty: back to the default).
void set_output_root(const QString &dir);
//! The default of output_root() on this system.
QString default_output_root();
//! The output folder of a new scenario: output_root()/output.
std::string new_scenario_output_directory();
/*!
 * \brief A scenario's output folder made absolute: a relative one is
 * relative to the scenario file's folder (\p scenario_path, if saved), else
 * to output_root(); an empty one stays empty (no files are written).
 */
std::string absolute_output_directory(const std::string &dir, const QString &scenario_path);

} // namespace gui
} // namespace amcset
