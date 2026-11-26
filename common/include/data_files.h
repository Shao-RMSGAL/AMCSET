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

#pragma once
//! \file data_files.h
//! \brief Data files that come with AMCSET (data/ in the source tree).

#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace amcset::data {

/*!
 * \brief Folders searched for AMCSET's data files, in order:
 * $AMCSET_DATA_DIR, data/ next to the executable and up to two folders above
 * it, share/amcset/ next to the executable's folder, and the data/ folder of
 * the source tree the program was built from.
 */
std::vector<std::filesystem::path> search_directories();

/*!
 * \brief The first existing file \p relative (e.g.
 * "srim2013/srim2013_stopping.amcset") in search_directories(); empty if
 * none has it.
 */
std::filesystem::path find(std::string_view relative);

//! Relative path of the SRIM-2013 stopping tables.
inline constexpr std::string_view srim_stopping_table = "srim2013/srim2013_stopping.amcset";

//! find(srim_stopping_table) as a string (empty if not found).
std::string bundled_srim_stopping_table();

//! Folder of the native engine's electron tables among the data files
//! (installed as share/amcset/electrons; common/assets in the source tree).
inline constexpr std::string_view electron_tables = "electrons";

/*!
 * \brief A table of the native engine's electron model (e.g.
 * "mott_scattering_parameters.csv"): electrons/<file> in
 * search_directories(), else common/assets/<file> of the source tree the
 * program was built from; empty if neither has it.
 */
std::filesystem::path find_electron_table(std::string_view file);

} // namespace amcset::data
