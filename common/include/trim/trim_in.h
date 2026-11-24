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

/*! \file trim_in.h
 *  \brief Reading and writing TRIM's batch input file TRIM.IN.
 *
 * TRIM.IN is the file TRIM reads in batch mode (TRIMAUTO = 1). Its layout is
 * documented in the SRIM manual ("TRIM.IN"); comment lines and value lines
 * alternate. AMCSET reads every field that has a meaning in SRIM mode and
 * writes TRIM.IN files that SRIM-2013 TRIM accepts, so the same scenario can
 * be run in both programs.
 */

#pragma once

#include "trim/trim.h"
#include "trim/trim_dat.h"

#include <filesystem>
#include <iosfwd>
#include <string>

namespace amcset {
namespace trim {

//! Extra TRIM.IN fields that have no TrimConfig equivalent.
struct TrimInExtras {
  int cascade_code = 1;          //!< 1..7 as in TRIM.IN
  int autosave = 0;              //!< Autosave interval
  int reminders = 0;
  int plot_type = 5;
  int stopping_version = 0;
  bool ranges_file = true;
  bool backscatter_file = true;
  bool transmit_file = true;
  bool sputter_file = true;
  int collisions_file = 0;
  int exyz_interval = 0;
};

/*!
 * \brief Parse a TRIM.IN file into a TrimConfig.
 *
 * Cascade types 4-7 set TrimConfig::source; their particles come from
 * TRIM.DAT, which read_trim_in() loads (parse_trim_in() leaves
 * TrimConfig::particles empty).
 *
 * Elements that TRIM lists once but that appear in several layers become one
 * TargetElement per layer sharing the same name, so their tallies are
 * combined into one column exactly as in TRIM.
 *
 * \param extras If not null, receives the fields without TrimConfig equivalent.
 * \throws std::runtime_error on malformed input (message names the line).
 */
TrimConfig parse_trim_in(std::istream &in, TrimInExtras *extras = nullptr);
//! parse_trim_in() from a file; for cascade types 4-7 also reads the
//! TRIM.DAT next to it. \throws std::runtime_error
TrimConfig read_trim_in(const std::filesystem::path &path, TrimInExtras *extras = nullptr);
//! TRIM.DAT (any letter case) in \p directory, or an empty path.
std::filesystem::path find_trim_dat(const std::filesystem::path &directory);

//! The TRIM.IN disk-file settings that match \p config (lists, COLLISON.txt,
//! EXYZ.txt); the other fields keep their defaults.
TrimInExtras extras_for(const TrimConfig &config);

//! Write \p config as a TRIM.IN that SRIM-2013 TRIM accepts.
void write_trim_in(std::ostream &out, const TrimConfig &config,
                   const TrimInExtras &extras = {});

} // namespace trim
} // namespace amcset
