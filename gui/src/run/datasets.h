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

/*! \file datasets.h
 *  \brief Result datasets the GUI shows that do not come from the SRIM-mode
 *  engine: the native engine's live histograms and saved results.json files.
 */

#pragma once

#include "output/table.h"

#include <QString>

namespace amcset {
namespace common {
struct Tallies;
}
namespace gui {

//! The native engine's tallies as tables (depth histograms per incident particle).
output::Dataset native_dataset(const common::Tallies &tallies, const QString &title);

/*!
 * \brief Read a results.json written by AMCSET (output::write_json).
 * \return false with \p error set if the file cannot be read.
 */
bool read_results_json(const QString &path, output::Dataset &out, QString *error);

} // namespace gui
} // namespace amcset
