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

/*! \file writers.h
 *  \brief CSV, JSON and Apache Parquet files from an output::Dataset.
 *
 * - CSV: one file per table, `csv/<table>.csv`, RFC 4180 quoting, headings
 *   "Label [unit]", optional "# ..." comment lines with the title,
 *   description and run metadata.
 * - JSON: one file, `results.json`, with the metadata, the configuration,
 *   the scalars (value, unit, dimension) and every table column-wise
 *   (column metadata plus one array of values per column).
 * - Parquet: one file per table, `parquet/<table>.parquet` (Parquet format 1, one
 *   row group, plain encoding, no compression, so any Parquet reader can
 *   open it without codecs). Real columns are DOUBLE, integer columns INT64,
 *   text columns UTF-8 strings; the units, labels and run metadata are in
 *   the file's key-value metadata ("amcset.units", "amcset.labels",
 *   "amcset.metadata", ...).
 */

#pragma once

#include "output/table.h"

#include <filesystem>
#include <iosfwd>
#include <string>
#include <vector>

namespace amcset {
namespace output {

//! Write one (already shaped) table as CSV.
void write_csv(std::ostream &out, const Table &table, const Dataset &context,
               const OutputSpec &spec);
//! Write the whole dataset as one JSON document (tables shaped by \p spec).
void write_json(std::ostream &out, const Dataset &data, const OutputSpec &spec);
//! Write one (already shaped) table as a Parquet file. \throws std::runtime_error
void write_parquet(const std::filesystem::path &path, const Table &table, const Dataset &context,
                   const OutputSpec &spec);

/*!
 * \brief Write \p data in every format of \p spec except Format::trim
 * (engines write their TRIM files themselves) to \p directory: CSV files
 * in `csv/`, Parquet files in `parquet/`, the JSON document as
 * `results.json`.
 * \return The files written.
 * \throws std::invalid_argument for an invalid spec, std::runtime_error on I/O errors.
 */
std::vector<std::filesystem::path> write_dataset(const std::filesystem::path &directory,
                                                 const Dataset &data, const OutputSpec &spec);

} // namespace output
} // namespace amcset
