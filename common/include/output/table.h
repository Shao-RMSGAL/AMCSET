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

/*! \file table.h
 *  \brief Results as tables of typed columns, and how users shape them.
 *
 * Engines describe their results as an output::Dataset: scalar results plus
 * tables whose numeric columns carry a physical dimension and a unit. The
 * writers (writers.h) turn a Dataset into CSV, JSON or Parquet files, and the
 * GUI plots its tables, so every result is available in every format.
 *
 * An OutputSpec says what to write and how: which formats, which tables and
 * columns, which rows (value ranges, zero rows, decimation), the unit of
 * every column (per dimension or per column) and the number format (general,
 * fixed or scientific notation, digits). apply() produces the shaped copy of
 * a table that the writers then print.
 */

#pragma once

#include "output/units.h"

#include <cstddef>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <vector>

namespace amcset {
namespace output {

//! Kind of values a column holds.
enum class ColumnType { real, integer, text };

//! One column of a table.
struct Column {
  std::string key;          //!< Identifier (letters, digits, '_'), unique within the table
  std::string label;        //!< Human-readable heading
  ColumnType type = ColumnType::real;
  Dimension dimension = Dimension::dimensionless;
  std::string unit = "1";   //!< Unit of #values (a name from units.h)
  std::string description;  //!< Optional longer explanation
  std::vector<double> values;    //!< real and integer columns
  std::vector<std::string> text; //!< text columns

  std::size_t size() const { return type == ColumnType::text ? text.size() : values.size(); }

  static Column real(std::string key, std::string label, Dimension d, std::string unit,
                     std::vector<double> values = {}, std::string description = {});
  static Column integer(std::string key, std::string label, std::vector<double> values = {},
                        std::string description = {});
  static Column strings(std::string key, std::string label, std::vector<std::string> text = {},
                        std::string description = {});
};

//! A table: columns of equal length.
struct Table {
  std::string key;          //!< Identifier, also the default file name (e.g. "range")
  std::string title;        //!< Heading
  std::string description;  //!< What the table contains and how it is normalised
  std::vector<Column> columns;

  std::size_t rows() const { return columns.empty() ? 0 : columns.front().size(); }
  //! Column by key, or nullptr.
  const Column *find(std::string_view key) const;
  Column *find(std::string_view key);
};

//! One scalar result.
struct Scalar {
  std::string key;
  std::string label;
  Dimension dimension = Dimension::dimensionless;
  std::string unit = "1";
  double value = 0.0;
  std::string description;
};

//! Everything a run produced.
struct Dataset {
  std::string title;
  //! Free-form facts about the run (program, version, ion, target, seed, ...).
  std::vector<std::pair<std::string, std::string>> metadata;
  //! The run configuration as JSON text (may be empty).
  std::string configuration_json;
  std::vector<Scalar> scalars;
  std::vector<Table> tables;

  const Table *find(std::string_view key) const;
};

// ---------------------------------------------------------------------------
// Shaping the output
// ---------------------------------------------------------------------------

//! How numbers are printed (CSV, TRIM-style text and JSON).
struct NumberFormat {
  enum class Notation {
    shortest,   //!< Shortest text that reads back exactly
    general,    //!< printf %g with `digits` significant digits (the default)
    fixed,      //!< printf %f with `digits` decimals
    scientific, //!< printf %E with `digits` decimals
    engineering //!< Scientific with the exponent a multiple of 3
  };
  //! Default: 12 significant digits, enough for any tally and free of the
  //! last-digit noise that unit conversions leave in "shortest" output.
  Notation notation = Notation::general;
  int digits = 12;
  bool uppercase = true;   //!< "E" rather than "e" in exponents

  bool operator==(const NumberFormat &) const = default;
};
std::string to_string(NumberFormat::Notation n);
//! "shortest", "general", "fixed", "scientific", "engineering". \throws std::invalid_argument
NumberFormat::Notation parse_notation(std::string_view s);
//! Print \p v with \p f ("nan", "inf", "-inf" for non-finite values).
std::string format_number(double v, const NumberFormat &f);
/*!
 * \brief format_number() into \p p, which must have room for
 * text::max_field (output/text.h) characters; returns the end. No
 * allocation (the writers call it for every cell).
 */
char *format_number_to(char *p, double v, const NumberFormat &f);

//! Keep only rows whose value in a column lies within [min, max].
struct RowFilter {
  std::string column;
  std::optional<double> min; //!< In min_unit, or the column's output unit if empty
  std::optional<double> max; //!< In max_unit, or the column's output unit if empty
  std::string min_unit;
  std::string max_unit;
};

//! Settings of one column.
struct ColumnSpec {
  std::string unit;                    //!< Output unit (empty = default)
  std::optional<NumberFormat> format;  //!< Number format (empty = table / global default)
  std::string label;                   //!< Replacement heading (empty = keep)
};

//! Settings of one table.
struct TableSpec {
  bool include = true;
  //! Columns to write, in this order (empty = all columns).
  std::vector<std::string> columns;
  //! Columns to leave out.
  std::vector<std::string> exclude;
  std::map<std::string, ColumnSpec> column_specs;
  std::vector<RowFilter> filters;
  //! Drop rows whose numeric columns (other than the first) are all zero.
  bool drop_zero_rows = false;
  //! Keep every n-th row (1 = all).
  std::size_t every = 1;
  //! Keep at most this many rows (0 = all), after the other filters.
  std::size_t max_rows = 0;
  std::optional<NumberFormat> format;
  //! File name (without extension) instead of the table key.
  std::string file_name;
};

//! The file formats the writers produce.
enum class Format {
  trim, //!< TRIM's own text files (RANGE.txt, VACANCY.txt, ...), where the engine has them
  csv,  //!< One CSV file per table, units in the headings
  json, //!< One JSON file with the metadata, scalars and all tables
  parquet //!< One Apache Parquet file per table
};
std::string to_string(Format f);
//! "trim", "txt", "csv", "json", "parquet". \throws std::invalid_argument
Format parse_format(std::string_view s);

//! What to write and how.
struct OutputSpec {
  std::set<Format> formats{Format::trim, Format::csv, Format::json};
  //! Output unit per dimension (absent = each column's own unit, which is TRIM's).
  std::map<Dimension, std::string> units;
  NumberFormat format;                    //!< Default number format
  std::map<std::string, TableSpec> tables; //!< By table key; "*" applies to all others
  char csv_delimiter = ',';
  //! Put "# ..." lines with the title, description and metadata at the top of CSV files.
  bool csv_comments = true;
  bool json_pretty = true;
  //! Parquet column names include the unit, e.g. "depth [A]".
  bool parquet_units_in_names = false;

  //! The spec of table \p key ("*" or a default spec if absent).
  const TableSpec &table(std::string_view key) const;
  /*!
   * \brief Check every unit, column and format name against \p data.
   * \throws std::invalid_argument naming the offending setting.
   */
  void validate(const Dataset &data) const;
};

/*!
 * \brief The copy of \p table that \p spec asks for: selected and reordered
 * columns, filtered rows, values converted to the output units (Column::unit
 * is updated), relabelled columns.
 * \throws std::invalid_argument for an unknown column or unit.
 */
Table apply(const Table &table, const OutputSpec &spec);
//! The number format of a column of a shaped table.
NumberFormat column_format(const Table &table, const Column &column, const OutputSpec &spec);
//! Scalars converted to the spec's per-dimension units.
std::vector<Scalar> apply(const std::vector<Scalar> &scalars, const OutputSpec &spec);

/*!
 * \brief Read an OutputSpec from JSON text, e.g.
 * \code{.json}
 * { "formats": ["csv", "json", "parquet"],
 *   "units": { "length": "nm", "energy": "keV" },
 *   "number_format": { "notation": "scientific", "digits": 4 },
 *   "tables": {
 *     "range": { "columns": ["depth", "ions"], "filters": [{"column": "depth", "max": 500}],
 *                "units": {"depth": "nm"}, "drop_zero_rows": true },
 *     "exits": { "include": false } } }
 * \endcode
 * \throws std::invalid_argument
 */
OutputSpec parse_output_spec(std::string_view json_text);
//! The JSON form of a spec (parse_output_spec reads it back).
std::string to_json(const OutputSpec &spec);

} // namespace output
} // namespace amcset
