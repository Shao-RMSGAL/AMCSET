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

/*! \file stream.h
 *  \brief A table written to disk piece by piece (CSV or Parquet), for
 *  outputs too large to hold in memory: the trajectory and MD event files.
 *
 * The columns are fixed when the file is opened; append() takes rows in
 * Tables with exactly those columns, and finish() completes the file. CSV
 * files have one header line of the column keys (no comment lines, so any
 * CSV reader takes them; the keys carry the units, e.g. "x_A") and numbers
 * in the shortest form that reads back exactly; Parquet
 * files are written as one row group per ~256 k rows (PLAIN encoding, no
 * compression, Parquet format 1; read by pandas/pyarrow/polars), with the
 * columns' units and descriptions in the key-value metadata as
 * output::write_parquet does.
 */

#pragma once

#include "output/table.h"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace amcset {
namespace output {

class StreamTableWriter {
public:
  enum class Format { csv, parquet };
  /*!
   * \brief Open \p path (created; parent folders too) with the columns of
   * \p schema (their values are ignored). The format follows the extension
   * (.csv or .parquet). \p metadata goes to CSV comment lines / Parquet
   * key-value metadata.
   * \throws std::invalid_argument (extension), std::runtime_error (file)
   */
  StreamTableWriter(const std::filesystem::path &path, const Table &schema,
                    std::vector<std::pair<std::string, std::string>> metadata = {});
  //! An unfinished file is closed as it is (a Parquet file without its
  //! footer cannot be read); finish() reports errors.
  ~StreamTableWriter() = default;
  StreamTableWriter(const StreamTableWriter &) = delete;
  StreamTableWriter &operator=(const StreamTableWriter &) = delete;

  //! Append the rows of \p rows (same column keys and types, in order).
  //! \throws std::invalid_argument (columns), std::runtime_error (writing)
  void append(const Table &rows);
  //! Write what is buffered and complete the file (Parquet footer).
  void finish();
  std::uint64_t rows_written() const noexcept { return rows_; }
  const std::filesystem::path &path() const noexcept { return path_; }
  Format format() const noexcept { return format_; }

  //! The Parquet writer's bookkeeping: where every row group's column chunks are.
  struct ParquetState {
    struct Chunk {
      std::int64_t offset = 0;
      std::int64_t size = 0;
    };
    struct RowGroup {
      std::int64_t rows = 0;
      std::int64_t bytes = 0;
      std::vector<Chunk> chunks;
    };
    std::int64_t offset = 0; //!< Bytes written so far
    std::vector<RowGroup> groups;
    std::vector<std::string> names;
  };

private:
  void parquet_open();
  void flush_parquet();
  void parquet_close();

  std::filesystem::path path_;
  Format format_ = Format::csv;
  Table schema_;
  std::vector<std::pair<std::string, std::string>> metadata_;
  std::ofstream out_;
  std::uint64_t rows_ = 0;
  bool finished_ = false;
  Table pending_; // Parquet rows not yet in a row group
  std::unique_ptr<ParquetState> parquet_;
};

} // namespace output
} // namespace amcset
