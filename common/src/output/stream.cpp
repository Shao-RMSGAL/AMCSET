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

/*! \file stream.cpp
 *  \brief StreamTableWriter: the file, the column checks and the CSV format
 *  (the Parquet row groups and footer are in parquet.cpp).
 */

#include "output/stream.h"
#include "output/text.h"

#include <charconv>
#include <cmath>
#include <stdexcept>
#include <system_error>

namespace amcset {
namespace output {

namespace {

//! Parquet row groups of at least this many rows (fewer only at the end).
constexpr std::size_t parquet_group_rows = 256 * 1024;

} // namespace

StreamTableWriter::StreamTableWriter(const std::filesystem::path &path, const Table &schema,
                                     std::vector<std::pair<std::string, std::string>> metadata)
    : path_(path), schema_(schema), metadata_(std::move(metadata)) {
  const auto ext = path.extension().string();
  if (ext == ".csv") {
    format_ = Format::csv;
  } else if (ext == ".parquet") {
    format_ = Format::parquet;
  } else {
    throw std::invalid_argument("'" + path.string() + "': the file name must end in .csv or .parquet");
  }
  for (auto &c : schema_.columns) {
    c.values.clear();
    c.text.clear();
  }
  if (path.has_parent_path()) {
    std::error_code ec;
    std::filesystem::create_directories(path.parent_path(), ec);
  }
  out_.open(path, std::ios::binary | std::ios::trunc);
  if (!out_) throw std::runtime_error("Cannot write " + path.string());
  if (format_ == Format::parquet) {
    parquet_open();
    return;
  }
  // CSV: one header line of the column keys.
  std::string header;
  for (std::size_t k = 0; k < schema_.columns.size(); ++k) {
    if (k) header += ',';
    header += schema_.columns[k].key;
  }
  header += '\n';
  out_.write(header.data(), static_cast<std::streamsize>(header.size()));
}

void StreamTableWriter::append(const Table &rows) {
  if (finished_) throw std::logic_error("StreamTableWriter::append after finish");
  if (rows.columns.size() != schema_.columns.size()) {
    throw std::invalid_argument("StreamTableWriter: " + std::to_string(rows.columns.size()) + " columns, expected " +
                                std::to_string(schema_.columns.size()));
  }
  const std::size_t n = rows.rows();
  for (std::size_t k = 0; k < rows.columns.size(); ++k) {
    const auto &c = rows.columns[k];
    if (c.key != schema_.columns[k].key || c.type != schema_.columns[k].type || c.size() != n) {
      throw std::invalid_argument("StreamTableWriter: column " + std::to_string(k) + " ('" + c.key +
                                  "') does not match the file's columns");
    }
  }
  rows_ += n;
  if (format_ == Format::parquet) {
    for (std::size_t k = 0; k < rows.columns.size(); ++k) {
      auto &dst = pending_.columns[k];
      const auto &src = rows.columns[k];
      if (src.type == ColumnType::text) {
        dst.text.insert(dst.text.end(), src.text.begin(), src.text.end());
      } else {
        dst.values.insert(dst.values.end(), src.values.begin(), src.values.end());
      }
    }
    if (pending_.rows() >= parquet_group_rows) flush_parquet();
    return;
  }
  // CSV.
  text::BufferedWriter w(out_);
  NumberFormat shortest;
  shortest.notation = NumberFormat::Notation::shortest;
  shortest.uppercase = false;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < rows.columns.size(); ++k) {
      const auto &c = rows.columns[k];
      if (k) w.put(',');
      switch (c.type) {
      case ColumnType::text: w.put(c.text[i]); break;
      case ColumnType::integer: {
        char *q = w.reserve(text::max_field);
        const double v = c.values[i];
        if (std::isfinite(v)) {
          q = std::to_chars(q, q + text::max_field, static_cast<long long>(std::llround(v))).ptr;
        }
        w.commit(q);
        break;
      }
      case ColumnType::real: {
        char *q = w.reserve(text::max_field);
        w.commit(format_number_to(q, c.values[i], shortest));
        break;
      }
      }
    }
    w.put('\n');
  }
  w.finish();
}

void StreamTableWriter::finish() {
  if (finished_) return;
  finished_ = true;
  if (format_ == Format::parquet) parquet_close();
  out_.flush();
  if (!out_) throw std::runtime_error("Error writing " + path_.string() + " (disk full?)");
  out_.close();
}

} // namespace output
} // namespace amcset
