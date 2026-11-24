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

/*! \file parquet.cpp
 *  \brief A small Apache Parquet writer for output tables.
 *
 * Writes Parquet format version 1 files as specified in the Apache Parquet
 * documentation (parquet.thrift): magic "PAR1", one row group whose column
 * chunks are sequences of uncompressed data pages (v1) with PLAIN-encoded
 * values of REQUIRED columns (so there are no repetition or definition
 * levels), then the FileMetaData structure serialised with the Thrift
 * compact protocol, its length and the closing magic. That subset is enough
 * for tables of doubles, 64-bit integers and strings and is read by every
 * Parquet implementation (Arrow/pandas/polars/DuckDB/Spark).
 *
 * StreamTableWriter's Parquet files (output/stream.h) are the same format
 * with several row groups, written one after another as rows arrive; the
 * footer listing them comes last.
 */

#include "output/stream.h"
#include "output/writers.h"

#include "config.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>

namespace amcset {
namespace output {

namespace {

// --- Thrift compact protocol (only what FileMetaData and PageHeader need) ---

enum CType : std::uint8_t {
  kStop = 0,
  kTrue = 1,
  kFalse = 2,
  kI32 = 5,
  kI64 = 6,
  kBinary = 8,
  kList = 9,
  kStruct = 12,
};

class Compact {
public:
  std::vector<std::uint8_t> bytes;

  void begin_struct() { last_.push_back(0); }
  void end_struct() {
    bytes.push_back(kStop);
    last_.pop_back();
  }
  void field_i32(int id, std::int32_t v) {
    header(id, kI32);
    varint(zigzag(v));
  }
  void field_i64(int id, std::int64_t v) {
    header(id, kI64);
    varint(zigzag(v));
  }
  void field_bool(int id, bool v) { header(id, v ? kTrue : kFalse); }
  void field_string(int id, const std::string &s) {
    header(id, kBinary);
    string(s);
  }
  void field_struct(int id) {
    header(id, kStruct);
    begin_struct();
  }
  void field_list(int id, CType element, std::size_t n) {
    header(id, kList);
    list_header(element, n);
  }
  // List elements.
  void element_i32(std::int32_t v) { varint(zigzag(v)); }
  void element_string(const std::string &s) { string(s); }

private:
  std::vector<int> last_;

  static std::uint64_t zigzag(std::int64_t v) {
    return (static_cast<std::uint64_t>(v) << 1) ^ static_cast<std::uint64_t>(v >> 63);
  }
  void varint(std::uint64_t v) {
    while (v >= 0x80) {
      bytes.push_back(static_cast<std::uint8_t>(v | 0x80));
      v >>= 7;
    }
    bytes.push_back(static_cast<std::uint8_t>(v));
  }
  void header(int id, CType t) {
    int &last = last_.back();
    const int delta = id - last;
    if (delta > 0 && delta <= 15) {
      bytes.push_back(static_cast<std::uint8_t>((delta << 4) | t));
    } else {
      bytes.push_back(t);
      varint(zigzag(static_cast<std::int16_t>(id)));
    }
    last = id;
  }
  void list_header(CType element, std::size_t n) {
    if (n < 15) {
      bytes.push_back(static_cast<std::uint8_t>((n << 4) | element));
    } else {
      bytes.push_back(static_cast<std::uint8_t>(0xF0 | element));
      varint(n);
    }
  }
  void string(const std::string &s) {
    varint(s.size());
    bytes.insert(bytes.end(), s.begin(), s.end());
  }
};

// Parquet enums (parquet.thrift).
constexpr std::int32_t kTypeInt64 = 2;
constexpr std::int32_t kTypeDouble = 5;
constexpr std::int32_t kTypeByteArray = 6;
constexpr std::int32_t kRequired = 0;
constexpr std::int32_t kConvertedUtf8 = 0;
constexpr std::int32_t kEncodingPlain = 0;
constexpr std::int32_t kEncodingRle = 3;
constexpr std::int32_t kCodecUncompressed = 0;
constexpr std::int32_t kPageData = 0;

void put_le64(std::vector<std::uint8_t> &b, std::uint64_t v) {
  for (int i = 0; i < 8; ++i) b.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}
void put_le32(std::vector<std::uint8_t> &b, std::uint32_t v) {
  for (int i = 0; i < 4; ++i) b.push_back(static_cast<std::uint8_t>(v >> (8 * i)));
}

std::int32_t physical_type(const Column &c) {
  switch (c.type) {
  case ColumnType::real: return kTypeDouble;
  case ColumnType::integer: return kTypeInt64;
  case ColumnType::text: return kTypeByteArray;
  }
  return kTypeDouble;
}

// PLAIN encoding of rows [first, last).
std::vector<std::uint8_t> plain_values(const Column &c, std::size_t first, std::size_t last) {
  std::vector<std::uint8_t> b;
  switch (c.type) {
  case ColumnType::real:
    b.reserve(8 * (last - first));
    for (std::size_t i = first; i < last; ++i) {
      std::uint64_t u;
      static_assert(sizeof u == sizeof(double));
      std::memcpy(&u, &c.values[i], sizeof u);
      put_le64(b, u);
    }
    break;
  case ColumnType::integer:
    b.reserve(8 * (last - first));
    for (std::size_t i = first; i < last; ++i) {
      const double v = c.values[i];
      const auto x = std::isfinite(v) ? static_cast<std::int64_t>(std::llround(v)) : std::int64_t{0};
      put_le64(b, static_cast<std::uint64_t>(x));
    }
    break;
  case ColumnType::text:
    for (std::size_t i = first; i < last; ++i) {
      const std::string &s = c.text[i];
      put_le32(b, static_cast<std::uint32_t>(s.size()));
      b.insert(b.end(), s.begin(), s.end());
    }
    break;
  }
  return b;
}

std::vector<std::uint8_t> page_header(std::size_t values, std::size_t data_bytes) {
  Compact t;
  t.begin_struct();
  t.field_i32(1, kPageData);
  t.field_i32(2, static_cast<std::int32_t>(data_bytes));
  t.field_i32(3, static_cast<std::int32_t>(data_bytes));
  t.field_struct(5); // DataPageHeader
  t.field_i32(1, static_cast<std::int32_t>(values));
  t.field_i32(2, kEncodingPlain);
  t.field_i32(3, kEncodingRle);
  t.field_i32(4, kEncodingRle);
  t.end_struct();
  t.end_struct();
  return t.bytes;
}

std::string json_escape(const std::string &s) {
  std::string o;
  for (unsigned char c : s) {
    if (c == '"' || c == '\\') {
      o += '\\';
      o += static_cast<char>(c);
    } else if (c < 0x20) {
      char buf[8];
      std::snprintf(buf, sizeof buf, "\\u%04x", c);
      o += buf;
    } else {
      o += static_cast<char>(c);
    }
  }
  return o;
}

} // namespace

void write_parquet(const std::filesystem::path &path, const Table &table, const Dataset &context,
                   const OutputSpec &spec) {
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("Cannot write " + path.string());
  const std::size_t rows = table.rows();
  for (const auto &c : table.columns) {
    if (c.size() != rows) {
      throw std::runtime_error("Parquet: column '" + c.key + "' of table '" + table.key +
                               "' has " + std::to_string(c.size()) + " values, expected " +
                               std::to_string(rows));
    }
  }
  std::vector<std::uint8_t> file = {'P', 'A', 'R', '1'};

  struct ChunkInfo {
    std::int64_t data_page_offset = 0;
    std::int64_t size = 0;
  };
  std::vector<ChunkInfo> chunks;
  constexpr std::size_t page_rows = 1u << 17;
  std::int64_t total_bytes = 0;
  for (const auto &c : table.columns) {
    ChunkInfo ci;
    ci.data_page_offset = static_cast<std::int64_t>(file.size());
    std::size_t first = 0;
    do {
      const std::size_t last = std::min(rows, first + page_rows);
      const auto data = plain_values(c, first, last);
      if (data.size() > 0x7fffffffu) throw std::runtime_error("Parquet: page too large");
      const auto header = page_header(last - first, data.size());
      file.insert(file.end(), header.begin(), header.end());
      file.insert(file.end(), data.begin(), data.end());
      first = last;
    } while (first < rows);
    ci.size = static_cast<std::int64_t>(file.size()) - ci.data_page_offset;
    total_bytes += ci.size;
    chunks.push_back(ci);
  }

  // Column names.
  std::vector<std::string> names;
  for (const auto &c : table.columns) {
    std::string n = c.key;
    if (spec.parquet_units_in_names && c.type == ColumnType::real && !c.unit.empty() && c.unit != "1") {
      n += " [" + c.unit + "]";
    }
    names.push_back(n);
  }

  // Key-value metadata.
  std::vector<std::pair<std::string, std::string>> kv;
  {
    std::string units = "{", labels = "{", dims = "{";
    for (std::size_t k = 0; k < table.columns.size(); ++k) {
      const auto &c = table.columns[k];
      const std::string sep = k ? "," : "";
      units += sep + "\"" + json_escape(names[k]) + "\":\"" + json_escape(c.unit) + "\"";
      labels += sep + "\"" + json_escape(names[k]) + "\":\"" + json_escape(c.label) + "\"";
      dims += sep + "\"" + json_escape(names[k]) + "\":\"" +
              (c.type == ColumnType::real ? to_string(c.dimension) : std::string()) + "\"";
    }
    kv.emplace_back("amcset.table", table.key);
    kv.emplace_back("amcset.title", table.title);
    kv.emplace_back("amcset.description", table.description);
    kv.emplace_back("amcset.units", units + "}");
    kv.emplace_back("amcset.labels", labels + "}");
    kv.emplace_back("amcset.dimensions", dims + "}");
    std::string meta = "{";
    for (std::size_t i = 0; i < context.metadata.size(); ++i) {
      meta += (i ? "," : "") + std::string("\"") + json_escape(context.metadata[i].first) + "\":\"" +
              json_escape(context.metadata[i].second) + "\"";
    }
    kv.emplace_back("amcset.metadata", meta + "}");
    kv.emplace_back("amcset.run_title", context.title);
  }

  Compact t;
  t.begin_struct(); // FileMetaData
  t.field_i32(1, 1);
  t.field_list(2, kStruct, table.columns.size() + 1);
  {
    // Root.
    t.begin_struct();
    t.field_string(4, "schema");
    t.field_i32(5, static_cast<std::int32_t>(table.columns.size()));
    t.end_struct();
    for (std::size_t k = 0; k < table.columns.size(); ++k) {
      const auto &c = table.columns[k];
      t.begin_struct();
      t.field_i32(1, physical_type(c));
      t.field_i32(3, kRequired);
      t.field_string(4, names[k]);
      if (c.type == ColumnType::text) t.field_i32(6, kConvertedUtf8);
      t.end_struct();
    }
  }
  t.field_i64(3, static_cast<std::int64_t>(rows));
  t.field_list(4, kStruct, 1);
  {
    t.begin_struct(); // RowGroup
    t.field_list(1, kStruct, table.columns.size());
    for (std::size_t k = 0; k < table.columns.size(); ++k) {
      const auto &c = table.columns[k];
      t.begin_struct(); // ColumnChunk
      t.field_i64(2, chunks[k].data_page_offset);
      t.field_struct(3); // ColumnMetaData
      t.field_i32(1, physical_type(c));
      t.field_list(2, static_cast<CType>(kI32), 2);
      t.element_i32(kEncodingPlain);
      t.element_i32(kEncodingRle);
      t.field_list(3, kBinary, 1);
      t.element_string(names[k]);
      t.field_i32(4, kCodecUncompressed);
      t.field_i64(5, static_cast<std::int64_t>(rows));
      t.field_i64(6, chunks[k].size);
      t.field_i64(7, chunks[k].size);
      t.field_i64(9, chunks[k].data_page_offset);
      t.end_struct(); // ColumnMetaData
      t.end_struct(); // ColumnChunk
    }
    t.field_i64(2, total_bytes);
    t.field_i64(3, static_cast<std::int64_t>(rows));
    t.end_struct();
  }
  t.field_list(5, kStruct, kv.size());
  for (const auto &p : kv) {
    t.begin_struct();
    t.field_string(1, p.first);
    t.field_string(2, p.second);
    t.end_struct();
  }
  t.field_string(6, std::string("AMCSET version ") + AMCSET_VERSION);
  t.end_struct();

  file.insert(file.end(), t.bytes.begin(), t.bytes.end());
  put_le32(file, static_cast<std::uint32_t>(t.bytes.size()));
  file.insert(file.end(), {'P', 'A', 'R', '1'});
  out.write(reinterpret_cast<const char *>(file.data()), static_cast<std::streamsize>(file.size()));
  if (!out) throw std::runtime_error("Error writing " + path.string());
}

// ---------------------------------------------------------------------------
// StreamTableWriter: Parquet with one row group per flush
// ---------------------------------------------------------------------------

void StreamTableWriter::parquet_open() {
  parquet_ = std::make_unique<ParquetState>();
  for (const auto &c : schema_.columns) parquet_->names.push_back(c.key);
  out_.write("PAR1", 4);
  parquet_->offset = 4;
  pending_ = schema_;
}

void StreamTableWriter::flush_parquet() {
  const std::size_t rows = pending_.rows();
  if (rows == 0) return;
  ParquetState::RowGroup g;
  g.rows = static_cast<std::int64_t>(rows);
  constexpr std::size_t page_rows = 1u << 17;
  for (const auto &c : pending_.columns) {
    ParquetState::Chunk ch;
    ch.offset = parquet_->offset;
    std::size_t first = 0;
    do {
      const std::size_t last = std::min(rows, first + page_rows);
      const auto data = plain_values(c, first, last);
      if (data.size() > 0x7fffffffu) throw std::runtime_error("Parquet: page too large");
      const auto header = page_header(last - first, data.size());
      out_.write(reinterpret_cast<const char *>(header.data()), static_cast<std::streamsize>(header.size()));
      out_.write(reinterpret_cast<const char *>(data.data()), static_cast<std::streamsize>(data.size()));
      parquet_->offset += static_cast<std::int64_t>(header.size() + data.size());
      first = last;
    } while (first < rows);
    ch.size = parquet_->offset - ch.offset;
    g.bytes += ch.size;
    g.chunks.push_back(ch);
  }
  if (!out_) throw std::runtime_error("Error writing " + path_.string());
  parquet_->groups.push_back(std::move(g));
  for (auto &c : pending_.columns) {
    c.values.clear();
    c.text.clear();
  }
}

void StreamTableWriter::parquet_close() {
  flush_parquet();
  std::int64_t total_rows = 0;
  for (const auto &g : parquet_->groups) total_rows += g.rows;
  const auto &cols = schema_.columns;
  std::vector<std::pair<std::string, std::string>> kv;
  {
    std::string units = "{", labels = "{", descriptions = "{";
    for (std::size_t k = 0; k < cols.size(); ++k) {
      const std::string sep = k ? "," : "";
      const std::string n = json_escape(cols[k].key);
      units += sep + "\"" + n + "\":\"" + json_escape(cols[k].unit) + "\"";
      labels += sep + "\"" + n + "\":\"" + json_escape(cols[k].label) + "\"";
      descriptions += sep + "\"" + n + "\":\"" + json_escape(cols[k].description) + "\"";
    }
    kv.emplace_back("amcset.table", schema_.key);
    kv.emplace_back("amcset.title", schema_.title);
    kv.emplace_back("amcset.description", schema_.description);
    kv.emplace_back("amcset.units", units + "}");
    kv.emplace_back("amcset.labels", labels + "}");
    kv.emplace_back("amcset.descriptions", descriptions + "}");
    for (const auto &m : metadata_) kv.emplace_back("amcset." + m.first, m.second);
  }
  Compact t;
  t.begin_struct(); // FileMetaData
  t.field_i32(1, 1);
  t.field_list(2, kStruct, cols.size() + 1);
  t.begin_struct();
  t.field_string(4, "schema");
  t.field_i32(5, static_cast<std::int32_t>(cols.size()));
  t.end_struct();
  for (std::size_t k = 0; k < cols.size(); ++k) {
    t.begin_struct();
    t.field_i32(1, physical_type(cols[k]));
    t.field_i32(3, kRequired);
    t.field_string(4, parquet_->names[k]);
    if (cols[k].type == ColumnType::text) t.field_i32(6, kConvertedUtf8);
    t.end_struct();
  }
  t.field_i64(3, total_rows);
  t.field_list(4, kStruct, parquet_->groups.size());
  for (const auto &g : parquet_->groups) {
    t.begin_struct(); // RowGroup
    t.field_list(1, kStruct, cols.size());
    for (std::size_t k = 0; k < cols.size(); ++k) {
      t.begin_struct(); // ColumnChunk
      t.field_i64(2, g.chunks[k].offset);
      t.field_struct(3); // ColumnMetaData
      t.field_i32(1, physical_type(cols[k]));
      t.field_list(2, static_cast<CType>(kI32), 2);
      t.element_i32(kEncodingPlain);
      t.element_i32(kEncodingRle);
      t.field_list(3, kBinary, 1);
      t.element_string(parquet_->names[k]);
      t.field_i32(4, kCodecUncompressed);
      t.field_i64(5, g.rows);
      t.field_i64(6, g.chunks[k].size);
      t.field_i64(7, g.chunks[k].size);
      t.field_i64(9, g.chunks[k].offset);
      t.end_struct(); // ColumnMetaData
      t.end_struct(); // ColumnChunk
    }
    t.field_i64(2, g.bytes);
    t.field_i64(3, g.rows);
    t.end_struct();
  }
  t.field_list(5, kStruct, kv.size());
  for (const auto &p : kv) {
    t.begin_struct();
    t.field_string(1, p.first);
    t.field_string(2, p.second);
    t.end_struct();
  }
  t.field_string(6, std::string("AMCSET version ") + AMCSET_VERSION);
  t.end_struct();
  out_.write(reinterpret_cast<const char *>(t.bytes.data()), static_cast<std::streamsize>(t.bytes.size()));
  std::vector<std::uint8_t> len;
  put_le32(len, static_cast<std::uint32_t>(t.bytes.size()));
  out_.write(reinterpret_cast<const char *>(len.data()), 4);
  out_.write("PAR1", 4);
}

} // namespace output
} // namespace amcset
