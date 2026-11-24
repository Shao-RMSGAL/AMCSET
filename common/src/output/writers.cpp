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

/*! \file writers.cpp
 *  \brief CSV and JSON writers and the dataset writer.
 */

#include "output/writers.h"
#include "output/text.h"

#include "parallel_tasks.h"

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <ostream>
#include <stdexcept>

namespace amcset {
namespace output {

namespace {

std::string csv_field(const std::string &s, char delim) {
  bool quote = s.find(delim) != std::string::npos || s.find('"') != std::string::npos ||
               s.find('\n') != std::string::npos || s.find('\r') != std::string::npos;
  if (!quote) return s;
  std::string out = "\"";
  for (char c : s) {
    if (c == '"') out += '"';
    out += c;
  }
  out += '"';
  return out;
}

std::string heading(const Column &c) {
  const std::string label = c.label.empty() ? c.key : c.label;
  if (c.type != ColumnType::real || c.unit.empty() || c.unit == "1") return label;
  return label + " [" + c.unit + "]";
}

//! An integer column's value ("%.0f") into \p p (room for text::max_field);
//! nothing if it is not finite.
char *integer_text_to(char *p, double v) {
  if (!std::isfinite(v)) return p;
  return text::put_fixed(p, v, 0, 0);
}

//! json_number() into \p p (room for text::max_field + 1).
char *json_number_to(char *p, double v, const NumberFormat &f) {
  if (!std::isfinite(v)) {
    std::memcpy(p, "null", 4);
    return p + 4;
  }
  // JSON needs a digit before the decimal point (".5" -> "0.5", "-.5" ->
  // "-0.5"); format_number never writes a leading '+'.
  char *q = p;
  char *e = format_number_to(q, v, f);
  const std::size_t at = (e > q && q[0] == '.') ? 0 : (e - q > 1 && q[0] == '-' && q[1] == '.') ? 1 : SIZE_MAX;
  if (at != SIZE_MAX) {
    std::memmove(q + at + 1, q + at, static_cast<std::size_t>(e - (q + at)));
    q[at] = '0';
    ++e;
  }
  return e;
}

void json_string(std::ostream &o, const std::string &s) {
  o << '"';
  for (unsigned char c : s) {
    switch (c) {
    case '"': o << "\\\""; break;
    case '\\': o << "\\\\"; break;
    case '\n': o << "\\n"; break;
    case '\r': o << "\\r"; break;
    case '\t': o << "\\t"; break;
    default:
      if (c < 0x20) {
        char buf[8];
        std::snprintf(buf, sizeof buf, "\\u%04x", c);
        o << buf;
      } else {
        o << static_cast<char>(c);
      }
    }
  }
  o << '"';
}

void json_number(std::ostream &o, double v, const NumberFormat &f) {
  char buf[text::max_field + 1];
  o.write(buf, json_number_to(buf, v, f) - buf);
}

std::ofstream open(const std::filesystem::path &p) {
  std::ofstream out(p, std::ios::binary);
  if (!out) throw std::runtime_error("Cannot write " + p.string());
  return out;
}

} // namespace

void write_csv(std::ostream &o, const Table &t, const Dataset &context, const OutputSpec &spec) {
  const char d = spec.csv_delimiter;
  if (spec.csv_comments) {
    o << "# " << (t.title.empty() ? t.key : t.title) << "\n";
    if (!t.description.empty()) {
      std::string line;
      for (char c : t.description) {
        if (c == '\n') {
          o << "# " << line << "\n";
          line.clear();
        } else {
          line += c;
        }
      }
      if (!line.empty()) o << "# " << line << "\n";
    }
    if (!context.title.empty()) o << "# run: " << context.title << "\n";
    for (const auto &kv : context.metadata) o << "# " << kv.first << ": " << kv.second << "\n";
    for (const auto &c : t.columns) {
      if (!c.description.empty()) o << "# " << heading(c) << ": " << c.description << "\n";
    }
  }
  for (std::size_t k = 0; k < t.columns.size(); ++k) {
    if (k) o << d;
    o << csv_field(heading(t.columns[k]), d);
  }
  o << "\n";
  std::vector<NumberFormat> formats;
  for (const auto &c : t.columns) formats.push_back(column_format(t, c, spec));
  const std::size_t n = t.rows();
  // Rows are assembled in a buffer and written in large chunks.
  text::BufferedWriter w(o);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < t.columns.size(); ++k) {
      if (k) w.put(d);
      const Column &c = t.columns[k];
      switch (c.type) {
      case ColumnType::text: w.put(csv_field(i < c.text.size() ? c.text[i] : std::string(), d)); break;
      case ColumnType::integer: {
        char *q = w.reserve(text::max_field);
        w.commit(integer_text_to(q, i < c.values.size() ? c.values[i] : 0.0));
        break;
      }
      case ColumnType::real: {
        char *q = w.reserve(text::max_field);
        w.commit(format_number_to(q, i < c.values.size() ? c.values[i] : 0.0, formats[k]));
        break;
      }
      }
    }
    w.put('\n');
  }
  w.finish();
}

void write_json(std::ostream &o, const Dataset &data, const OutputSpec &spec) {
  const bool pretty = spec.json_pretty;
  const char *nl = pretty ? "\n" : "";
  const auto ind = [&](int n) -> std::string { return pretty ? std::string(2 * n, ' ') : std::string(); };
  const char *sep = pretty ? ": " : ":";
  o << "{" << nl;
  o << ind(1) << "\"format\"" << sep << "\"amcset-results\"," << nl;
  o << ind(1) << "\"format_version\"" << sep << "1," << nl;
  o << ind(1) << "\"title\"" << sep;
  json_string(o, data.title);
  o << "," << nl << ind(1) << "\"metadata\"" << sep << "{";
  for (std::size_t i = 0; i < data.metadata.size(); ++i) {
    o << (i ? "," : "") << nl << ind(2);
    json_string(o, data.metadata[i].first);
    o << sep;
    json_string(o, data.metadata[i].second);
  }
  o << nl << ind(1) << "}," << nl;
  if (!data.configuration_json.empty()) {
    o << ind(1) << "\"configuration\"" << sep << data.configuration_json << "," << nl;
  }
  o << ind(1) << "\"scalars\"" << sep << "[";
  const auto scalars = apply(data.scalars, spec);
  for (std::size_t i = 0; i < scalars.size(); ++i) {
    const auto &s = scalars[i];
    o << (i ? "," : "") << nl << ind(2) << "{\"key\"" << sep;
    json_string(o, s.key);
    o << ", \"label\"" << sep;
    json_string(o, s.label);
    o << ", \"value\"" << sep;
    json_number(o, s.value, spec.format);
    o << ", \"unit\"" << sep;
    json_string(o, s.unit);
    o << ", \"dimension\"" << sep;
    json_string(o, to_string(s.dimension));
    if (!s.description.empty()) {
      o << ", \"description\"" << sep;
      json_string(o, s.description);
    }
    o << "}";
  }
  o << nl << ind(1) << "]," << nl;
  o << ind(1) << "\"tables\"" << sep << "[";
  bool first_table = true;
  for (const auto &raw : data.tables) {
    if (!spec.table(raw.key).include) continue;
    const Table t = apply(raw, spec);
    o << (first_table ? "" : ",") << nl << ind(2) << "{" << nl;
    first_table = false;
    o << ind(3) << "\"key\"" << sep;
    json_string(o, t.key);
    o << "," << nl << ind(3) << "\"title\"" << sep;
    json_string(o, t.title);
    o << "," << nl << ind(3) << "\"description\"" << sep;
    json_string(o, t.description);
    o << "," << nl << ind(3) << "\"rows\"" << sep << t.rows() << "," << nl;
    o << ind(3) << "\"columns\"" << sep << "[";
    for (std::size_t k = 0; k < t.columns.size(); ++k) {
      const Column &c = t.columns[k];
      o << (k ? "," : "") << nl << ind(4) << "{\"key\"" << sep;
      json_string(o, c.key);
      o << ", \"label\"" << sep;
      json_string(o, c.label);
      o << ", \"type\"" << sep
        << (c.type == ColumnType::real ? "\"real\"" : c.type == ColumnType::integer ? "\"integer\"" : "\"text\"");
      if (c.type == ColumnType::real) {
        o << ", \"unit\"" << sep;
        json_string(o, c.unit);
        o << ", \"dimension\"" << sep;
        json_string(o, to_string(c.dimension));
      }
      if (!c.description.empty()) {
        o << ", \"description\"" << sep;
        json_string(o, c.description);
      }
      o << "}";
    }
    o << nl << ind(3) << "]," << nl;
    o << ind(3) << "\"data\"" << sep << "{";
    for (std::size_t k = 0; k < t.columns.size(); ++k) {
      const Column &c = t.columns[k];
      const NumberFormat f = column_format(t, c, spec);
      o << (k ? "," : "") << nl << ind(4);
      json_string(o, c.key);
      o << sep << "[";
      const std::size_t n = c.size();
      if (c.type == ColumnType::text) {
        for (std::size_t i = 0; i < n; ++i) {
          if (i) o << ",";
          json_string(o, c.text[i]);
        }
      } else {
        // Numbers: assembled in a buffer, written in large chunks.
        text::BufferedWriter w(o);
        for (std::size_t i = 0; i < n; ++i) {
          char *q = w.reserve(text::max_field + 2);
          if (i) *q++ = ',';
          if (c.type == ColumnType::integer) {
            char *e = integer_text_to(q, c.values[i]);
            if (e == q) {
              std::memcpy(q, "null", 4);
              e = q + 4;
            }
            w.commit(e);
          } else {
            w.commit(json_number_to(q, c.values[i], f));
          }
        }
        w.finish();
      }
      o << "]";
    }
    o << nl << ind(3) << "}" << nl << ind(2) << "}";
  }
  o << nl << ind(1) << "]" << nl << "}" << nl;
}

std::vector<std::filesystem::path> write_dataset(const std::filesystem::path &dir,
                                                 const Dataset &data, const OutputSpec &spec) {
  spec.validate(data);
  std::vector<std::filesystem::path> files;
  const bool any = spec.formats.count(Format::csv) || spec.formats.count(Format::json) ||
                   spec.formats.count(Format::parquet);
  if (!any) return files;
  std::filesystem::create_directories(dir);
  if (spec.formats.count(Format::csv)) std::filesystem::create_directories(dir / "csv");
  if (spec.formats.count(Format::parquet)) std::filesystem::create_directories(dir / "parquet");
  // The files are independent: one task per table (CSV and Parquet) and one
  // for the JSON document, run concurrently. The list keeps the table order.
  std::vector<std::function<void()>> tasks;
  const bool tables = spec.formats.count(Format::csv) || spec.formats.count(Format::parquet);
  for (const auto &raw : data.tables) {
    const TableSpec &ts = spec.table(raw.key);
    if (!ts.include || !tables) continue;
    const std::string base = ts.file_name.empty() ? raw.key : ts.file_name;
    if (spec.formats.count(Format::csv)) files.push_back(dir / "csv" / (base + ".csv"));
    if (spec.formats.count(Format::parquet)) files.push_back(dir / "parquet" / (base + ".parquet"));
    tasks.push_back([&raw, &spec, &data, &dir, base] {
      const Table t = apply(raw, spec);
      if (spec.formats.count(Format::csv)) {
        const auto p = dir / "csv" / (base + ".csv");
        auto o = open(p);
        write_csv(o, t, data, spec);
        if (!o) throw std::runtime_error("Error writing " + p.string());
      }
      if (spec.formats.count(Format::parquet)) {
        write_parquet(dir / "parquet" / (base + ".parquet"), t, data, spec);
      }
    });
  }
  if (spec.formats.count(Format::json)) {
    const auto p = dir / "results.json";
    files.push_back(p);
    tasks.push_back([p, &data, &spec] {
      auto o = open(p);
      write_json(o, data, spec);
      if (!o) throw std::runtime_error("Error writing " + p.string());
    });
  }
  common::run_tasks(tasks);
  return files;
}

} // namespace output
} // namespace amcset
