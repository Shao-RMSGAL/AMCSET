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

/*! \file table.cpp
 *  \brief Output tables, number formats and output specifications.
 */

#include "output/table.h"
#include "output/text.h"

#include <boost/json.hpp>

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>

namespace amcset {
namespace output {

namespace json = boost::json;

// ---------------------------------------------------------------------------
// Columns, tables, datasets
// ---------------------------------------------------------------------------

Column Column::real(std::string key, std::string label, Dimension d, std::string unit,
                    std::vector<double> values, std::string description) {
  Column c;
  c.key = std::move(key);
  c.label = std::move(label);
  c.type = ColumnType::real;
  c.dimension = d;
  c.unit = std::move(unit);
  (void)find_unit(c.unit, d); // the unit must exist and fit the dimension
  c.values = std::move(values);
  c.description = std::move(description);
  return c;
}

Column Column::integer(std::string key, std::string label, std::vector<double> values,
                       std::string description) {
  Column c;
  c.key = std::move(key);
  c.label = std::move(label);
  c.type = ColumnType::integer;
  c.values = std::move(values);
  c.description = std::move(description);
  return c;
}

Column Column::strings(std::string key, std::string label, std::vector<std::string> text,
                       std::string description) {
  Column c;
  c.key = std::move(key);
  c.label = std::move(label);
  c.type = ColumnType::text;
  c.unit.clear();
  c.text = std::move(text);
  c.description = std::move(description);
  return c;
}

const Column *Table::find(std::string_view k) const {
  for (const auto &c : columns) {
    if (c.key == k) return &c;
  }
  return nullptr;
}

Column *Table::find(std::string_view k) {
  for (auto &c : columns) {
    if (c.key == k) return &c;
  }
  return nullptr;
}

const Table *Dataset::find(std::string_view k) const {
  for (const auto &t : tables) {
    if (t.key == k) return &t;
  }
  return nullptr;
}

// ---------------------------------------------------------------------------
// Numbers
// ---------------------------------------------------------------------------

std::string to_string(NumberFormat::Notation n) {
  switch (n) {
  case NumberFormat::Notation::shortest: return "shortest";
  case NumberFormat::Notation::general: return "general";
  case NumberFormat::Notation::fixed: return "fixed";
  case NumberFormat::Notation::scientific: return "scientific";
  case NumberFormat::Notation::engineering: return "engineering";
  }
  return "shortest";
}

NumberFormat::Notation parse_notation(std::string_view s) {
  using N = NumberFormat::Notation;
  if (s == "shortest" || s == "default" || s == "auto") return N::shortest;
  if (s == "general" || s == "g") return N::general;
  if (s == "fixed" || s == "f" || s == "decimal") return N::fixed;
  if (s == "scientific" || s == "e" || s == "sci") return N::scientific;
  if (s == "engineering" || s == "eng") return N::engineering;
  throw std::invalid_argument("unknown number notation '" + std::string(s) +
                              "' (shortest, general, fixed, scientific, engineering)");
}

std::string format_number(double v, const NumberFormat &f) {
  char buf[text::max_field];
  return std::string(buf, format_number_to(buf, v, f));
}

char *format_number_to(char *p, double v, const NumberFormat &f) {
  const auto copy = [&p](const char *s, std::size_t n) {
    std::memcpy(p, s, n);
    return p + n;
  };
  if (std::isnan(v)) return copy("nan", 3);
  if (std::isinf(v)) return v > 0 ? copy("inf", 3) : copy("-inf", 4);
  const int d = std::clamp(f.digits, 0, 30);
  char *const end = p + text::max_field;
  // std::to_chars with a precision prints what printf prints with it
  // (%.*g, %.*f, %.*e), without printf's cost (OutputText tests).
  const auto upper = [&f](char *b, char *e) {
    if (f.uppercase) {
      for (char *q = b; q != e; ++q) {
        if (*q == 'e') *q = 'E';
      }
    }
    return e;
  };
  using N = NumberFormat::Notation;
  switch (f.notation) {
  case N::shortest: return upper(p, std::to_chars(p, end, v).ptr);
  case N::general: return upper(p, std::to_chars(p, end, v, std::chars_format::general, std::max(d, 1)).ptr);
  case N::fixed: return std::to_chars(p, end, v, std::chars_format::fixed, d).ptr;
  case N::scientific: return upper(p, std::to_chars(p, end, v, std::chars_format::scientific, d).ptr);
  case N::engineering: break;
  }
  char buf[96];
  {
    if (v == 0.0) {
      std::snprintf(buf, sizeof buf, f.uppercase ? "%.*fE+00" : "%.*fe+00", d, 0.0);
      return copy(buf, std::strlen(buf));
    }
    int e = static_cast<int>(std::floor(std::log10(std::abs(v))));
    e -= ((e % 3) + 3) % 3;
    double m = v / std::pow(10.0, e);
    // Rounding can push the mantissa to 1000.
    char probe[64];
    std::snprintf(probe, sizeof probe, "%.*f", d, std::abs(m));
    if (std::atof(probe) >= 1000.0) {
      e += 3;
      m /= 1000.0;
    }
    std::snprintf(buf, sizeof buf, f.uppercase ? "%.*fE%+03d" : "%.*fe%+03d", d, m, e);
    return copy(buf, std::strlen(buf));
  }
}

// ---------------------------------------------------------------------------
// Formats
// ---------------------------------------------------------------------------

std::string to_string(Format f) {
  switch (f) {
  case Format::trim: return "trim";
  case Format::csv: return "csv";
  case Format::json: return "json";
  case Format::parquet: return "parquet";
  }
  return "csv";
}

Format parse_format(std::string_view s) {
  if (s == "trim" || s == "txt" || s == "srim") return Format::trim;
  if (s == "csv") return Format::csv;
  if (s == "json") return Format::json;
  if (s == "parquet" || s == "pq") return Format::parquet;
  throw std::invalid_argument("unknown output format '" + std::string(s) +
                              "' (trim, csv, json, parquet)");
}

// ---------------------------------------------------------------------------
// Applying a spec
// ---------------------------------------------------------------------------

const TableSpec &OutputSpec::table(std::string_view key) const {
  static const TableSpec defaults;
  auto it = tables.find(std::string(key));
  if (it != tables.end()) return it->second;
  it = tables.find("*");
  return it != tables.end() ? it->second : defaults;
}

namespace {

// Output unit of a column under a spec.
const Unit *output_unit(const Column &c, const OutputSpec &spec, const TableSpec &ts) {
  if (c.type != ColumnType::real) return nullptr;
  const auto cs = ts.column_specs.find(c.key);
  if (cs != ts.column_specs.end() && !cs->second.unit.empty()) {
    return &find_unit(cs->second.unit, c.dimension);
  }
  const auto du = spec.units.find(c.dimension);
  if (du != spec.units.end() && !du->second.empty()) return &find_unit(du->second, c.dimension);
  return &find_unit(c.unit, c.dimension);
}

std::vector<double> converted(const Column &c, const Unit *to) {
  if (to == nullptr) return c.values;
  const double f = conversion_factor(find_unit(c.unit, c.dimension), *to);
  std::vector<double> v(c.values);
  if (f != 1.0) {
    for (double &x : v) x *= f;
  }
  return v;
}

void check_column(const Table &t, const std::string &key, const char *what) {
  if (t.find(key) == nullptr) {
    std::string known;
    for (const auto &c : t.columns) known += (known.empty() ? "" : ", ") + c.key;
    throw std::invalid_argument("table '" + t.key + "': " + what + " '" + key +
                                "' is not a column (columns: " + known + ")");
  }
}

} // namespace

Table apply(const Table &table, const OutputSpec &spec) {
  const TableSpec &ts = spec.table(table.key);
  // Columns, in order.
  std::vector<const Column *> cols;
  if (!ts.columns.empty()) {
    for (const auto &k : ts.columns) {
      check_column(table, k, "column");
      cols.push_back(table.find(k));
    }
  } else {
    for (const auto &c : table.columns) cols.push_back(&c);
  }
  for (const auto &k : ts.exclude) check_column(table, k, "excluded column");
  cols.erase(std::remove_if(cols.begin(), cols.end(),
                            [&](const Column *c) {
                              return std::find(ts.exclude.begin(), ts.exclude.end(), c->key) !=
                                     ts.exclude.end();
                            }),
             cols.end());
  for (const auto &kv : ts.column_specs) check_column(table, kv.first, "column setting for");

  // Rows.
  const std::size_t n = table.rows();
  std::vector<char> keep(n, 1);
  for (const auto &f : ts.filters) {
    check_column(table, f.column, "filter column");
    const Column &c = *table.find(f.column);
    if (c.type == ColumnType::text) {
      throw std::invalid_argument("table '" + table.key + "': cannot filter the text column '" +
                                  c.key + "' by value");
    }
    const Unit *u = output_unit(c, spec, ts);
    const auto v = converted(c, u);
    const auto bound = [&](const std::optional<double> &b, const std::string &unit) {
      if (!b || unit.empty() || u == nullptr) return b;
      return std::optional<double>(*b * conversion_factor(find_unit(unit, c.dimension), *u));
    };
    const auto lo = bound(f.min, f.min_unit);
    const auto hi = bound(f.max, f.max_unit);
    for (std::size_t i = 0; i < n && i < v.size(); ++i) {
      if ((lo && v[i] < *lo) || (hi && v[i] > *hi)) keep[i] = 0;
    }
  }
  if (ts.drop_zero_rows) {
    for (std::size_t i = 0; i < n; ++i) {
      if (!keep[i]) continue;
      bool any = false;
      for (std::size_t k = 1; k < cols.size() && !any; ++k) {
        const Column &c = *cols[k];
        if (c.type != ColumnType::text && i < c.values.size() && c.values[i] != 0.0) any = true;
      }
      if (!any && cols.size() > 1) keep[i] = 0;
    }
  }
  std::vector<std::size_t> rows;
  {
    const std::size_t every = std::max<std::size_t>(ts.every, 1);
    std::size_t seen = 0;
    for (std::size_t i = 0; i < n; ++i) {
      if (!keep[i]) continue;
      if (seen++ % every != 0) continue;
      rows.push_back(i);
      if (ts.max_rows != 0 && rows.size() >= ts.max_rows) break;
    }
  }

  Table out;
  out.key = table.key;
  out.title = table.title;
  out.description = table.description;
  for (const Column *c : cols) {
    Column o;
    o.key = c->key;
    o.label = c->label;
    o.type = c->type;
    o.dimension = c->dimension;
    o.description = c->description;
    const auto cs = ts.column_specs.find(c->key);
    if (cs != ts.column_specs.end() && !cs->second.label.empty()) o.label = cs->second.label;
    if (c->type == ColumnType::text) {
      o.unit.clear();
      o.text.reserve(rows.size());
      for (auto i : rows) o.text.push_back(i < c->text.size() ? c->text[i] : std::string());
    } else {
      const Unit *u = output_unit(*c, spec, ts);
      const auto v = converted(*c, u);
      o.unit = u != nullptr ? u->name : c->unit;
      o.values.reserve(rows.size());
      for (auto i : rows) o.values.push_back(i < v.size() ? v[i] : 0.0);
    }
    out.columns.push_back(std::move(o));
  }
  return out;
}

NumberFormat column_format(const Table &table, const Column &column, const OutputSpec &spec) {
  const TableSpec &ts = spec.table(table.key);
  const auto cs = ts.column_specs.find(column.key);
  if (cs != ts.column_specs.end() && cs->second.format) return *cs->second.format;
  if (ts.format) return *ts.format;
  return spec.format;
}

std::vector<Scalar> apply(const std::vector<Scalar> &scalars, const OutputSpec &spec) {
  std::vector<Scalar> out = scalars;
  for (auto &s : out) {
    const auto du = spec.units.find(s.dimension);
    if (du == spec.units.end() || du->second.empty()) continue;
    const Unit &to = find_unit(du->second, s.dimension);
    s.value *= conversion_factor(find_unit(s.unit, s.dimension), to);
    s.unit = to.name;
  }
  return out;
}

void OutputSpec::validate(const Dataset &data) const {
  for (const auto &kv : units) (void)find_unit(kv.second, kv.first);
  for (const auto &kv : tables) {
    if (kv.first == "*") continue;
    const Table *t = data.find(kv.first);
    if (t == nullptr) {
      std::string known;
      for (const auto &x : data.tables) known += (known.empty() ? "" : ", ") + x.key;
      throw std::invalid_argument("output spec: unknown table '" + kv.first + "' (tables: " +
                                  known + ")");
    }
  }
  for (const auto &t : data.tables) {
    if (!table(t.key).include) continue;
    (void)apply(t, *this); // checks columns, units and filters
  }
}

// ---------------------------------------------------------------------------
// JSON form
// ---------------------------------------------------------------------------

namespace {

std::string str(const json::value &v, const std::string &where) {
  if (!v.is_string()) throw std::invalid_argument(where + ": expected a string");
  return std::string(v.as_string());
}

double num(const json::value &v, const std::string &where) {
  double d = 0.0;
  if (v.is_double()) d = v.as_double();
  else if (v.is_int64()) d = static_cast<double>(v.as_int64());
  else if (v.is_uint64()) d = static_cast<double>(v.as_uint64());
  else throw std::invalid_argument(where + ": expected a number");
  if (!std::isfinite(d)) throw std::invalid_argument(where + ": not a finite number");
  return d;
}

//! A whole number in [lo, hi], checked before the integer conversion.
std::size_t whole(const json::value &v, const std::string &where, double lo, double hi) {
  const double d = num(v, where);
  if (!(d >= lo && d <= hi) || d != std::floor(d)) {
    throw std::invalid_argument(where + ": expected a whole number from " + std::to_string(static_cast<long long>(lo)) +
                                " to " + std::to_string(static_cast<long long>(hi)));
  }
  return static_cast<std::size_t>(d);
}

/*!
 * A table's file name: one plain name, no folders. It is joined to the
 * output folder, so "..", separators, drive letters and characters that are
 * not allowed in Windows file names are rejected.
 */
std::string file_name(const json::value &v, const std::string &where) {
  const std::string s = str(v, where);
  bool ok = !s.empty() && s.size() <= 200 && s != "." && s != ".." && s.back() != '.' && s.back() != ' ';
  for (unsigned char c : s) {
    if (c < 32 || std::string_view("<>:\"/\\|?*").find(static_cast<char>(c)) != std::string_view::npos) ok = false;
  }
  if (!ok) {
    throw std::invalid_argument(where + ": '" + s +
                                "' is not a plain file name (no folders, no < > : \" / \\ | ? *)");
  }
  return s;
}

bool flag(const json::value &v, const std::string &where) {
  if (!v.is_bool()) throw std::invalid_argument(where + ": expected true or false");
  return v.as_bool();
}

// Returns a pointer (not a reference) so that callers passing a temporary
// 'where' string do not look like dangling references to the compiler.
const json::object *objp(const json::value &v, const std::string &where) {
  if (!v.is_object()) throw std::invalid_argument(where + ": expected an object");
  return &v.as_object();
}

std::vector<std::string> strings(const json::value &v, const std::string &where) {
  std::vector<std::string> out;
  if (v.is_string()) {
    out.push_back(std::string(v.as_string()));
    return out;
  }
  if (!v.is_array()) throw std::invalid_argument(where + ": expected a list of strings");
  for (std::size_t i = 0; i < v.as_array().size(); ++i) {
    out.push_back(str(v.as_array()[i], where + "[" + std::to_string(i) + "]"));
  }
  return out;
}

void reject_unknown(const json::object &o, std::initializer_list<const char *> known,
                    const std::string &where) {
  for (const auto &kv : o) {
    bool ok = false;
    for (const char *k : known) ok = ok || kv.key() == k;
    if (!ok) {
      std::string list;
      for (const char *k : known) list += (list.empty() ? "" : ", ") + std::string(k);
      throw std::invalid_argument(where + ": unknown setting '" + std::string(kv.key()) +
                                  "' (expected " + list + ")");
    }
  }
}

NumberFormat parse_format_object(const json::value &v, const std::string &where) {
  NumberFormat f;
  if (v.is_string()) {
    f.notation = parse_notation(std::string(v.as_string()));
    return f;
  }
  const auto &o = *objp(v, where);
  reject_unknown(o, {"notation", "digits", "uppercase"}, where);
  if (const auto *x = o.if_contains("notation")) f.notation = parse_notation(str(*x, where + ".notation"));
  if (const auto *x = o.if_contains("digits")) f.digits = static_cast<int>(whole(*x, where + ".digits", 0, 30));
  if (const auto *x = o.if_contains("uppercase")) f.uppercase = flag(*x, where + ".uppercase");
  return f;
}

// "100 nm" -> (100, "nm"); 100 -> (100, "").
std::pair<double, std::string> bound_value(const json::value &v, const std::string &where) {
  if (!v.is_string()) return {num(v, where), {}};
  const std::string s(v.as_string());
  std::size_t pos = 0;
  double x = 0.0;
  try {
    x = std::stod(s, &pos);
  } catch (const std::exception &) {
    throw std::invalid_argument(where + ": expected a number or \"number unit\"");
  }
  if (!std::isfinite(x)) throw std::invalid_argument(where + ": not a finite number");
  std::string unit = s.substr(pos);
  const auto b = unit.find_first_not_of(" \t");
  unit = b == std::string::npos ? std::string() : unit.substr(b);
  if (!unit.empty()) (void)find_unit(unit);
  return {x, unit};
}

json::object format_object(const NumberFormat &f) {
  json::object o;
  o["notation"] = to_string(f.notation);
  o["digits"] = f.digits;
  o["uppercase"] = f.uppercase;
  return o;
}

} // namespace

OutputSpec parse_output_spec(std::string_view text) {
  json::value root;
  try {
    json::parse_options opt;
    opt.allow_comments = true;
    opt.allow_trailing_commas = true;
    root = json::parse(text, {}, opt);
  } catch (const std::exception &e) {
    throw std::invalid_argument(std::string("output spec: invalid JSON: ") + e.what());
  }
  OutputSpec spec;
  const auto &o = *objp(root, "output spec");
  reject_unknown(o, {"formats", "units", "number_format", "csv", "json", "parquet", "tables",
                     "$schema", "comment"},
                 "output spec");
  if (const auto *v = o.if_contains("formats")) {
    spec.formats.clear();
    for (const auto &f : strings(*v, "formats")) spec.formats.insert(parse_format(f));
  }
  if (const auto *v = o.if_contains("units")) {
    for (const auto &kv : *objp(*v, "units")) {
      const Dimension d = parse_dimension(std::string(kv.key()));
      const std::string u = str(kv.value(), "units." + std::string(kv.key()));
      (void)find_unit(u, d);
      spec.units[d] = u;
    }
  }
  if (const auto *v = o.if_contains("number_format")) spec.format = parse_format_object(*v, "number_format");
  if (const auto *v = o.if_contains("csv")) {
    const auto &c = *objp(*v, "csv");
    reject_unknown(c, {"delimiter", "comments"}, "csv");
    if (const auto *x = c.if_contains("delimiter")) {
      const auto d = str(*x, "csv.delimiter");
      if (d == "tab" || d == "\\t" || d == "\t") spec.csv_delimiter = '\t';
      else if (d.size() == 1) spec.csv_delimiter = d[0];
      else throw std::invalid_argument("csv.delimiter: one character or \"tab\"");
    }
    if (const auto *x = c.if_contains("comments")) spec.csv_comments = flag(*x, "csv.comments");
  }
  if (const auto *v = o.if_contains("json")) {
    const auto &c = *objp(*v, "json");
    reject_unknown(c, {"pretty"}, "json");
    if (const auto *x = c.if_contains("pretty")) spec.json_pretty = flag(*x, "json.pretty");
  }
  if (const auto *v = o.if_contains("parquet")) {
    const auto &c = *objp(*v, "parquet");
    reject_unknown(c, {"units_in_names"}, "parquet");
    if (const auto *x = c.if_contains("units_in_names")) spec.parquet_units_in_names = flag(*x, "parquet.units_in_names");
  }
  if (const auto *v = o.if_contains("tables")) {
    for (const auto &kv : *objp(*v, "tables")) {
      const std::string key(kv.key());
      const std::string w = "tables." + key;
      TableSpec ts;
      if (kv.value().is_bool()) {
        ts.include = kv.value().as_bool();
        spec.tables[key] = ts;
        continue;
      }
      const auto &t = *objp(kv.value(), w);
      reject_unknown(t, {"include", "columns", "exclude", "units", "labels", "formats",
                         "number_format", "filters", "drop_zero_rows", "every", "max_rows",
                         "file_name"},
                     w);
      if (const auto *x = t.if_contains("include")) ts.include = flag(*x, w + ".include");
      if (const auto *x = t.if_contains("columns")) ts.columns = strings(*x, w + ".columns");
      if (const auto *x = t.if_contains("exclude")) ts.exclude = strings(*x, w + ".exclude");
      if (const auto *x = t.if_contains("units")) {
        for (const auto &u : *objp(*x, w + ".units")) {
          ts.column_specs[std::string(u.key())].unit = str(u.value(), w + ".units." + std::string(u.key()));
          (void)find_unit(ts.column_specs[std::string(u.key())].unit);
        }
      }
      if (const auto *x = t.if_contains("labels")) {
        for (const auto &u : *objp(*x, w + ".labels")) {
          ts.column_specs[std::string(u.key())].label = str(u.value(), w + ".labels." + std::string(u.key()));
        }
      }
      if (const auto *x = t.if_contains("formats")) {
        for (const auto &u : *objp(*x, w + ".formats")) {
          ts.column_specs[std::string(u.key())].format =
              parse_format_object(u.value(), w + ".formats." + std::string(u.key()));
        }
      }
      if (const auto *x = t.if_contains("number_format")) ts.format = parse_format_object(*x, w + ".number_format");
      if (const auto *x = t.if_contains("filters")) {
        if (!x->is_array()) throw std::invalid_argument(w + ".filters: expected a list");
        for (std::size_t i = 0; i < x->as_array().size(); ++i) {
          const std::string fw = w + ".filters[" + std::to_string(i) + "]";
          const auto &f = *objp(x->as_array()[i], fw);
          reject_unknown(f, {"column", "min", "max"}, fw);
          RowFilter rf;
          if (const auto *c = f.if_contains("column")) rf.column = str(*c, fw + ".column");
          else throw std::invalid_argument(fw + ": 'column' is required");
          if (const auto *c = f.if_contains("min")) {
            auto b = bound_value(*c, fw + ".min");
            rf.min = b.first;
            rf.min_unit = b.second;
          }
          if (const auto *c = f.if_contains("max")) {
            auto b = bound_value(*c, fw + ".max");
            rf.max = b.first;
            rf.max_unit = b.second;
          }
          ts.filters.push_back(std::move(rf));
        }
      }
      if (const auto *x = t.if_contains("drop_zero_rows")) ts.drop_zero_rows = flag(*x, w + ".drop_zero_rows");
      if (const auto *x = t.if_contains("every")) ts.every = whole(*x, w + ".every", 1, 1e15);
      if (const auto *x = t.if_contains("max_rows")) ts.max_rows = whole(*x, w + ".max_rows", 0, 1e15);
      if (const auto *x = t.if_contains("file_name")) ts.file_name = file_name(*x, w + ".file_name");
      spec.tables[key] = std::move(ts);
    }
  }
  return spec;
}

std::string to_json(const OutputSpec &spec) {
  json::object o;
  json::array formats;
  for (auto f : spec.formats) formats.push_back(json::value(to_string(f)));
  o["formats"] = std::move(formats);
  if (!spec.units.empty()) {
    json::object u;
    for (const auto &kv : spec.units) u[to_string(kv.first)] = kv.second;
    o["units"] = std::move(u);
  }
  o["number_format"] = format_object(spec.format);
  {
    json::object c;
    c["delimiter"] = spec.csv_delimiter == '\t' ? std::string("tab") : std::string(1, spec.csv_delimiter);
    c["comments"] = spec.csv_comments;
    o["csv"] = std::move(c);
  }
  {
    json::object c;
    c["pretty"] = spec.json_pretty;
    o["json"] = std::move(c);
  }
  {
    json::object c;
    c["units_in_names"] = spec.parquet_units_in_names;
    o["parquet"] = std::move(c);
  }
  if (!spec.tables.empty()) {
    json::object tables;
    for (const auto &kv : spec.tables) {
      const TableSpec &ts = kv.second;
      json::object t;
      t["include"] = ts.include;
      if (!ts.columns.empty()) {
        json::array a;
        for (const auto &c : ts.columns) a.push_back(json::value(c));
        t["columns"] = std::move(a);
      }
      if (!ts.exclude.empty()) {
        json::array a;
        for (const auto &c : ts.exclude) a.push_back(json::value(c));
        t["exclude"] = std::move(a);
      }
      json::object units, labels, fmts;
      for (const auto &cs : ts.column_specs) {
        if (!cs.second.unit.empty()) units[cs.first] = cs.second.unit;
        if (!cs.second.label.empty()) labels[cs.first] = cs.second.label;
        if (cs.second.format) fmts[cs.first] = format_object(*cs.second.format);
      }
      if (!units.empty()) t["units"] = std::move(units);
      if (!labels.empty()) t["labels"] = std::move(labels);
      if (!fmts.empty()) t["formats"] = std::move(fmts);
      if (ts.format) t["number_format"] = format_object(*ts.format);
      if (!ts.filters.empty()) {
        json::array a;
        for (const auto &f : ts.filters) {
          json::object fo;
          fo["column"] = f.column;
          const auto b = [](double v, const std::string &u) -> json::value {
            if (u.empty()) return json::value(v);
            return json::value(format_number(v, NumberFormat{}) + " " + u);
          };
          if (f.min) fo["min"] = b(*f.min, f.min_unit);
          if (f.max) fo["max"] = b(*f.max, f.max_unit);
          a.push_back(std::move(fo));
        }
        t["filters"] = std::move(a);
      }
      if (ts.drop_zero_rows) t["drop_zero_rows"] = true;
      if (ts.every > 1) t["every"] = ts.every;
      if (ts.max_rows > 0) t["max_rows"] = ts.max_rows;
      if (!ts.file_name.empty()) t["file_name"] = ts.file_name;
      tables[kv.first] = std::move(t);
    }
    o["tables"] = std::move(tables);
  }
  return json::serialize(o);
}

} // namespace output
} // namespace amcset
