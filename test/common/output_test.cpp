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

/*! \file output_test.cpp
 *  \brief Tests of the output tables, units, specs and file writers.
 */

#include "output/table.h"
#include "output/text.h"
#include "output/units.h"
#include "output/writers.h"

#include <gtest/gtest.h>
#include "output/stream.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>

using namespace amcset::output;

namespace {

Table sample_table() {
  Table t;
  t.key = "range";
  t.title = "Sample";
  t.columns.push_back(Column::real("depth", "Depth", Dimension::length, "A", {100, 200, 300, 400}));
  t.columns.push_back(Column::real("ions", "Ions", Dimension::inverse_length, "1/cm", {0, 2e4, 0, 1e4}));
  t.columns.push_back(Column::real("loss", "Loss", Dimension::energy_per_length, "eV/A", {10, 20, 30, 40}));
  t.columns.push_back(Column::strings("tag", "Tag, with comma", {"a", "b\"q", "c", "d"}));
  return t;
}

std::filesystem::path temp_dir(const std::string &name) {
  auto p = std::filesystem::temp_directory_path() / ("amcset_output_test_" + name);
  std::filesystem::remove_all(p);
  std::filesystem::create_directories(p);
  return p;
}

} // namespace

TEST(OutputUnits, ConvertWithinADimension) {
  EXPECT_NEAR(conversion_factor("A", "nm"), 0.1, 1e-15);
  EXPECT_NEAR(conversion_factor("um", "A"), 1e4, 1e-9);
  EXPECT_NEAR(conversion_factor("eV/A", "keV/um"), 10.0, 1e-12);
  EXPECT_NEAR(conversion_factor("1/A", "1/cm"), 1e8, 1e-4);
  EXPECT_NEAR(conversion_factor("keV", "eV"), 1e3, 1e-9);
  EXPECT_NEAR(conversion_factor("deg", "rad"), 0.017453292519943295, 1e-15);
  EXPECT_NEAR(conversion_factor("1", "%"), 100.0, 1e-12);
  EXPECT_EQ(&find_unit("Angstrom"), &find_unit("A"));
  EXPECT_EQ(find_unit("(atoms/cm3)/(atoms/cm2)").name, "1/cm");
}

TEST(OutputUnits, RejectMismatchedDimensions) {
  EXPECT_THROW(conversion_factor("A", "eV"), std::invalid_argument);
  EXPECT_THROW(find_unit("keV", Dimension::length), std::invalid_argument);
  EXPECT_THROW(find_unit("parsec"), std::invalid_argument);
  for (auto d : {Dimension::dimensionless, Dimension::length, Dimension::inverse_length,
                 Dimension::energy, Dimension::energy_per_length, Dimension::inverse_energy,
                 Dimension::angle, Dimension::time}) {
    EXPECT_FALSE(units_of(d).empty());
    EXPECT_EQ(parse_dimension(to_string(d)), d);
    EXPECT_DOUBLE_EQ(base_unit(d).per_base, 1.0);
  }
}

TEST(OutputFormat, NumberNotations) {
  NumberFormat f;
  f.notation = NumberFormat::Notation::scientific;
  f.digits = 3;
  EXPECT_EQ(format_number(12345.678, f), "1.235E+04");
  f.uppercase = false;
  EXPECT_EQ(format_number(12345.678, f), "1.235e+04");
  f.notation = NumberFormat::Notation::fixed;
  f.digits = 2;
  EXPECT_EQ(format_number(3.14159, f), "3.14");
  f.notation = NumberFormat::Notation::engineering;
  f.digits = 2;
  f.uppercase = true;
  EXPECT_EQ(format_number(12345.678, f), "12.35E+03");
  EXPECT_EQ(format_number(0.00012, f), "120.00E-06");
  EXPECT_EQ(format_number(999999.0, f), "1.00E+06");
  f.notation = NumberFormat::Notation::shortest;
  EXPECT_EQ(format_number(0.1, f), "0.1");
  f.notation = NumberFormat::Notation::general;
  f.digits = 12;
  EXPECT_EQ(format_number(14.999999999999998, f), "15");
  EXPECT_EQ(format_number(std::nan(""), f), "nan");
}

TEST(OutputSpec, ShapesTables) {
  const Table t = sample_table();
  OutputSpec spec;
  spec.units[Dimension::length] = "nm";
  TableSpec ts;
  ts.columns = {"depth", "ions", "loss"};
  ts.column_specs["loss"].unit = "keV/um";
  ts.column_specs["ions"].label = "Ion density";
  ts.drop_zero_rows = false;
  RowFilter f;
  f.column = "depth";
  f.min = 15.0; // nm
  ts.filters.push_back(f);
  spec.tables["range"] = ts;
  const Table out = apply(t, spec);
  ASSERT_EQ(out.columns.size(), 3u);
  EXPECT_EQ(out.columns[0].unit, "nm");
  EXPECT_EQ(out.rows(), 3u);
  EXPECT_DOUBLE_EQ(out.columns[0].values[0], 20.0);
  EXPECT_DOUBLE_EQ(out.columns[2].values[0], 200.0); // 20 eV/A = 200 keV/um
  EXPECT_EQ(out.columns[1].label, "Ion density");

  // Zero rows and decimation; a filter bound with its own unit.
  spec.tables["range"].drop_zero_rows = true;
  spec.tables["range"].columns = {"depth", "ions"};
  spec.tables["range"].filters.clear();
  EXPECT_EQ(apply(t, spec).rows(), 2u);
  spec.tables["range"].drop_zero_rows = false;
  spec.tables["range"].every = 2;
  EXPECT_EQ(apply(t, spec).rows(), 2u);
  spec.tables["range"].every = 1;
  RowFilter g;
  g.column = "depth";
  g.max = 0.025;
  g.max_unit = "um";
  spec.tables["range"].filters = {g};
  EXPECT_EQ(apply(t, spec).rows(), 2u);

  spec.tables["range"].exclude = {"nope"};
  EXPECT_THROW(apply(t, spec), std::invalid_argument);
}

TEST(OutputSpec, JsonRoundTrip) {
  const char *text = R"({
    "formats": ["csv", "parquet"],
    "units": {"length": "um", "energy": "MeV"},
    "number_format": {"notation": "scientific", "digits": 4},
    "csv": {"delimiter": "tab", "comments": false},
    "tables": {
      "range": {"columns": ["depth", "ions"], "units": {"depth": "nm"},
                "filters": [{"column": "depth", "min": "100 nm", "max": 2}],
                "drop_zero_rows": true, "every": 2, "max_rows": 10,
                "number_format": "fixed", "file_name": "ranges"},
      "exits": false
    }
  })";
  const OutputSpec s = parse_output_spec(text);
  EXPECT_EQ(s.formats.size(), 2u);
  EXPECT_TRUE(s.formats.count(Format::parquet));
  EXPECT_EQ(s.units.at(Dimension::length), "um");
  EXPECT_EQ(s.format.notation, NumberFormat::Notation::scientific);
  EXPECT_EQ(s.csv_delimiter, '\t');
  EXPECT_FALSE(s.table("exits").include);
  const auto &r = s.table("range");
  EXPECT_EQ(r.filters.at(0).min_unit, "nm");
  EXPECT_EQ(*r.filters.at(0).max, 2.0);
  EXPECT_EQ(r.every, 2u);
  const OutputSpec back = parse_output_spec(to_json(s));
  EXPECT_EQ(to_json(back), to_json(s));
  EXPECT_THROW(parse_output_spec(R"({"units": {"length": "keV"}})"), std::invalid_argument);
  EXPECT_THROW(parse_output_spec(R"({"formats": ["xlsx"]})"), std::invalid_argument);
  EXPECT_THROW(parse_output_spec(R"({"bogus": 1})"), std::invalid_argument);
}

TEST(OutputWriters, CsvQuotesAndUnits) {
  const Table t = sample_table();
  Dataset d;
  d.title = "test";
  OutputSpec spec;
  spec.csv_comments = false;
  std::ostringstream o;
  write_csv(o, apply(t, spec), d, spec);
  const std::string s = o.str();
  EXPECT_EQ(s.substr(0, s.find('\n')), "Depth [A],Ions [1/cm],Loss [eV/A],\"Tag, with comma\"");
  EXPECT_NE(s.find("\"b\"\"q\""), std::string::npos);
}

TEST(OutputWriters, DatasetFilesAndParquetLayout) {
  Dataset d;
  d.title = "test run";
  d.metadata = {{"ion", "Fe"}};
  d.scalars.push_back(Scalar{"range_mean", "Range", Dimension::length, "A", 1234.5, ""});
  d.tables.push_back(sample_table());
  Table ints;
  ints.key = "list";
  ints.columns.push_back(Column::integer("ion", "Ion", {1, 2, 3}));
  ints.columns.push_back(Column::real("e", "E", Dimension::energy, "eV", {1.5, 2.5, 3.5}));
  d.tables.push_back(ints);
  OutputSpec spec;
  spec.formats = {Format::csv, Format::json, Format::parquet};
  spec.units[Dimension::length] = "nm";
  const auto dir = temp_dir("dataset");
  const auto files = write_dataset(dir, d, spec);
  EXPECT_EQ(files.size(), 5u);
  for (const char *f : {"csv/range.csv", "csv/list.csv", "parquet/range.parquet",
                        "parquet/list.parquet", "results.json"}) {
    EXPECT_TRUE(std::filesystem::exists(dir / f)) << f;
  }
  // JSON carries the converted scalar.
  std::ifstream js(dir / "results.json");
  std::stringstream jt;
  jt << js.rdbuf();
  EXPECT_NE(jt.str().find("\"value\": 123.45, \"unit\": \"nm\""), std::string::npos) << jt.str();
  // Parquet: magic at both ends, footer length inside the file.
  std::ifstream pq(dir / "parquet/range.parquet", std::ios::binary);
  std::string bytes((std::istreambuf_iterator<char>(pq)), std::istreambuf_iterator<char>());
  ASSERT_GT(bytes.size(), 12u);
  EXPECT_EQ(bytes.substr(0, 4), "PAR1");
  EXPECT_EQ(bytes.substr(bytes.size() - 4), "PAR1");
  std::uint32_t footer = 0;
  for (int i = 0; i < 4; ++i) {
    footer |= static_cast<std::uint32_t>(static_cast<unsigned char>(bytes[bytes.size() - 8 + i])) << (8 * i);
  }
  EXPECT_LT(footer + 12u, bytes.size());
  EXPECT_NE(bytes.find("amcset.units"), std::string::npos);
  EXPECT_NE(bytes.find("\"depth\":\"nm\""), std::string::npos);
}

// ---------------------------------------------------------------------------
// Fast number text (output/text.h, format_number_to): exactly printf's text
// ---------------------------------------------------------------------------

namespace {

//! Values that stress the conversions: random bit patterns over the whole
//! exponent range, exact decimal ties, zeros, denormals, extremes, integers.
std::vector<double> hostile_values(std::size_t random_count) {
  std::vector<double> v = {0.0, -0.0, 1.0, -1.0, 0.5, 1.5, 2.5, -2.5, 0.125, 0.375, 1.0625, 0.05, 0.15,
                           0.25, 0.35, 1e23, 9.5, 99.5, 999.5, 9999.5, 0.00015, 1.00005, 2.00005,
                           123456.5, 5e-324, -5e-324, 2.2250738585072014e-308, 1.7976931348623157e308,
                           -1.7976931348623157e308, 1e-300, 1e300, 65263.6, 1e15, 1e16, 1e17,
                           9.99995e10, 9.99995, 0.999995, 9.999949999e-5};
  std::mt19937_64 rng(20261123);
  for (std::size_t i = 0; i < random_count; ++i) {
    std::uint64_t bits = rng();
    double d = 0.0;
    std::memcpy(&d, &bits, sizeof d);
    if (std::isfinite(d)) v.push_back(d);
    // Values of the size the output files hold, also on short decimal grids
    // (many ties and near-ties).
    const double m = std::ldexp(static_cast<double>(rng() % 1000000), static_cast<int>(rng() % 60) - 40);
    v.push_back((rng() & 1) ? m : -m);
    v.push_back(static_cast<double>(rng() % 200000) / 16.0);
  }
  return v;
}

// The reference: printf itself, with formats that are string literals at
// every call site.
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
#endif
std::string printf_text(const char *format, int width, int precision, double v) {
  char buf[512];
  if (width >= 0) {
    std::snprintf(buf, sizeof buf, format, width, precision, v);
  } else {
    std::snprintf(buf, sizeof buf, format, precision, v);
  }
  return buf;
}
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

} // namespace

TEST(OutputText, SciAndFixedFieldsMatchPrintf) {
  using namespace amcset::output::text;
  char buf[max_field];
  int checked = 0;
  for (const double v : hostile_values(30000)) {
    for (const auto [w, d] : {std::pair{11, 4}, std::pair{12, 4}, std::pair{13, 6}, std::pair{10, 3},
                              std::pair{0, 0}, std::pair{0, 11}, std::pair{0, 17}}) {
      const std::string e(buf, put_sci(buf, v, w, d));
      ASSERT_EQ(e, printf_text("%*.*E", w, d, v)) << "value " << v << " width " << w << " precision " << d;
      const std::string el(buf, put_sci(buf, v, w, d, false));
      ASSERT_EQ(el, printf_text("%*.*e", w, d, v)) << "value " << v;
      if (std::abs(v) < 1e30) {
        const std::string f(buf, put_fixed(buf, v, w == 0 ? 0 : 9, d));
        ASSERT_EQ(f, printf_text("%*.*f", w == 0 ? 0 : 9, d, v)) << "value " << v << " precision " << d;
      }
      ++checked;
    }
  }
  // The largest values in fixed notation fit the field.
  const std::string big(buf, put_fixed(buf, -1.7976931348623157e308, 0, 30));
  EXPECT_EQ(big, printf_text("%*.*f", 0, 30, -1.7976931348623157e308));
  // Non-finite values are printf's own text.
  for (const double v : {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
    EXPECT_EQ(std::string(buf, put_sci(buf, v, 11, 4)), printf_text("%*.*E", 11, 4, v));
    EXPECT_EQ(std::string(buf, put_fixed(buf, v, 9, 6)), printf_text("%*.*f", 9, 6, v));
  }
  EXPECT_GT(checked, 100000);
}

TEST(OutputText, UnsignedFieldsMatchPrintf) {
  using namespace amcset::output::text;
  char buf[max_field];
  std::mt19937_64 rng(7);
  for (int i = 0; i < 20000; ++i) {
    const unsigned long long v = i < 100 ? static_cast<unsigned long long>(i) : rng() >> (rng() % 64);
    char ref[64];
    std::snprintf(ref, sizeof ref, "%07llu", v);
    ASSERT_EQ(std::string(buf, put_unsigned(buf, v, 7, '0')), ref);
    std::snprintf(ref, sizeof ref, "%6llu", v);
    ASSERT_EQ(std::string(buf, put_unsigned(buf, v, 6)), ref);
  }
}

TEST(OutputText, FormatNumberMatchesPrintf) {
  using N = NumberFormat::Notation;
  char buf[text::max_field];
  for (const double v : hostile_values(20000)) {
    for (const int d : {0, 1, 4, 6, 12, 17}) {
      for (const bool up : {true, false}) {
        NumberFormat f;
        f.digits = d;
        f.uppercase = up;
        f.notation = N::general;
        ASSERT_EQ(std::string(buf, format_number_to(buf, v, f)),
                  printf_text(up ? "%.*G" : "%.*g", -1, std::max(d, 1), v))
            << "value " << v << " digits " << d;
        f.notation = N::scientific;
        ASSERT_EQ(std::string(buf, format_number_to(buf, v, f)), printf_text(up ? "%.*E" : "%.*e", -1, d, v))
            << "value " << v << " digits " << d;
        if (std::abs(v) < 1e30) {
          f.notation = N::fixed;
          ASSERT_EQ(std::string(buf, format_number_to(buf, v, f)), printf_text("%.*f", -1, d, v))
              << "value " << v << " digits " << d;
        }
        // The string form is the same text.
        ASSERT_EQ(format_number(v, f), std::string(buf, format_number_to(buf, v, f)));
      }
    }
  }
  NumberFormat f;
  EXPECT_EQ(std::string(buf, format_number_to(buf, std::numeric_limits<double>::quiet_NaN(), f)), "nan");
  EXPECT_EQ(std::string(buf, format_number_to(buf, -std::numeric_limits<double>::infinity(), f)), "-inf");
}

TEST(OutputText, BufferedWriterWritesEverythingInOrder) {
  std::ostringstream o;
  {
    amcset::output::text::BufferedWriter w(o, 1000); // small: many flushes
    std::string expected;
    for (int i = 0; i < 5000; ++i) {
      char *q = w.reserve(64);
      q = amcset::output::text::put_unsigned(q, static_cast<unsigned long long>(i), 5, '0');
      *q++ = '\n';
      w.commit(q);
      char ref[16];
      std::snprintf(ref, sizeof ref, "%05d\n", i);
      expected += ref;
    }
    const std::string long_text(70000, 'x'); // longer than the buffer
    w.put(long_text);
    expected += long_text;
    w.finish();
    EXPECT_EQ(o.str(), expected);
  }
}

// The streaming writer: CSV with one header line of keys, Parquet with
// several row groups (its footer lists them; the bytes end in "PAR1").
TEST(OutputStream, CsvAndParquetPieceByPiece) {
  const auto dir = std::filesystem::temp_directory_path() / "amcset_test_stream";
  std::filesystem::remove_all(dir);
  Table schema;
  schema.key = "pieces";
  schema.columns.push_back(Column::integer("n", "N"));
  schema.columns.push_back(Column::real("x_A", "x", Dimension::length, "A"));
  schema.columns.push_back(Column::strings("kind", "Kind"));
  const auto piece = [&](std::size_t first, std::size_t count) {
    Table t = schema;
    for (std::size_t i = first; i < first + count; ++i) {
      t.columns[0].values.push_back(static_cast<double>(i));
      t.columns[1].values.push_back(0.1 * static_cast<double>(i));
      t.columns[2].text.push_back(i % 2 ? "odd" : "even");
    }
    return t;
  };
  {
    StreamTableWriter w(dir / "a.csv", schema);
    w.append(piece(0, 3));
    w.append(piece(3, 2));
    w.finish();
    EXPECT_EQ(w.rows_written(), 5u);
  }
  std::ifstream in(dir / "a.csv");
  std::stringstream ss;
  ss << in.rdbuf();
  EXPECT_EQ(ss.str(), "n,x_A,kind\n0,0,even\n1,0.1,odd\n2,0.2,even\n3,0.30000000000000004,odd\n4,0.4,even\n");
  {
    StreamTableWriter w(dir / "b.parquet", schema, {{"note", "test"}});
    for (int k = 0; k < 3; ++k) w.append(piece(static_cast<std::size_t>(k) * 200000, 200000)); // > one row group
    w.finish();
  }
  std::ifstream pin(dir / "b.parquet", std::ios::binary);
  std::stringstream ps;
  ps << pin.rdbuf();
  const std::string bytes = ps.str();
  ASSERT_GT(bytes.size(), 600000u * 8u);
  EXPECT_EQ(bytes.substr(0, 4), "PAR1");
  EXPECT_EQ(bytes.substr(bytes.size() - 4), "PAR1");
  EXPECT_NE(bytes.find("amcset.note"), std::string::npos);
  EXPECT_NE(bytes.find("\"x_A\":\"A\""), std::string::npos);
  EXPECT_THROW(StreamTableWriter(dir / "c.txt", schema), std::invalid_argument);
  StreamTableWriter w(dir / "d.csv", schema);
  Table wrong = schema;
  wrong.columns.pop_back();
  EXPECT_THROW(w.append(wrong), std::invalid_argument);
}
