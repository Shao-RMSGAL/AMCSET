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

/*! \file run_config.cpp
 *  \brief JSON scenario files, unit parsing and engine construction.
 */

#include "run_config.h"

#include "amcset_common.h"
#include "crystal/cell_io.h"
#include "crystal/grain_io.h"
#include "element_data.h"
#include "isotope_data.h"
#include "output/table.h"
#include "trim/trim_dat.h"

// Boost.JSON is used header-only: its implementation is compiled into this
// translation unit (and only this one).
#include <boost/json/src.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace amcset {
namespace config {

using namespace common;
namespace json = boost::json;

// ---------------------------------------------------------------------------
// Quantity parsing
// ---------------------------------------------------------------------------

namespace {

std::string lower(std::string_view s) {
  std::string out;
  for (char c : s) out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  return out;
}

std::string trim(std::string_view s) {
  std::size_t a = 0;
  std::size_t b = s.size();
  while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
  while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
  return std::string(s.substr(a, b - a));
}

//! Split "1.5 MeV" / "1.5MeV" into (1.5, "MeV").
std::pair<double, std::string> split_quantity(std::string_view text, std::string_view what) {
  const std::string t = trim(text);
  if (t.empty()) throw std::invalid_argument(std::string(what) + ": empty value");
  std::size_t used = 0;
  double v = 0.0;
  try {
    v = std::stod(t, &used);
  } catch (const std::exception &) {
    throw std::invalid_argument(std::string(what) + ": '" + t + "' does not start with a number");
  }
  if (!std::isfinite(v)) throw std::invalid_argument(std::string(what) + ": not a finite number");
  return {v, trim(std::string_view(t).substr(used))};
}

template <class Q>
Q parse_with(std::string_view text, std::string_view what,
             const std::map<std::string, Q> &units, const char *expected) {
  const auto [v, unit] = split_quantity(text, what);
  if (unit.empty()) {
    throw std::invalid_argument(std::string(what) + ": '" + std::string(text) +
                                "' needs a unit (" + expected + ")");
  }
  const auto it = units.find(lower(unit));
  if (it == units.end()) {
    throw std::invalid_argument(std::string(what) + ": unknown unit '" + unit +
                                "' (expected " + expected + ")");
  }
  return v * it->second;
}

const std::map<std::string, energy_quantity> energy_units = {
    {"ev", electron_volt},          {"kev", kilo_electron_volt},
    {"mev", mega_electron_volt},    {"gev", giga_electron_volt},
    {"j", joule}};
const std::map<std::string, length_quantity> length_units = {
    {"a", angstrom},       {"angstrom", angstrom}, {"angstroms", angstrom},
    {"\xC3\x85", angstrom}, {"nm", nanometer},       {"um", micrometer},
    {"micron", micrometer}, {"\xC2\xB5m", micrometer}, {"mm", 1e-3 * meter},
    {"cm", centimeter},    {"m", meter}};
const std::map<std::string, mass_quantity> mass_units = {
    {"amu", atomic_mass_unit}, {"u", atomic_mass_unit}, {"da", atomic_mass_unit},
    {"kg", kilogram},          {"g", gram}};
const std::map<std::string, mass_density_quantity> density_units = {
    {"g/cm3", gram_per_cubic_centimeter},  {"g/cm^3", gram_per_cubic_centimeter},
    {"g/cc", gram_per_cubic_centimeter},   {"g/ml", gram_per_cubic_centimeter},
    {"kg/m3", kg_per_cubic_meter},         {"kg/m^3", kg_per_cubic_meter}};
const std::map<std::string, angle_quantity> angle_units = {
    {"deg", degree}, {"degree", degree}, {"degrees", degree}, {"\xC2\xB0", degree},
    {"rad", radian}, {"radian", radian}, {"radians", radian}};

std::string number(double v) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.10g", v);
  return buf;
}

} // namespace

energy_quantity parse_energy(std::string_view text, std::string_view what) {
  return parse_with(text, what, energy_units, "eV, keV, MeV, GeV");
}

std::vector<energy_quantity> parse_energy_list(const std::string &text, const std::string &what) {
  std::vector<std::string> parts;
  std::string cur;
  for (char ch : text) {
    if (ch == ',') {
      parts.push_back(cur);
      cur.clear();
    } else {
      cur += ch;
    }
  }
  parts.push_back(cur);
  const auto trim_ws = [](std::string s) {
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
    while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
    return s;
  };
  // A bare number takes the unit of the next entry that has one.
  std::vector<energy_quantity> out(parts.size());
  std::string unit;
  for (std::size_t k = parts.size(); k-- > 0;) {
    std::string p = trim_ws(parts[k]);
    if (p.empty()) throw std::invalid_argument(what + ": an empty entry in '" + text + "'");
    // The unit is what follows the number (strtod reads "1e3" as a number).
    char *end = nullptr;
    std::strtod(p.c_str(), &end);
    const std::string rest = trim_ws(std::string(end));
    if (!rest.empty()) {
      unit = rest;
    } else {
      if (unit.empty()) throw std::invalid_argument(what + ": '" + p + "' has no unit (e.g. \"0.5,1,2 keV\")");
      p += " " + unit;
    }
    out[k] = parse_energy(p, what);
  }
  return out;
}
length_quantity parse_length(std::string_view text, std::string_view what) {
  return parse_with(text, what, length_units, "A, nm, um, mm, cm, m");
}
mass_quantity parse_mass(std::string_view text, std::string_view what) {
  return parse_with(text, what, mass_units, "amu, g, kg");
}
mass_density_quantity parse_density(std::string_view text, std::string_view what) {
  return parse_with(text, what, density_units, "g/cm3, kg/m3");
}
angle_quantity parse_angle(std::string_view text, std::string_view what) {
  return parse_with(text, what, angle_units, "deg, rad");
}
temperature_quantity parse_temperature(std::string_view text, std::string_view what) {
  const auto [v, unit] = split_quantity(text, what);
  const std::string u = lower(unit);
  if (u == "k" || u == "kelvin") return v * kelvin;
  // Celsius is an offset, not a factor.
  if (u == "c" || u == "degc" || u == "celsius" || u == "\xC2\xB0""c") return (v + 273.15) * kelvin;
  if (u.empty()) {
    throw std::invalid_argument(std::string(what) + ": '" + std::string(text) + "' needs a unit (K, C)");
  }
  throw std::invalid_argument(std::string(what) + ": unknown unit '" + unit + "' (expected K, C)");
}

std::string format_energy(energy_quantity e) {
  const double ev = in_units(e, electron_volt);
  const double a = std::abs(ev);
  if (a >= 1e9) return number(ev / 1e9) + " GeV";
  if (a >= 1e6) return number(ev / 1e6) + " MeV";
  if (a >= 1e3) return number(ev / 1e3) + " keV";
  return number(ev) + " eV";
}
std::string format_length(length_quantity l) {
  const double a = in_units(l, angstrom);
  if (std::abs(a) >= 1e7) return number(a / 1e7) + " mm";
  if (std::abs(a) >= 1e4) return number(a / 1e4) + " um";
  return number(a) + " A";
}
std::string format_mass(mass_quantity m) { return number((m / atomic_mass_unit).value()) + " amu"; }
std::string format_density(mass_density_quantity d) {
  return number(in_units(d, gram_per_cubic_centimeter)) + " g/cm3";
}
std::string format_angle(angle_quantity a) { return number(in_units(a, degree)) + " deg"; }
std::string format_temperature(temperature_quantity t) { return number(in_units(t, kelvin)) + " K"; }

// ---------------------------------------------------------------------------
// JSON reading helpers
// ---------------------------------------------------------------------------

namespace {

class Reader {
public:
  Reader(const json::object &o, std::string path) : o_(o), path_(std::move(path)) {}

  bool has(std::string_view k) const { return o_.contains(k); }
  const std::string &path() const { return path_; }
  std::string where(std::string_view k) const {
    return path_.empty() ? std::string(k) : path_ + "." + std::string(k);
  }
  const json::value *get(std::string_view k) const {
    const auto it = o_.find(k);
    return it == o_.end() ? nullptr : &it->value();
  }
  std::string str(std::string_view k, std::string def = {}) const {
    const auto *v = get(k);
    if (!v || v->is_null()) return def;
    if (!v->is_string()) throw std::invalid_argument(where(k) + ": expected a string");
    return std::string(v->as_string());
  }
  //! A finite number (a JSON number, or a string that is one number and
  //! nothing else: "12abc", "nan" and "inf" are errors).
  double num(std::string_view k, double def) const {
    const auto *v = get(k);
    if (!v || v->is_null()) return def;
    double d = 0.0;
    bool ok = true;
    if (v->is_double()) {
      d = v->as_double();
    } else if (v->is_int64()) {
      d = static_cast<double>(v->as_int64());
    } else if (v->is_uint64()) {
      d = static_cast<double>(v->as_uint64());
    } else if (v->is_string()) {
      const std::string s = trim(v->as_string());
      try {
        std::size_t used = 0;
        d = std::stod(s, &used);
        ok = used == s.size();
      } catch (const std::exception &) {
        ok = false;
      }
    } else {
      ok = false;
    }
    if (!ok) throw std::invalid_argument(where(k) + ": expected a number");
    if (!std::isfinite(d)) throw std::invalid_argument(where(k) + ": not a finite number");
    return d;
  }
  /*!
   * A whole number in [0, \p max] (checked before the integer conversion:
   * out-of-range doubles have no defined integer value, and narrowing to a
   * smaller type would wrap silently, e.g. a Z of 2^32 + 26).
   */
  std::uint64_t count(std::string_view k, std::uint64_t def,
                      std::uint64_t max = std::uint64_t{1} << 53) const {
    const double d = num(k, static_cast<double>(def));
    if (!(d >= 0.0) || d != std::floor(d)) {
      throw std::invalid_argument(where(k) + ": expected a non-negative integer");
    }
    if (d > static_cast<double>(max)) {
      throw std::invalid_argument(where(k) + ": at most " + std::to_string(max));
    }
    return static_cast<std::uint64_t>(d);
  }
  bool flag(std::string_view k, bool def) const {
    const auto *v = get(k);
    if (!v || v->is_null()) return def;
    if (v->is_bool()) return v->as_bool();
    if (v->is_string()) {
      const auto s = lower(v->as_string());
      if (s == "true" || s == "yes" || s == "on" || s == "1") return true;
      if (s == "false" || s == "no" || s == "off" || s == "0") return false;
    }
    if (v->is_int64()) return v->as_int64() != 0;
    throw std::invalid_argument(where(k) + ": expected true or false");
  }
  //! A quantity given as "<number> <unit>" under \p k, or as a plain number
  //! under \p k_unit (e.g. "energy_keV") in unit \p unit.
  template <class Q, class F>
  Q quantity(std::string_view k, std::string_view k_unit, Q unit, Q def, F parse) const {
    if (const auto *v = get(k); v && !v->is_null()) {
      if (v->is_string()) return parse(std::string_view(v->as_string()), where(k));
      throw std::invalid_argument(where(k) + ": expected a string with a unit, e.g. \"1 MeV\"");
    }
    if (!k_unit.empty() && has(k_unit)) return num(k_unit, 0.0) * unit;
    return def;
  }
  const json::object &object(std::string_view k) const {
    const auto *v = get(k);
    if (!v || !v->is_object()) {
      static const json::object empty;
      if (v && !v->is_null()) throw std::invalid_argument(where(k) + ": expected an object");
      return empty;
    }
    return v->as_object();
  }
  const json::array *array(std::string_view k) const {
    const auto *v = get(k);
    if (!v || v->is_null()) return nullptr;
    if (!v->is_array()) throw std::invalid_argument(where(k) + ": expected an array");
    return &v->as_array();
  }
  void reject_unknown(const std::vector<std::string> &known) const {
    for (const auto &kv : o_) {
      if (std::find(known.begin(), known.end(), std::string(kv.key())) == known.end()) {
        throw std::invalid_argument(where(kv.key()) + ": unknown key");
      }
    }
  }
  void reject_unknown(std::initializer_list<std::string_view> known) const {
    for (const auto &kv : o_) {
      bool ok = false;
      for (auto k : known) ok = ok || kv.key() == k;
      if (!ok) {
        throw std::invalid_argument(where(kv.key()) + ": unknown key");
      }
    }
  }

private:
  const json::object &o_;
  std::string path_;
};

unsigned element_z(const Reader &r) {
  if (r.has("element")) {
    const auto name = r.str("element");
    const auto z = Element::find_z(name);
    if (!z) throw std::invalid_argument(r.where("element") + ": unknown element '" + name + "'");
    if (r.has("z") && r.count("z", 0, Element::max_z) != *z) {
      throw std::invalid_argument(r.where("z") + ": does not match element '" + name + "'");
    }
    return *z;
  }
  if (r.has("z")) {
    const auto z = static_cast<unsigned>(r.count("z", 0, Element::max_z));
    if (z < 1) throw std::invalid_argument(r.where("z") + ": must be 1 to " + std::to_string(Element::max_z));
    return z;
  }
  throw std::invalid_argument(r.where("element") + ": missing (give \"element\" or \"z\")");
}

std::string resolve(const std::string &p, const std::filesystem::path &base) {
  if (p.empty() || base.empty()) return p;
  const std::filesystem::path path(p);
  return path.is_absolute() ? p : (base / path).lexically_normal().string();
}

// ---------------------------------------------------------------------------
// Crystals (a layer's "crystal" object; not TRIM)
// ---------------------------------------------------------------------------

//! Run \p f; library errors get the key \p where in front.
template <class F>
auto with_where(const std::string &where, F &&f) {
  try {
    return f();
  } catch (const std::invalid_argument &e) {
    throw std::invalid_argument(where + ": " + e.what());
  } catch (const std::runtime_error &e) {
    throw std::runtime_error(where + ": " + e.what());
  }
}

double json_number(const json::value &v, const std::string &where) {
  double d = 0.0;
  if (v.is_double()) {
    d = v.as_double();
  } else if (v.is_int64()) {
    d = static_cast<double>(v.as_int64());
  } else if (v.is_uint64()) {
    d = static_cast<double>(v.as_uint64());
  } else {
    throw std::invalid_argument(where + ": expected a number");
  }
  if (!std::isfinite(d)) throw std::invalid_argument(where + ": not a finite number");
  return d;
}

//! An array of \p min to \p max numbers.
std::vector<double> number_list(const json::array *a, const std::string &where, std::size_t min, std::size_t max) {
  const std::string expected =
      min == max ? std::to_string(min) : std::to_string(min) + " to " + std::to_string(max);
  if (!a) throw std::invalid_argument(where + ": expected an array of " + expected + " numbers");
  if (a->size() < min || a->size() > max) {
    throw std::invalid_argument(where + ": expected " + expected + " numbers, got " + std::to_string(a->size()));
  }
  std::vector<double> out;
  for (std::size_t i = 0; i < a->size(); ++i) out.push_back(json_number((*a)[i], where + "." + std::to_string(i)));
  return out;
}

//! Whole numbers as JSON integers (Miller indices read "1", not "1.0").
json::value index_value(double d) {
  if (d == std::floor(d) && std::abs(d) < 1e15) return json::value(static_cast<std::int64_t>(d));
  return json::value(d);
}

/*!
 * The elements on a site or sublattice: "Fe"; {"Fe": 0.9, "Cr": 0.1};
 * {"element": "Fe", "fraction": 0.9}; or an array of the last form (or of
 * symbols). An empty array is an empty site.
 */
std::vector<crystal::Occupant> parse_occupants(const json::value &v, const std::string &where) {
  std::vector<crystal::Occupant> out;
  const auto add = [&](const std::string &name, double fraction, const std::string &w) {
    const auto z = Element::find_z(name);
    if (!z) throw std::invalid_argument(w + ": unknown element '" + name + "'");
    if (!(fraction > 0.0 && fraction <= 1.0)) {
      throw std::invalid_argument(w + ": the fraction of " + name + " must be above 0 and at most 1, got " + number(fraction));
    }
    for (const auto &o : out) {
      if (o.z == *z) throw std::invalid_argument(w + ": " + name + " is listed twice");
    }
    out.push_back({*z, fraction});
  };
  const auto add_object = [&](const json::object &o, const std::string &w) {
    if (o.contains("element") || o.contains("z")) {
      const Reader r(o, w);
      r.reject_unknown({"element", "z", "fraction"});
      add(std::string(Element::by_z(element_z(r)).symbol()), r.num("fraction", 1.0), w);
    } else {
      if (o.empty()) throw std::invalid_argument(w + ": expected {\"Fe\": 0.9, ...} or {\"element\": \"Fe\", \"fraction\": 0.9}");
      for (const auto &kv : o) add(std::string(kv.key()), json_number(kv.value(), w + "." + std::string(kv.key())), w);
    }
  };
  if (v.is_string()) {
    add(std::string(v.as_string()), 1.0, where);
  } else if (v.is_object()) {
    add_object(v.as_object(), where);
  } else if (v.is_array()) {
    const auto &a = v.as_array();
    for (std::size_t i = 0; i < a.size(); ++i) {
      const std::string w = where + "." + std::to_string(i);
      if (a[i].is_string()) {
        add(std::string(a[i].as_string()), 1.0, w);
      } else if (a[i].is_object()) {
        add_object(a[i].as_object(), w);
      } else {
        throw std::invalid_argument(w + ": expected an element symbol or {\"element\": ..., \"fraction\": ...}");
      }
    }
  } else {
    throw std::invalid_argument(where + ": expected an element symbol, an object or an array");
  }
  double sum = 0.0;
  for (const auto &o : out) sum += o.fraction;
  if (sum > 1.0 + 1e-9) throw std::invalid_argument(where + ": the fractions add up to " + number(sum) + ", more than 1");
  return out;
}

CrystalDescription parse_crystal_object(const Reader &cr, const std::filesystem::path &base_dir) {
  const std::string &where = cr.path();
  cr.reject_unknown({"name", "structure", "a", "a_A", "c_over_a", "u", "sublattices", "interstitials", "file",
                     "format", "elements", "lammps_types", "atom_style", "lattice", "sites", "orientation",
                     "temperature", "temperature_K", "vibrations", "polycrystal"});
  const auto pl = [](std::string_view s, const std::string &w) { return parse_length(s, w); };
  const auto pa = [](std::string_view s, const std::string &w) { return parse_angle(s, w); };
  const auto pt = [](std::string_view s, const std::string &w) { return parse_temperature(s, w); };
  const bool by_structure = cr.has("structure");
  const bool by_file = cr.has("file");
  const bool by_sites = cr.has("lattice") || cr.has("sites");
  if (static_cast<int>(by_structure) + static_cast<int>(by_file) + static_cast<int>(by_sites) != 1) {
    throw std::invalid_argument(where + ": give exactly one of \"structure\" (a built-in structure), \"file\" (a "
                                        "cell file) or \"lattice\" with \"sites\"");
  }
  const auto only_with = [&](std::initializer_list<const char *> keys, const char *owner) {
    for (const char *k : keys) {
      if (cr.has(k)) throw std::invalid_argument(cr.where(k) + ": only with \"" + owner + "\"");
    }
  };
  CrystalDescription out;
  crystal::UnitCell &cell = out.crystal.cell;
  if (by_structure) {
    only_with({"format", "elements", "lammps_types", "atom_style"}, "file");
    crystal::StructureSpec spec;
    spec.name = cr.str("structure");
    const bool has_a = cr.has("a") || cr.has("a_A");
    const length_quantity a = cr.quantity("a", "a_A", angstrom, 0.0 * angstrom, pl);
    if (has_a && !(a > 0.0 * angstrom)) throw std::invalid_argument(cr.where("a") + ": must be positive");
    spec.a_A = in_units(a, angstrom);
    out.lattice_from_density = !has_a;
    if (cr.has("c_over_a")) {
      spec.c_over_a = cr.num("c_over_a", 0.0);
      if (!(spec.c_over_a > 0.0)) throw std::invalid_argument(cr.where("c_over_a") + ": must be positive");
    }
    spec.u = cr.num("u", spec.u);
    const json::array *subs = cr.array("sublattices");
    if (!subs) {
      throw std::invalid_argument(cr.where("sublattices") + ": missing (the elements of each sublattice, e.g. "
                                                            "[\"Fe\"] or [{\"Fe\": 0.9, \"Cr\": 0.1}])");
    }
    for (std::size_t i = 0; i < subs->size(); ++i) {
      spec.sublattices.push_back(parse_occupants((*subs)[i], cr.where("sublattices") + "." + std::to_string(i)));
    }
    if (const json::array *ins = cr.array("interstitials")) {
      for (std::size_t i = 0; i < ins->size(); ++i) {
        const std::string iw = cr.where("interstitials") + "." + std::to_string(i);
        if (!(*ins)[i].is_object()) throw std::invalid_argument(iw + ": expected an object");
        const Reader ir((*ins)[i].as_object(), iw);
        ir.reject_unknown({"sites", "occupants"});
        crystal::InterstitialSpec is;
        is.sites = lower(ir.str("sites"));
        if (is.sites.empty()) throw std::invalid_argument(ir.where("sites") + ": missing (octahedral or tetrahedral)");
        if (const auto *occ = ir.get("occupants")) is.occupants = parse_occupants(*occ, ir.where("occupants"));
        spec.interstitials.push_back(std::move(is));
      }
    }
    cell = with_where(where, [&] { return crystal::build_structure(spec); });
  } else if (by_file) {
    only_with({"a", "a_A", "c_over_a", "u", "sublattices", "interstitials"}, "structure");
    const std::string file = resolve(cr.str("file"), base_dir);
    if (file.empty()) throw std::invalid_argument(cr.where("file") + ": empty");
    std::optional<crystal::CellFormat> format;
    if (cr.has("format")) {
      format = with_where(cr.where("format"), [&] { return crystal::parse_cell_format(cr.str("format")); });
    }
    crystal::CellReadOptions opt;
    if (const json::array *els = cr.array("elements")) {
      for (std::size_t i = 0; i < els->size(); ++i) {
        const std::string ew = cr.where("elements") + "." + std::to_string(i);
        if (!(*els)[i].is_string()) throw std::invalid_argument(ew + ": expected an element symbol");
        const auto z = Element::find_z(std::string((*els)[i].as_string()));
        if (!z) throw std::invalid_argument(ew + ": unknown element '" + std::string((*els)[i].as_string()) + "'");
        opt.elements.push_back(*z);
      }
    }
    if (cr.has("lammps_types")) {
      for (const auto &kv : cr.object("lammps_types")) {
        const std::string tw = cr.where("lammps_types") + "." + std::string(kv.key());
        int type = 0;
        try {
          std::size_t used = 0;
          type = std::stoi(std::string(kv.key()), &used);
          if (used != kv.key().size() || type < 1) throw std::invalid_argument("type");
        } catch (const std::exception &) {
          throw std::invalid_argument(tw + ": the key must be an atom type number (1, 2, ...)");
        }
        if (!kv.value().is_string()) throw std::invalid_argument(tw + ": expected an element symbol");
        const auto z = Element::find_z(std::string(kv.value().as_string()));
        if (!z) throw std::invalid_argument(tw + ": unknown element '" + std::string(kv.value().as_string()) + "'");
        opt.lammps_types[type] = *z;
      }
    }
    opt.lammps_atom_style = lower(cr.str("atom_style", ""));
    cell = with_where(where, [&] { return crystal::read_cell(file, format, opt); });
  } else {
    only_with({"a", "a_A", "c_over_a", "u", "sublattices", "interstitials"}, "structure");
    only_with({"format", "elements", "lammps_types", "atom_style"}, "file");
    if (!cr.has("lattice")) throw std::invalid_argument(cr.where("lattice") + ": missing");
    const Reader lr(cr.object("lattice"), cr.where("lattice"));
    lr.reject_unknown({"a", "b", "c", "alpha", "beta", "gamma", "vectors_A"});
    if (lr.has("vectors_A")) {
      for (const char *k : {"a", "b", "c", "alpha", "beta", "gamma"}) {
        if (lr.has(k)) throw std::invalid_argument(lr.where(k) + ": give vectors_A or the cell parameters, not both");
      }
      const json::array *rows = lr.array("vectors_A");
      if (!rows || rows->size() != 3) throw std::invalid_argument(lr.where("vectors_A") + ": expected three vectors");
      crystal::Mat3 m;
      for (std::size_t i = 0; i < 3; ++i) {
        const std::string rw = lr.where("vectors_A") + "." + std::to_string(i);
        if (!(*rows)[i].is_array()) throw std::invalid_argument(rw + ": expected [x, y, z] in Angstrom");
        const auto v = number_list(&(*rows)[i].as_array(), rw, 3, 3);
        m[static_cast<int>(i)] = crystal::Vec3{v[0], v[1], v[2]};
      }
      cell.lattice = with_where(lr.where("vectors_A"), [&] { return crystal::Lattice::from_vectors_A(m); });
    } else {
      if (!lr.has("a")) throw std::invalid_argument(lr.where("a") + ": missing (or give vectors_A)");
      const length_quantity a = lr.quantity("a", "", angstrom, 0.0 * angstrom, pl);
      const length_quantity b = lr.quantity("b", "", angstrom, a, pl);
      const length_quantity c = lr.quantity("c", "", angstrom, a, pl);
      const angle_quantity alpha = lr.quantity("alpha", "", degree, 90.0 * degree, pa);
      const angle_quantity beta = lr.quantity("beta", "", degree, 90.0 * degree, pa);
      const angle_quantity gamma = lr.quantity("gamma", "", degree, 90.0 * degree, pa);
      cell.lattice = with_where(lr.path(), [&] { return crystal::Lattice::from_parameters(a, b, c, alpha, beta, gamma); });
    }
    const json::array *sites = cr.array("sites");
    if (!sites || sites->empty()) throw std::invalid_argument(cr.where("sites") + ": at least one site is required");
    for (std::size_t i = 0; i < sites->size(); ++i) {
      const std::string sw = cr.where("sites") + "." + std::to_string(i);
      if (!(*sites)[i].is_object()) throw std::invalid_argument(sw + ": expected an object");
      const Reader sr((*sites)[i].as_object(), sw);
      sr.reject_unknown({"label", "position", "kind", "occupants", "element", "fraction"});
      crystal::Site site;
      const auto p = number_list(sr.array("position"), sr.where("position"), 3, 3);
      site.position = crystal::wrap_fractional(crystal::Vec3{p[0], p[1], p[2]});
      site.kind = with_where(sr.where("kind"), [&] { return crystal::parse_site_kind(lower(sr.str("kind", "lattice"))); });
      if (sr.has("occupants")) {
        if (sr.has("element") || sr.has("fraction")) {
          throw std::invalid_argument(sw + ": give \"occupants\" or \"element\" (with \"fraction\"), not both");
        }
        site.occupants = parse_occupants(*sr.get("occupants"), sr.where("occupants"));
      } else if (sr.has("element")) {
        json::object one;
        one["element"] = sr.str("element");
        one["fraction"] = sr.num("fraction", 1.0);
        site.occupants = parse_occupants(json::value(std::move(one)), sw);
      }
      site.label = sr.str("label", "");
      if (site.label.empty()) {
        site.label = (site.occupants.empty() ? std::string("X")
                                             : std::string(Element::by_z(site.occupants.front().z).symbol())) +
                     std::to_string(i + 1);
      }
      cell.sites.push_back(std::move(site));
    }
  }
  cell.name = cr.str("name", cell.name);

  if (cr.has("orientation")) {
    const Reader orr(cr.object("orientation"), cr.where("orientation"));
    orr.reject_unknown({"surface", "reference", "euler_deg", "twist", "twist_deg"});
    auto &o = out.crystal.orientation;
    if (orr.has("euler_deg")) {
      if (orr.has("surface") || orr.has("reference")) {
        throw std::invalid_argument(orr.where("euler_deg") + ": give euler_deg or surface/reference, not both");
      }
      const auto e = number_list(orr.array("euler_deg"), orr.where("euler_deg"), 3, 3);
      o.kind = crystal::OrientationKind::euler;
      o.euler_deg = {e[0], e[1], e[2]};
    } else {
      if (orr.has("surface")) o.surface = number_list(orr.array("surface"), orr.where("surface"), 3, 4);
      if (orr.has("reference")) o.reference = number_list(orr.array("reference"), orr.where("reference"), 3, 4);
    }
    o.twist = orr.quantity("twist", "twist_deg", degree, 0.0 * degree, pa);
  }
  out.crystal.temperature = cr.quantity("temperature", "temperature_K", kelvin, out.crystal.temperature, pt);
  if (const json::array *vs = cr.array("vibrations")) {
    for (std::size_t i = 0; i < vs->size(); ++i) {
      const std::string vw = cr.where("vibrations") + "." + std::to_string(i);
      if (!(*vs)[i].is_object()) throw std::invalid_argument(vw + ": expected an object");
      const Reader vr((*vs)[i].as_object(), vw);
      vr.reject_unknown({"element", "z", "debye_temperature", "debye_temperature_K", "rms", "rms_A"});
      crystal::Vibration v;
      v.z = element_z(vr);
      if (vr.has("debye_temperature") || vr.has("debye_temperature_K")) {
        v.debye_temperature = vr.quantity("debye_temperature", "debye_temperature_K", kelvin, 0.0 * kelvin, pt);
      }
      if (vr.has("rms") || vr.has("rms_A")) v.rms = vr.quantity("rms", "rms_A", angstrom, 0.0 * angstrom, pl);
      if (!v.debye_temperature && !v.rms) {
        throw std::invalid_argument(vw + ": give debye_temperature (e.g. \"470 K\") or rms (e.g. \"0.06 A\")");
      }
      out.crystal.vibrations.push_back(v);
    }
  }
  if (const auto *pv = cr.get("polycrystal"); pv != nullptr && !pv->is_null()) {
    // A polycrystal: generated grains, or grains from a file (read here).
    const Reader pr(cr.object("polycrystal"), cr.where("polycrystal"));
    const std::string &pw = pr.path();
    pr.reject_unknown({"source", "grain_size", "grain_size_A", "columnar", "texture", "spread", "spread_deg", "seed",
                       "period", "boundary_width", "boundary_width_A", "file", "scale", "scale_A", "misorientation",
                       "misorientation_deg", "phase"});
    crystal::Polycrystal p;
    p.file = pr.str("file");
    try {
      if (pr.has("source")) p.source = crystal::parse_polycrystal_source(lower(pr.str("source")));
      else if (!p.file.empty()) p.source = crystal::grain_source_from_path(p.file);
    } catch (const std::invalid_argument &e) {
      throw std::invalid_argument(pr.where("source") + ": " + e.what());
    }
    p.grain_size = pr.quantity("grain_size", "grain_size_A", angstrom, p.grain_size, pl);
    p.columnar = pr.flag("columnar", p.columnar);
    if (pr.has("texture")) {
      try {
        p.texture = crystal::parse_texture(lower(pr.str("texture")));
      } catch (const std::invalid_argument &e) {
        throw std::invalid_argument(pr.where("texture") + ": " + e.what());
      }
    }
    p.spread = pr.quantity("spread", "spread_deg", degree, p.spread, pa);
    const double seed = pr.num("seed", static_cast<double>(p.seed));
    if (!(seed >= 0.0) || seed != std::floor(seed) || seed > 9.0e15) {
      throw std::invalid_argument(pr.where("seed") + ": a whole number from 0");
    }
    p.seed = static_cast<std::uint64_t>(seed);
    if (const auto *period = pr.get("period"); period != nullptr && !period->is_null()) {
      // One length for both lateral directions, or [y, z].
      const auto one = [&](const json::value &v, const std::string &w) {
        if (v.is_string()) return parse_length(std::string(v.as_string()), w);
        if (v.is_number()) return v.to_number<double>() * angstrom;
        throw std::invalid_argument(w + ": expected a length");
      };
      if (period->is_array()) {
        const auto &a = period->as_array();
        if (a.size() != 2) throw std::invalid_argument(pr.where("period") + ": one length, or [y, z]");
        p.period_y = one(a[0], pr.where("period") + ".0");
        p.period_z = one(a[1], pr.where("period") + ".1");
      } else {
        p.period_y = p.period_z = one(*period, pr.where("period"));
      }
    }
    p.boundary_width = pr.quantity("boundary_width", "boundary_width_A", angstrom, p.boundary_width, pl);
    p.scale = pr.quantity("scale", "scale_A", angstrom, p.scale, pl);
    p.misorientation = pr.quantity("misorientation", "misorientation_deg", degree, p.misorientation, pa);
    const double phase = pr.num("phase", 0.0);
    if (!(phase >= 0.0) || phase != std::floor(phase) || phase > 1.0e6) {
      throw std::invalid_argument(pr.where("phase") + ": 0 (every phase) or a phase number");
    }
    p.phase = static_cast<int>(phase);
    crystal::check_polycrystal(p, pw);
    if (p.source != crystal::Polycrystal::Source::voronoi) {
      const std::filesystem::path path = std::filesystem::path(p.file).is_absolute() ? std::filesystem::path(p.file)
                                                                                      : base_dir / p.file;
      try {
        p.map = crystal::read_grain_map(p, path, cell);
      } catch (const std::invalid_argument &e) {
        throw std::invalid_argument(pr.where("file") + ": " + e.what());
      }
    }
    out.crystal.polycrystal = std::move(p);
  }
  return out;
}

/*!
 * A layer's crystal: sizes a built-in structure from the density when it has
 * no lattice constant, then makes the layer's elements, stoichiometry and
 * density the unit cell's. What the scenario says about them must agree.
 */
void read_layer_crystal(trim::TargetLayer &layer, const Reader &cr, const std::filesystem::path &base_dir,
                        const std::string &lw, bool density_given, const std::vector<bool> &stoichiometry_given) {
  CrystalDescription d = parse_crystal_object(cr, base_dir);
  const std::string &cw = cr.path();
  if (layer.gas) throw std::invalid_argument(cw + ": a gas layer cannot be a crystal");
  const auto composition = d.crystal.cell.composition();
  if (composition.empty()) throw std::invalid_argument(cw + ": the unit cell has no atoms");
  std::map<unsigned, mass_quantity> masses;
  for (std::size_t k = 0; k < layer.elements.size(); ++k) {
    const auto &e = layer.elements[k];
    if (composition.count(e.z) == 0) {
      throw std::invalid_argument(lw + ".elements." + std::to_string(k) + ": " +
                                  std::string(Element::by_z(e.z).symbol()) + " is not in the crystal's unit cell");
    }
    masses[e.z] = e.mass;
  }
  const bool from_density = d.lattice_from_density;
  if (from_density) {
    mass_density_quantity rho = density_given ? layer.density : 0.0 * gram_per_cubic_centimeter;
    if (!density_given && composition.size() == 1) rho = Element::by_z(composition.begin()->first).solid_density();
    if (!(rho > 0.0 * gram_per_cubic_centimeter)) {
      throw std::invalid_argument(cw + ": give the lattice constant (\"a\") or the layer's density");
    }
    d.crystal.cell.scale_to_density(rho, masses);
  }
  // Stoichiometries written in the scenario must be the cell's.
  const bool all_given = !layer.elements.empty() && layer.elements.size() == composition.size() &&
                         std::all_of(stoichiometry_given.begin(), stoichiometry_given.end(), [](bool b) { return b; });
  if (all_given) {
    double given_sum = 0.0;
    for (const auto &e : layer.elements) given_sum += e.stoichiometry;
    const double cell_sum = d.crystal.cell.atoms_per_cell();
    for (std::size_t k = 0; k < layer.elements.size() && given_sum > 0.0; ++k) {
      const auto &e = layer.elements[k];
      const double given = e.stoichiometry / given_sum;
      const double in_cell = composition.at(e.z) / cell_sum;
      if (std::abs(given - in_cell) > 0.005) {
        throw std::invalid_argument(lw + ".elements." + std::to_string(k) + ".stoichiometry: the fraction of " +
                                    std::string(Element::by_z(e.z).symbol()) + " is " + number(given) +
                                    " here and " + number(in_cell) + " in the crystal's unit cell; remove the "
                                    "stoichiometry or correct the crystal");
      }
    }
  }
  const double rho_given = in_units(layer.density, gram_per_cubic_centimeter);
  layer.crystal = std::move(d.crystal);
  trim::apply_crystal(layer);
  if (!from_density && density_given) {
    const double rho_cell = in_units(layer.density, gram_per_cubic_centimeter);
    if (std::abs(rho_given - rho_cell) > 0.02 * rho_cell) {
      throw std::invalid_argument(lw + ".density: " + number(rho_given) + " g/cm3 disagrees with the crystal's unit "
                                  "cell (" + number(rho_cell) + " g/cm3 from its lattice and atoms); remove one of the "
                                  "two (a built-in structure without \"a\" takes its lattice constant from the density)");
    }
  }
}

json::object crystal_object(const crystal::Crystal &c) {
  json::object co;
  if (!c.cell.name.empty()) co["name"] = c.cell.name;
  json::object lat;
  json::array rows;
  for (int i = 0; i < 3; ++i) {
    json::array r;
    for (int j = 0; j < 3; ++j) r.emplace_back(c.cell.lattice.vectors_A()[i][j]);
    rows.emplace_back(std::move(r));
  }
  lat["vectors_A"] = std::move(rows);
  co["lattice"] = std::move(lat);
  json::array sites;
  for (const auto &s : c.cell.sites) {
    json::object so;
    so["label"] = s.label;
    json::array p;
    for (int j = 0; j < 3; ++j) p.emplace_back(s.position[j]);
    so["position"] = std::move(p);
    so["kind"] = crystal::to_string(s.kind);
    // One element: "element" (and "fraction" below 1); several: a map;
    // none: an empty list.
    if (s.occupants.size() == 1) {
      so["element"] = std::string(Element::by_z(s.occupants.front().z).symbol());
      if (s.occupants.front().fraction != 1.0) so["fraction"] = s.occupants.front().fraction;
    } else if (s.occupants.empty()) {
      so["occupants"] = json::array();
    } else {
      json::object occ;
      for (const auto &o : s.occupants) occ[std::string(Element::by_z(o.z).symbol())] = o.fraction;
      so["occupants"] = std::move(occ);
    }
    sites.emplace_back(std::move(so));
  }
  co["sites"] = std::move(sites);
  json::object oo;
  if (c.orientation.kind == crystal::OrientationKind::euler) {
    json::array e;
    for (double a : c.orientation.euler_deg) e.emplace_back(a);
    oo["euler_deg"] = std::move(e);
  } else {
    json::array surface;
    for (double h : c.orientation.surface) surface.emplace_back(index_value(h));
    oo["surface"] = std::move(surface);
    if (!c.orientation.reference.empty()) {
      json::array reference;
      for (double u : c.orientation.reference) reference.emplace_back(index_value(u));
      oo["reference"] = std::move(reference);
    }
  }
  oo["twist"] = format_angle(c.orientation.twist);
  co["orientation"] = std::move(oo);
  co["temperature"] = format_temperature(c.temperature);
  if (!c.vibrations.empty()) {
    json::array vs;
    for (const auto &v : c.vibrations) {
      json::object vo;
      vo["element"] = std::string(Element::by_z(v.z).symbol());
      if (v.debye_temperature) vo["debye_temperature"] = format_temperature(*v.debye_temperature);
      if (v.rms) vo["rms"] = format_length(*v.rms);
      vs.emplace_back(std::move(vo));
    }
    co["vibrations"] = std::move(vs);
  }
  if (c.polycrystal) {
    const auto &p = *c.polycrystal;
    json::object po;
    po["source"] = crystal::to_string(p.source);
    if (p.source == crystal::Polycrystal::Source::voronoi) {
      po["grain_size"] = format_length(p.grain_size);
      if (p.columnar) po["columnar"] = true;
      po["texture"] = crystal::to_string(p.texture);
      if (p.spread != 0.0 * degree) po["spread"] = format_angle(p.spread);
      po["seed"] = p.seed;
    } else {
      po["file"] = p.file;
      if (p.source == crystal::Polycrystal::Source::neper) po["scale"] = format_length(p.scale);
      if (p.source == crystal::Polycrystal::Source::ebsd) {
        po["misorientation"] = format_angle(p.misorientation);
        if (p.phase != 0) po["phase"] = p.phase;
      }
    }
    if (p.period_y > 0.0 * angstrom || p.period_z > 0.0 * angstrom) {
      json::array period;
      period.emplace_back(format_length(p.period_y));
      period.emplace_back(format_length(p.period_z));
      po["period"] = std::move(period);
    }
    if (p.boundary_width > 0.0 * angstrom) po["boundary_width"] = format_length(p.boundary_width);
    co["polycrystal"] = std::move(po);
  }
  return co;
}

physics::ScatteringMethod parse_scattering(const std::string &s, const std::string &where) {
  const auto t = lower(s);
  if (t == "magic") return physics::ScatteringMethod::magic;
  if (t == "integral" || t == "exact" || t == "gauss-mehler") return physics::ScatteringMethod::integral;
  throw std::invalid_argument(where + ": expected magic or integral");
}

trim::TrimConfig::PairPotentialUse parse_pair_use(const std::string &s, const std::string &where) {
  const auto t = lower(s);
  if (t == "none" || t == "universal" || t == "off") return trim::TrimConfig::PairPotentialUse::none;
  if (t == "recoils") return trim::TrimConfig::PairPotentialUse::recoils;
  if (t == "all") return trim::TrimConfig::PairPotentialUse::all;
  throw std::invalid_argument(where + ": expected none, recoils or all");
}

std::string to_string(trim::TrimConfig::PairPotentialUse u) {
  switch (u) {
  case trim::TrimConfig::PairPotentialUse::none: return "none";
  case trim::TrimConfig::PairPotentialUse::recoils: return "recoils";
  case trim::TrimConfig::PairPotentialUse::all: return "all";
  }
  return "recoils";
}

trim::TrimConfig::StragglingDistribution parse_straggling(const std::string &s,
                                                          const std::string &where) {
  const auto t = lower(s);
  if (t == "gamma") return trim::TrimConfig::StragglingDistribution::gamma;
  if (t == "truncated-gaussian" || t == "truncated_gaussian" || t == "gaussian") {
    return trim::TrimConfig::StragglingDistribution::truncated_gaussian;
  }
  throw std::invalid_argument(where + ": expected gamma or truncated-gaussian");
}

} // namespace

std::string to_string(Mode m) { return m == Mode::native ? "native" : "srim"; }

Mode parse_mode(std::string_view s) {
  const auto t = lower(s);
  if (t == "srim" || t == "trim" || t == "srim-mode") return Mode::srim;
  if (t == "native" || t == "amcset") return Mode::native;
  throw std::invalid_argument("mode: expected srim or native, got '" + std::string(s) + "'");
}

// ---------------------------------------------------------------------------
// parse_config
// ---------------------------------------------------------------------------

RunConfig parse_config(std::string_view json_text, const std::filesystem::path &base_dir) {
  json::value root;
  try {
    json::parse_options opt;
    opt.allow_comments = true;
    opt.allow_trailing_commas = true;
    root = json::parse(json_text, {}, opt);
  } catch (const std::exception &e) {
    throw std::invalid_argument(std::string("Invalid JSON: ") + e.what());
  }
  if (!root.is_object()) throw std::invalid_argument("The scenario must be a JSON object");
  const Reader r(root.as_object(), "");
  r.reject_unknown({"mode", "title", "ion", "layers", "target", "ions", "first_ion", "seed", "threads",
                    "damage", "source", "stopping", "srim_data", "physics", "output", "native",
                    "auto_depth", "$schema", "comment", "description"});

  RunConfig c;
  auto &t = c.trim;
  c.mode = parse_mode(r.str("mode", "srim"));
  t.title = r.str("title", t.title);

  // --- Ion ------------------------------------------------------------------
  {
    const Reader ion(r.object("ion"), "ion");
    if (r.object("ion").empty()) throw std::invalid_argument("ion: missing");
    ion.reject_unknown({"element", "z", "mass", "mass_amu", "energy", "energy_keV", "angle",
                        "angle_deg", "mass_number", "spot"});
    t.ion.z = element_z(ion);
    t.ion.mass = ion.quantity("mass", "mass_amu", atomic_mass_unit,
                              Element::by_z(t.ion.z).most_abundant_isotope_mass(),
                              [](std::string_view s, const std::string &w) { return parse_mass(s, w); });
    t.ion.energy = ion.quantity("energy", "energy_keV", kilo_electron_volt, 0.0 * electron_volt,
                                [](std::string_view s, const std::string &w) { return parse_energy(s, w); });
    if (!(t.ion.energy > 0.0 * electron_volt)) throw std::invalid_argument("ion.energy: missing or not positive");
    t.ion.angle = ion.quantity("angle", "angle_deg", degree, 0.0 * degree,
                               [](std::string_view s, const std::string &w) { return parse_angle(s, w); });
    c.native.ion_mass_number = static_cast<unsigned>(ion.count("mass_number", 0, 300));
    if (ion.has("spot")) {
      // Where the ions meet the surface (crystal transport; not TRIM).
      const Reader sp(ion.object("spot"), "ion.spot");
      sp.reject_unknown({"shape", "width", "height", "size", "diameter", "fwhm"});
      auto &spot = t.ion.spot;
      if (sp.has("shape")) {
        spot.shape = with_where(sp.where("shape"), [&] { return trim::parse_beam_spot_shape(sp.str("shape")); });
      }
      const auto pl = [](std::string_view s, const std::string &w) { return parse_length(s, w); };
      // "size" (also "diameter", "fwhm") sets both; "width" and "height" each.
      for (const char *both : {"size", "diameter", "fwhm"}) {
        if (!sp.has(both)) continue;
        spot.width = sp.quantity(both, "", angstrom, spot.width, pl);
        spot.height = spot.width;
      }
      spot.width = sp.quantity("width", "", angstrom, spot.width, pl);
      spot.height = sp.quantity("height", "", angstrom, spot.height, pl);
    }
  }

  // --- Target ---------------------------------------------------------------
  const json::array *layers = r.array("layers");
  if (!layers) layers = r.array("target");
  if (!layers || layers->empty()) throw std::invalid_argument("layers: at least one layer is required");
  for (std::size_t i = 0; i < layers->size(); ++i) {
    const std::string lw = "layers." + std::to_string(i);
    if (!(*layers)[i].is_object()) throw std::invalid_argument(lw + ": expected an object");
    const Reader lr((*layers)[i].as_object(), lw);
    lr.reject_unknown({"name", "thickness", "width_A", "density", "density_g_cm3", "gas",
                       "compound_correction", "elements", "crystal"});
    trim::TargetLayer layer;
    layer.width = lr.quantity("thickness", "width_A", angstrom, 0.0 * angstrom,
                              [](std::string_view s, const std::string &w) { return parse_length(s, w); });
    layer.gas = lr.flag("gas", false);
    layer.compound_correction = lr.num("compound_correction", 1.0);
    const json::array *els = lr.array("elements");
    // With a crystal the elements may be left out: they are the unit cell's.
    const bool has_crystal = lr.has("crystal") && !lr.get("crystal")->is_null();
    if ((!els || els->empty()) && !has_crystal) {
      throw std::invalid_argument(lw + ".elements: at least one element is required");
    }
    std::vector<bool> stoichiometry_given;
    for (std::size_t k = 0; els && k < els->size(); ++k) {
      const std::string ew = lw + ".elements." + std::to_string(k);
      if (!(*els)[k].is_object()) throw std::invalid_argument(ew + ": expected an object");
      const Reader er((*els)[k].as_object(), ew);
      er.reject_unknown({"element", "z", "name", "stoichiometry", "fraction", "mass", "mass_amu",
                         "Ed", "Ed_eV", "Eb", "Eb_eV", "Es", "Es_eV", "trim_atom", "damage_efficiency"});
      const unsigned z = element_z(er);
      auto el = trim::TargetElement::defaults(z, 1.0);
      el.stoichiometry = er.num("stoichiometry", er.num("fraction", 1.0));
      stoichiometry_given.push_back(er.has("stoichiometry") || er.has("fraction"));
      el.name = er.str("name", el.name);
      const auto pe = [](std::string_view s, const std::string &w) { return parse_energy(s, w); };
      el.mass = er.quantity("mass", "mass_amu", atomic_mass_unit, el.mass,
                            [](std::string_view s, const std::string &w) { return parse_mass(s, w); });
      el.displacement_energy = er.quantity("Ed", "Ed_eV", electron_volt, el.displacement_energy, pe);
      el.lattice_binding_energy = er.quantity("Eb", "Eb_eV", electron_volt, el.lattice_binding_energy, pe);
      el.surface_binding_energy = er.quantity("Es", "Es_eV", electron_volt, el.surface_binding_energy, pe);
      el.trim_atom = static_cast<unsigned>(er.count("trim_atom", 0, 100000));
      if (er.has("damage_efficiency")) {
        // arc-dpa / rpa constants (Nordlund et al. 2018); all four are needed.
        const Reader de(er.object("damage_efficiency"), ew + ".damage_efficiency");
        de.reject_unknown({"arc_b", "arc_c", "rpa_b", "rpa_c"});
        for (const char *key : {"arc_b", "arc_c", "rpa_b", "rpa_c"}) {
          if (!de.has(key)) throw std::invalid_argument(ew + ".damage_efficiency." + key + ": missing (give all four)");
        }
        trim::DamageEfficiency d;
        d.arc_b = de.num("arc_b", 0.0);
        d.arc_c = de.num("arc_c", 1.0);
        d.rpa_b = de.quantity("rpa_b", "", electron_volt, 0.0 * electron_volt, pe);
        d.rpa_c = de.num("rpa_c", 1.0);
        el.damage_efficiency = d;
      }
      layer.elements.push_back(el);
    }
    layer.name = lr.str("name", "");
    const bool density_given = lr.has("density") || lr.has("density_g_cm3");
    // Default density: the element's own for single-element layers.
    mass_density_quantity def_density = 0.0 * gram_per_cubic_centimeter;
    if (layer.elements.size() == 1) {
      const auto &el = Element::by_z(layer.elements.front().z);
      def_density = layer.gas ? el.gas_density() : el.solid_density();
    }
    layer.density = lr.quantity("density", "density_g_cm3", gram_per_cubic_centimeter, def_density,
                                [](std::string_view s, const std::string &w) { return parse_density(s, w); });
    if (has_crystal) {
      read_layer_crystal(layer, Reader(lr.object("crystal"), lw + ".crystal"), base_dir, lw, density_given,
                         stoichiometry_given);
    } else if (!(layer.density > 0.0 * gram_per_cubic_centimeter)) {
      throw std::invalid_argument(lw + ".density: required for compounds");
    }
    if (layer.name.empty()) {
      for (const auto &e : layer.elements) layer.name += std::string(Element::by_z(e.z).symbol());
    }
    t.layers.push_back(std::move(layer));
  }

  // --- Calculation ----------------------------------------------------------
  t.ion_count = r.count("ions", t.ion_count);
  t.first_ion = r.count("first_ion", t.first_ion);
  t.seed = static_cast<std::uint32_t>(r.count("seed", t.seed, 4294967295u));
  t.threads = static_cast<unsigned>(r.count("threads", 0, trim::TrimConfig::max_threads));
  t.damage = trim::parse_damage_calculation(r.str("damage", trim::to_string(t.damage)));
  {
    // Automatic target depth (trim::size_target) and the escape warning.
    const Reader a(r.object("auto_depth"), "auto_depth");
    a.reject_unknown({"enabled", "test_ions", "margin", "warn_transmitted"});
    auto &ad = t.auto_depth;
    ad.enabled = a.flag("enabled", ad.enabled);
    ad.test_ions = a.count("test_ions", ad.test_ions, 10000000);
    if (ad.test_ions == 0) throw std::invalid_argument("auto_depth.test_ions: must be at least 1");
    ad.margin = a.num("margin", ad.margin);
    if (!(ad.margin >= 0.01 && ad.margin <= 0.5)) {
      throw std::invalid_argument("auto_depth.margin: the fraction at the back of the target in which the "
                                  "deepest ion must stop, 0.01 to 0.5 (default 0.1)");
    }
    ad.warn_transmitted = a.flag("warn_transmitted", ad.warn_transmitted);
  }

  // --- Particle source (TRIM.DAT) ----------------------------------------------
  if (r.has("source")) {
    const Reader sr(r.object("source"), "source");
    sr.reject_unknown({"type", "file", "particles", "title"});
    const std::string type = lower(sr.str("type", "beam"));
    if (type == "beam") {
      t.source = trim::TrimConfig::Source::beam;
    } else if (type == "ions" || type == "trim.dat ions" || type == "particle_ions") {
      t.source = trim::TrimConfig::Source::particle_ions;
    } else if (type == "recoils" || type == "trim.dat recoils" || type == "particle_recoils" ||
               type == "neutrons") {
      t.source = trim::TrimConfig::Source::particle_recoils;
    } else {
      throw std::invalid_argument(sr.where("type") + ": expected beam, ions or recoils");
    }
    if (t.source != trim::TrimConfig::Source::beam) {
      const std::string file = resolve(sr.str("file", ""), base_dir);
      if (!file.empty()) {
        auto list = trim::read_trim_dat(file);
        t.particles = std::move(list.particles);
        t.particles_title = list.title;
      }
      if (const json::array *ps = sr.array("particles")) {
        for (std::size_t i = 0; i < ps->size(); ++i) {
          const std::string pw = "source.particles." + std::to_string(i);
          if (!(*ps)[i].is_object()) throw std::invalid_argument(pw + ": expected an object");
          const Reader pr((*ps)[i].as_object(), pw);
          pr.reject_unknown({"name", "element", "z", "energy", "x", "y", "z_position", "direction"});
          trim::SourceParticle sp;
          sp.name = pr.str("name", "P" + std::to_string(i + 1));
          sp.z = element_z(pr);
          const auto pe = [](std::string_view s, const std::string &w) { return parse_energy(s, w); };
          const auto pl = [](std::string_view s, const std::string &w) { return parse_length(s, w); };
          sp.energy = pr.quantity("energy", "", electron_volt, 0.0 * electron_volt, pe);
          sp.x = pr.quantity("x", "", angstrom, 0.0 * angstrom, pl);
          sp.y = pr.quantity("y", "", angstrom, 0.0 * angstrom, pl);
          sp.z_position = pr.quantity("z_position", "", angstrom, 0.0 * angstrom, pl);
          if (const json::array *dir = pr.array("direction")) {
            if (dir->size() != 3) throw std::invalid_argument(pw + ".direction: expected [cos x, cos y, cos z]");
            double v[3];
            for (int k = 0; k < 3; ++k) {
              const auto &x = (*dir)[static_cast<std::size_t>(k)];
              if (x.is_double()) v[k] = x.as_double();
              else if (x.is_int64()) v[k] = static_cast<double>(x.as_int64());
              else throw std::invalid_argument(pw + ".direction: expected numbers");
            }
            sp.cos_x = v[0];
            sp.cos_y = v[1];
            sp.cos_z = v[2];
          }
          t.particles.push_back(std::move(sp));
        }
      }
      t.particles_title = sr.str("title", t.particles_title);
      if (t.particles.empty()) {
        throw std::invalid_argument("source: a particle source needs 'file' (TRIM.DAT) or 'particles'");
      }
      if (!r.has("ions")) t.ion_count = 0;
    }
  }

  {
    const auto &so = r.object("stopping");
    const Reader sr(so, "stopping");
    sr.reject_unknown({"model", "table"});
    t.stopping_model = physics::parse_electronic_stopping_model(
        sr.str("model", physics::to_string(t.stopping_model)));
    t.stopping_table = resolve(sr.str("table", ""), base_dir);
  }
  t.srim_data_directory = resolve(r.str("srim_data", ""), base_dir);

  {
    const Reader p(r.object("physics"), "physics");
    p.reject_unknown({"transport", "pair_potentials", "scattering", "electronic_straggling",
                      "straggling_distribution", "free_flight_path", "min_energy_transfer",
                      "final_energy", "ion_stop_energy", "recoil_stop_energy",
                      "max_impact_parameter_ws", "kp_subtract_binding_energy",
                      "surface_atom_depth", "max_recoils_per_ion", "bethe_bloch", "trim_bugs",
                      "fast_math", "crystal"});
    const auto pe = [](std::string_view s, const std::string &w) { return parse_energy(s, w); };
    if (p.has("transport")) {
      t.transport = with_where(p.where("transport"), [&] { return trim::parse_transport(p.str("transport")); });
    }
    if (p.has("crystal")) {
      // Settings of the crystal transport (not TRIM).
      const Reader cp(p.object("crystal"), "physics.crystal");
      cp.reject_unknown({"max_impact_parameter", "simultaneous_distance", "local_stopping_fraction",
                         "local_stopping_decay", "thermal_vibrations"});
      const auto pl = [](std::string_view s, const std::string &w) { return parse_length(s, w); };
      auto &cr = t.crystal;
      if (const auto *v = cp.get("max_impact_parameter"); v && v->is_string() && lower(v->as_string()) == "auto") {
        cr.max_impact_parameter = 0.0 * angstrom;
      } else {
        cr.max_impact_parameter = cp.quantity("max_impact_parameter", "", angstrom, cr.max_impact_parameter, pl);
      }
      cr.simultaneous_distance = cp.quantity("simultaneous_distance", "", angstrom, cr.simultaneous_distance, pl);
      cr.local_stopping_fraction = cp.num("local_stopping_fraction", cr.local_stopping_fraction);
      cr.local_stopping_decay = cp.num("local_stopping_decay", cr.local_stopping_decay);
      cr.thermal_vibrations = cp.flag("thermal_vibrations", cr.thermal_vibrations);
    }
    if (p.has("pair_potentials")) t.pair_potentials = parse_pair_use(p.str("pair_potentials"), p.where("pair_potentials"));
    if (p.has("scattering")) t.scattering = parse_scattering(p.str("scattering"), p.where("scattering"));
    t.electronic_straggling = p.flag("electronic_straggling", t.electronic_straggling);
    if (p.has("straggling_distribution")) {
      t.straggling_distribution = parse_straggling(p.str("straggling_distribution"), p.where("straggling_distribution"));
    }
    t.free_flight_path = p.flag("free_flight_path", t.free_flight_path);
    t.min_energy_transfer = p.quantity("min_energy_transfer", "", electron_volt, t.min_energy_transfer, pe);
    t.final_energy = p.quantity("final_energy", "", electron_volt, t.final_energy, pe);
    t.ion_stop_energy = p.quantity("ion_stop_energy", "", electron_volt, t.ion_stop_energy, pe);
    t.recoil_stop_energy = p.quantity("recoil_stop_energy", "", electron_volt, t.recoil_stop_energy, pe);
    t.max_impact_parameter_ws = p.num("max_impact_parameter_ws", t.max_impact_parameter_ws);
    t.kp_subtract_binding_energy = p.flag("kp_subtract_binding_energy", t.kp_subtract_binding_energy);
    t.surface_atom_depth = p.quantity("surface_atom_depth", "", angstrom, t.surface_atom_depth,
                                      [](std::string_view s, const std::string &w) { return parse_length(s, w); });
    t.max_recoils_per_ion = p.count("max_recoils_per_ion", t.max_recoils_per_ion);
    t.fast_math = p.flag("fast_math", t.fast_math);
    if (const auto *tb = p.get("trim_bugs"); tb && !tb->is_null()) {
      // "trim" (reproduce every bug of TRIM.exe, the default), "fixed", or an
      // object of flags (true = as TRIM).
      if (tb->is_string()) {
        const auto v = lower(tb->as_string());
        if (v == "trim" || v == "reproduce") t.bugs.set_all(true);
        else if (v == "fixed" || v == "fix") t.bugs.set_all(false);
        else throw std::invalid_argument("physics.trim_bugs: expected \"trim\", \"fixed\" or an object");
      } else {
        const Reader b(p.object("trim_bugs"), "physics.trim_bugs");
        std::vector<std::string> names;
        for (const auto &m : trim::TrimBugs::members()) names.emplace_back(m.first);
        b.reject_unknown(names);
        for (const auto &m : trim::TrimBugs::members()) t.bugs.*(m.second) = b.flag(m.first, t.bugs.*(m.second));
      }
    }
    const Reader bb(p.object("bethe_bloch"), "physics.bethe_bloch");
    bb.reject_unknown({"effective_charge", "density_effect"});
    t.bethe_bloch.effective_charge = bb.flag("effective_charge", t.bethe_bloch.effective_charge);
    t.bethe_bloch.density_effect = bb.flag("density_effect", t.bethe_bloch.density_effect);
  }

  {
    const Reader o(r.object("output"), "output");
    o.reject_unknown({"directory", "depth_window", "bins", "range_3d", "exit_lists",
                      "max_exit_records", "collision_details", "trajectories", "exyz_interval",
                      "files", "trim_single_precision", "collision_recoils", "trajectory_file",
                      "md_events", "dpa_fluence"});
    t.output_directory = resolve(o.str("directory", ""), base_dir);
    if (const auto *w = o.array("depth_window")) {
      if (w->size() != 2 || !(*w)[0].is_string() || !(*w)[1].is_string()) {
        throw std::invalid_argument("output.depth_window: expected [\"start\", \"end\"] with units");
      }
      t.window_start = parse_length(std::string_view((*w)[0].as_string()), "output.depth_window[0]");
      t.window_end = parse_length(std::string_view((*w)[1].as_string()), "output.depth_window[1]");
    }
    t.depth_bins = o.count("bins", t.depth_bins, trim::TrimConfig::max_depth_bins);
    t.record_ion_positions = o.flag("range_3d", t.record_ion_positions);
    t.record_exiting_atoms = o.flag("exit_lists", t.record_exiting_atoms);
    t.max_exit_records = o.count("max_exit_records", t.max_exit_records);
    t.collision_details_ions = o.count("collision_details", t.collision_details_ions);
    t.collision_recoil_details = o.flag("collision_recoils", t.collision_recoil_details);
    if (const auto *tr = o.get("trajectories"); tr != nullptr && tr->is_string()) {
      if (lower(std::string(tr->as_string())) != "all") {
        throw std::invalid_argument("output.trajectories: a number of ions or \"all\", got '" +
                                    std::string(tr->as_string()) + "'");
      }
      t.trajectory_ions = trim::TrimConfig::all_ions;
    } else {
      t.trajectory_ions = o.count("trajectories", t.trajectory_ions);
    }
    if (const auto *f = o.get("dpa_fluence"); f != nullptr && !f->is_null()) {
      // Ions per cm2: a number, or a string with "/cm2", "ions/cm2" or "/m2".
      if (f->is_number()) {
        t.dpa_fluence_per_cm2 = o.num("dpa_fluence", t.dpa_fluence_per_cm2);
      } else if (f->is_string()) {
        const std::string text(f->as_string());
        char *end = nullptr;
        const double v = std::strtod(text.c_str(), &end);
        std::string unit(end);
        unit.erase(std::remove(unit.begin(), unit.end(), ' '), unit.end());
        if (end == text.c_str() || !std::isfinite(v)) {
          throw std::invalid_argument("output.dpa_fluence: '" + text + "' does not start with a number");
        }
        if (unit.empty() || unit == "/cm2" || unit == "ions/cm2" || unit == "cm-2") {
          t.dpa_fluence_per_cm2 = v;
        } else if (unit == "/m2" || unit == "ions/m2" || unit == "m-2") {
          t.dpa_fluence_per_cm2 = v * 1e-4;
        } else {
          throw std::invalid_argument("output.dpa_fluence: unknown unit '" + unit + "' (ions/cm2 or ions/m2)");
        }
      } else {
        throw std::invalid_argument("output.dpa_fluence: ions per cm2, e.g. 1e15 or \"1e15 /cm2\"");
      }
    }
    const auto pe = [](std::string_view s, const std::string &w) { return parse_energy(s, w); };
    if (o.has("trajectory_file")) {
      // Export of the recorded trajectories (docs/OUTPUT.md "Trajectory file").
      const Reader tf(o.object("trajectory_file"), "output.trajectory_file");
      tf.reject_unknown({"file", "points", "min_energy"});
      auto &te = t.trajectory_export;
      te.file = tf.str("file", "");
      const std::string points = lower(tf.str("points", "all"));
      if (points == "all") {
        te.births_only = false;
      } else if (points == "births") {
        te.births_only = true;
      } else {
        throw std::invalid_argument("output.trajectory_file.points: \"all\" or \"births\", got '" + points + "'");
      }
      te.min_energy = tf.quantity("min_energy", "", electron_volt, te.min_energy, pe);
    }
    if (o.has("md_events")) {
      // Hybrid BCA-MD recoil events (docs/OUTPUT.md "MD event file").
      const Reader m(o.object("md_events"), "output.md_events");
      m.reject_unknown({"file", "E_c", "E_min", "mode"});
      auto &md = t.md_events;
      md.file = m.str("file", "");
      if (m.has("E_c")) {
        md.transition_energies.clear();
        if (const auto *v = m.get("E_c"); v != nullptr && v->is_array()) {
          const auto &a = v->as_array();
          for (std::size_t i = 0; i < a.size(); ++i) {
            const std::string w = "output.md_events.E_c[" + std::to_string(i) + "]";
            if (!a[i].is_string()) throw std::invalid_argument(w + ": expected a string with a unit, e.g. \"1 keV\"");
            md.transition_energies.push_back(parse_energy(std::string_view(a[i].as_string()), w));
          }
        } else {
          md.transition_energies.push_back(m.quantity("E_c", "", electron_volt, 0.0 * electron_volt, pe));
        }
      }
      if (const auto *v = m.get("E_min"); v != nullptr && v->is_string() && lower(std::string(v->as_string())) == "ed") {
        md.min_energy = 0.0 * electron_volt;
      } else {
        md.min_energy = m.quantity("E_min", "", electron_volt, md.min_energy, pe);
      }
      const std::string mode = lower(m.str("mode", md.handoff ? "handoff" : "record"));
      if (mode == "record") {
        md.handoff = false;
      } else if (mode == "handoff") {
        md.handoff = true;
      } else {
        throw std::invalid_argument("output.md_events.mode: \"record\" or \"handoff\", got '" + mode + "'");
      }
    }
    t.single_precision_tables = o.flag("trim_single_precision", t.single_precision_tables);
    t.exyz_interval = o.quantity("exyz_interval", "", electron_volt, t.exyz_interval,
                                 [](std::string_view s, const std::string &w) { return parse_energy(s, w); });
    if (o.has("files")) {
      try {
        t.outputs = output::parse_output_spec(json::serialize(json::value(o.object("files"))));
      } catch (const std::invalid_argument &e) {
        throw std::invalid_argument(std::string("output.files: ") + e.what());
      }
    }
  }

  {
    const Reader n(r.object("native"), "native");
    n.reject_unknown({"electron", "damage_cascade", "screening_model", "electronic_stopping_model",
                      "custom_electronic_coefficient", "custom_electronic_exponent",
                      "custom_recoil_electronic_coefficient", "custom_recoil_electronic_exponent",
                      "custom_nuclear_screening_scale", "ion_stop_energy", "displacement_energy",
                      "electron_stop_energy", "angle_divisors", "flying_distances", "range", "bins",
                      "per_bombardment_output", "electron_screening_table", "mott_table"});
    auto &nv = c.native;
    const auto pe = [](std::string_view s, const std::string &w) { return parse_energy(s, w); };
    nv.electron = n.flag("electron", nv.electron);
    nv.damage_cascade = n.flag("damage_cascade", nv.damage_cascade);
    nv.screening_model = lower(n.str("screening_model", nv.screening_model));
    nv.electronic_stopping_model = lower(n.str("electronic_stopping_model", nv.electronic_stopping_model));
    nv.custom_electronic_coefficient = n.num("custom_electronic_coefficient", nv.custom_electronic_coefficient);
    nv.custom_electronic_exponent = n.num("custom_electronic_exponent", nv.custom_electronic_exponent);
    nv.custom_recoil_electronic_coefficient = n.num("custom_recoil_electronic_coefficient", nv.custom_recoil_electronic_coefficient);
    nv.custom_recoil_electronic_exponent = n.num("custom_recoil_electronic_exponent", nv.custom_recoil_electronic_exponent);
    nv.custom_nuclear_screening_scale = n.num("custom_nuclear_screening_scale", nv.custom_nuclear_screening_scale);
    nv.ion_stop_energy = n.quantity("ion_stop_energy", "", electron_volt, nv.ion_stop_energy, pe);
    nv.displacement_energy = n.quantity("displacement_energy", "", electron_volt, nv.displacement_energy, pe);
    nv.electron_stop_energy = n.quantity("electron_stop_energy", "", electron_volt, nv.electron_stop_energy, pe);
    // (Sizes of the native engine's tables: bounded so that a typo cannot
    // ask for terabytes.)
    nv.angle_divisors = n.count("angle_divisors", nv.angle_divisors, 10000000);
    nv.flying_distances = n.count("flying_distances", nv.flying_distances, 10000000);
    nv.range = n.quantity("range", "", angstrom, nv.range,
                          [](std::string_view s, const std::string &w) { return parse_length(s, w); });
    nv.bins = n.count("bins", nv.bins, 10000000);
    nv.per_bombardment_output = n.flag("per_bombardment_output", nv.per_bombardment_output);
    nv.electron_screening_table = n.str("electron_screening_table", nv.electron_screening_table);
    nv.mott_table = n.str("mott_table", nv.mott_table);
  }
  return c;
}

RunConfig load_config(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("Cannot open scenario file " + path.string());
  std::stringstream ss;
  ss << in.rdbuf();
  try {
    return parse_config(ss.str(), path.parent_path());
  } catch (const std::invalid_argument &e) {
    throw std::invalid_argument(path.string() + ": " + e.what());
  }
}

// ---------------------------------------------------------------------------
// to_json
// ---------------------------------------------------------------------------

namespace {

json::object to_object(const RunConfig &c) {
  const auto &t = c.trim;
  json::object root;
  root["mode"] = to_string(c.mode);
  root["title"] = t.title;
  json::object ion;
  ion["element"] = std::string(Element::by_z(t.ion.z).symbol());
  ion["mass"] = format_mass(t.ion.mass);
  ion["energy"] = format_energy(t.ion.energy);
  ion["angle"] = format_angle(t.ion.angle);
  if (!(t.ion.spot == trim::BeamSpot{})) {
    // (Written only when set: the default spot changes nothing in a file.)
    json::object spot;
    spot["shape"] = trim::to_string(t.ion.spot.shape);
    spot["width"] = format_length(t.ion.spot.width);
    spot["height"] = format_length(t.ion.spot.height);
    ion["spot"] = std::move(spot);
  }
  if (c.native.ion_mass_number != 0) ion["mass_number"] = c.native.ion_mass_number;
  root["ion"] = std::move(ion);
  json::array layers;
  for (const auto &l : t.layers) {
    json::object lo;
    lo["name"] = l.name;
    lo["thickness"] = format_length(l.width);
    lo["density"] = format_density(l.density);
    lo["gas"] = l.gas;
    lo["compound_correction"] = l.compound_correction;
    json::array els;
    for (const auto &e : l.elements) {
      json::object eo;
      eo["element"] = std::string(Element::by_z(e.z).symbol());
      if (!e.name.empty() && e.name != Element::by_z(e.z).symbol()) eo["name"] = e.name;
      eo["stoichiometry"] = e.stoichiometry;
      eo["mass"] = format_mass(e.mass);
      eo["Ed"] = format_energy(e.displacement_energy);
      eo["Eb"] = format_energy(e.lattice_binding_energy);
      eo["Es"] = format_energy(e.surface_binding_energy);
      if (e.trim_atom > 0) eo["trim_atom"] = e.trim_atom;
      if (e.damage_efficiency) {
        json::object de;
        de["arc_b"] = e.damage_efficiency->arc_b;
        de["arc_c"] = e.damage_efficiency->arc_c;
        de["rpa_b"] = format_energy(e.damage_efficiency->rpa_b);
        de["rpa_c"] = e.damage_efficiency->rpa_c;
        eo["damage_efficiency"] = std::move(de);
      }
      els.push_back(std::move(eo));
    }
    lo["elements"] = std::move(els);
    if (l.crystal) lo["crystal"] = crystal_object(*l.crystal);
    layers.push_back(std::move(lo));
  }
  root["layers"] = std::move(layers);
  root["ions"] = t.ion_count;
  if (t.first_ion != 0) root["first_ion"] = t.first_ion;
  {
    const trim::TrimConfig::AutoDepth d;
    const auto &ad = t.auto_depth;
    if (ad.enabled != d.enabled || ad.test_ions != d.test_ions || ad.margin != d.margin ||
        ad.warn_transmitted != d.warn_transmitted) {
      json::object a;
      a["enabled"] = ad.enabled;
      a["test_ions"] = ad.test_ions;
      a["margin"] = ad.margin;
      a["warn_transmitted"] = ad.warn_transmitted;
      root["auto_depth"] = std::move(a);
    }
  }
  root["seed"] = t.seed;
  root["threads"] = t.threads;
  root["damage"] = trim::to_string(t.damage);
  if (t.source != trim::TrimConfig::Source::beam) {
    json::object so;
    so["type"] = t.source == trim::TrimConfig::Source::particle_ions ? "ions" : "recoils";
    if (!t.particles_title.empty()) so["title"] = t.particles_title;
    json::array ps;
    for (const auto &sp : t.particles) {
      json::object po;
      po["name"] = sp.name;
      po["element"] = std::string(Element::by_z(sp.z).symbol());
      po["energy"] = format_energy(sp.energy);
      po["x"] = format_length(sp.x);
      po["y"] = format_length(sp.y);
      po["z_position"] = format_length(sp.z_position);
      // (Arrays are filled element by element: GCC's LTO warns, wrongly,
      // inside Boost.JSON's initializer-list construction.)
      json::array direction;
      direction.emplace_back(sp.cos_x);
      direction.emplace_back(sp.cos_y);
      direction.emplace_back(sp.cos_z);
      po["direction"] = std::move(direction);
      ps.push_back(std::move(po));
    }
    so["particles"] = std::move(ps);
    root["source"] = std::move(so);
  }
  json::object stopping;
  stopping["model"] = physics::to_string(t.stopping_model);
  stopping["table"] = t.stopping_table;
  root["stopping"] = std::move(stopping);
  root["srim_data"] = t.srim_data_directory;
  json::object p;
  p["transport"] = trim::to_string(t.transport);
  if (t.transport == trim::TrimConfig::Transport::crystal || !(t.crystal == trim::TrimConfig::CrystalPhysics{})) {
    json::object cr;
    if (t.crystal.max_impact_parameter > 0.0 * angstrom) {
      cr["max_impact_parameter"] = format_length(t.crystal.max_impact_parameter);
    } else {
      cr["max_impact_parameter"] = "auto";
    }
    cr["simultaneous_distance"] = format_length(t.crystal.simultaneous_distance);
    cr["local_stopping_fraction"] = t.crystal.local_stopping_fraction;
    cr["local_stopping_decay"] = t.crystal.local_stopping_decay;
    cr["thermal_vibrations"] = t.crystal.thermal_vibrations;
    p["crystal"] = std::move(cr);
  }
  p["pair_potentials"] = to_string(t.pair_potentials);
  p["scattering"] = t.scattering == physics::ScatteringMethod::magic ? "magic" : "integral";
  p["electronic_straggling"] = t.electronic_straggling;
  p["straggling_distribution"] =
      t.straggling_distribution == trim::TrimConfig::StragglingDistribution::gamma ? "gamma"
                                                                                   : "truncated-gaussian";
  p["free_flight_path"] = t.free_flight_path;
  p["min_energy_transfer"] = format_energy(t.min_energy_transfer);
  p["final_energy"] = format_energy(t.final_energy);
  p["ion_stop_energy"] = format_energy(t.ion_stop_energy);
  p["recoil_stop_energy"] = format_energy(t.recoil_stop_energy);
  p["max_impact_parameter_ws"] = t.max_impact_parameter_ws;
  p["kp_subtract_binding_energy"] = t.kp_subtract_binding_energy;
  p["surface_atom_depth"] = format_length(t.surface_atom_depth);
  p["max_recoils_per_ion"] = t.max_recoils_per_ion;
  p["fast_math"] = t.fast_math;
  json::object bb;
  bb["effective_charge"] = t.bethe_bloch.effective_charge;
  bb["density_effect"] = t.bethe_bloch.density_effect;
  p["bethe_bloch"] = std::move(bb);
  if (t.bugs.all_reproduced()) {
    p["trim_bugs"] = "trim";
  } else {
    json::object tb;
    for (const auto &m : trim::TrimBugs::members()) tb[m.first] = t.bugs.*(m.second);
    p["trim_bugs"] = std::move(tb);
  }
  root["physics"] = std::move(p);
  json::object o;
  o["directory"] = t.output_directory;
  {
    json::array window;
    window.emplace_back(format_length(t.window_start));
    window.emplace_back(format_length(t.window_end));
    o["depth_window"] = std::move(window);
  }
  o["bins"] = t.depth_bins;
  o["range_3d"] = t.record_ion_positions;
  o["exit_lists"] = t.record_exiting_atoms;
  o["max_exit_records"] = t.max_exit_records;
  o["collision_details"] = t.collision_details_ions;
  o["collision_recoils"] = t.collision_recoil_details;
  if (t.dpa_fluence_per_cm2 != 1e15) o["dpa_fluence"] = t.dpa_fluence_per_cm2;
  if (t.trajectory_ions == trim::TrimConfig::all_ions) {
    o["trajectories"] = "all";
  } else {
    o["trajectories"] = t.trajectory_ions;
  }
  if (!t.trajectory_export.file.empty()) {
    json::object tf;
    tf["file"] = t.trajectory_export.file;
    tf["points"] = t.trajectory_export.births_only ? "births" : "all";
    tf["min_energy"] = format_energy(t.trajectory_export.min_energy);
    o["trajectory_file"] = std::move(tf);
  }
  if (t.md_events.enabled()) {
    json::object m;
    m["file"] = t.md_events.file;
    json::array ec;
    for (const auto &e : t.md_events.transition_energies) ec.emplace_back(format_energy(e));
    m["E_c"] = std::move(ec);
    m["E_min"] = t.md_events.min_energy > 0.0 * electron_volt ? format_energy(t.md_events.min_energy) : std::string("Ed");
    m["mode"] = t.md_events.handoff ? "handoff" : "record";
    o["md_events"] = std::move(m);
  }
  if (t.exyz_interval > 0.0 * electron_volt) o["exyz_interval"] = format_energy(t.exyz_interval);
  o["trim_single_precision"] = t.single_precision_tables;
  o["files"] = json::parse(output::to_json(t.outputs));
  root["output"] = std::move(o);
  const auto &nv = c.native;
  json::object n;
  n["electron"] = nv.electron;
  n["damage_cascade"] = nv.damage_cascade;
  n["screening_model"] = nv.screening_model;
  n["electronic_stopping_model"] = nv.electronic_stopping_model;
  n["custom_electronic_coefficient"] = nv.custom_electronic_coefficient;
  n["custom_electronic_exponent"] = nv.custom_electronic_exponent;
  n["custom_recoil_electronic_coefficient"] = nv.custom_recoil_electronic_coefficient;
  n["custom_recoil_electronic_exponent"] = nv.custom_recoil_electronic_exponent;
  n["custom_nuclear_screening_scale"] = nv.custom_nuclear_screening_scale;
  n["ion_stop_energy"] = format_energy(nv.ion_stop_energy);
  n["displacement_energy"] = format_energy(nv.displacement_energy);
  n["electron_stop_energy"] = format_energy(nv.electron_stop_energy);
  n["angle_divisors"] = nv.angle_divisors;
  n["flying_distances"] = nv.flying_distances;
  n["range"] = format_length(nv.range);
  n["bins"] = nv.bins;
  n["per_bombardment_output"] = nv.per_bombardment_output;
  n["electron_screening_table"] = nv.electron_screening_table;
  n["mott_table"] = nv.mott_table;
  root["native"] = std::move(n);
  return root;
}

void pretty(std::ostream &os, const json::value &jv, int indent = 0) {
  const std::string pad(static_cast<std::size_t>(indent), ' ');
  switch (jv.kind()) {
  case json::kind::object: {
    const auto &obj = jv.get_object();
    if (obj.empty()) {
      os << "{}";
      return;
    }
    os << "{\n";
    std::size_t i = 0;
    for (const auto &kv : obj) {
      os << pad << "  " << json::serialize(json::value(kv.key())) << ": ";
      pretty(os, kv.value(), indent + 2);
      if (++i < obj.size()) os << ",";
      os << "\n";
    }
    os << pad << "}";
    return;
  }
  case json::kind::array: {
    const auto &arr = jv.get_array();
    bool simple = true;
    for (const auto &v : arr) simple = simple && !v.is_object() && !v.is_array();
    if (simple) {
      // On one line; numbers in the shortest form (0.5, not 5E-1).
      os << "[";
      for (std::size_t i = 0; i < arr.size(); ++i) {
        if (i) os << ",";
        pretty(os, arr[i], indent);
      }
      os << "]";
      return;
    }
    os << "[\n";
    for (std::size_t i = 0; i < arr.size(); ++i) {
      os << pad << "  ";
      pretty(os, arr[i], indent + 2);
      if (i + 1 < arr.size()) os << ",";
      os << "\n";
    }
    os << pad << "]";
    return;
  }
  case json::kind::double_: {
    // Shortest representation that reads back exactly ("1", "0.91",
    // "1e-05"), instead of Boost.JSON's "1E0".
    const double d = jv.get_double();
    if (!std::isfinite(d)) {
      os << "null";
      return;
    }
    char buf[32];
    for (int precision = 1; precision <= 17; ++precision) {
      std::snprintf(buf, sizeof buf, "%.*g", precision, d);
      if (std::strtod(buf, nullptr) == d) break;
    }
    std::string text(buf);
    // Keep it a JSON number that is visibly floating point.
    if (text.find_first_of(".eEn") == std::string::npos) text += ".0";
    os << text;
    return;
  }
  default:
    os << json::serialize(jv);
  }
}

} // namespace

std::string to_json(const RunConfig &config) {
  std::ostringstream os;
  pretty(os, json::value(to_object(config)));
  os << "\n";
  return os.str();
}

CrystalDescription parse_crystal(std::string_view json_text, const std::filesystem::path &base_dir,
                                 const std::string &where) {
  json::value root;
  try {
    json::parse_options opt;
    opt.allow_comments = true;
    opt.allow_trailing_commas = true;
    root = json::parse(json_text, {}, opt);
  } catch (const std::exception &e) {
    throw std::invalid_argument(where + ": invalid JSON: " + e.what());
  }
  if (!root.is_object()) throw std::invalid_argument(where + ": expected an object");
  return parse_crystal_object(Reader(root.as_object(), where), base_dir);
}

std::string crystal_to_json(const crystal::Crystal &crystal) {
  std::ostringstream os;
  pretty(os, json::value(crystal_object(crystal)));
  os << "\n";
  return os.str();
}

void save_config(const RunConfig &config, const std::filesystem::path &path) {
  if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
  std::ofstream out(path, std::ios::binary);
  if (!out) throw std::runtime_error("Cannot write " + path.string());
  out << to_json(config);
}

// ---------------------------------------------------------------------------
// Overrides
// ---------------------------------------------------------------------------

void apply_override(RunConfig &config, std::string_view assignment) {
  const auto eq = assignment.find('=');
  if (eq == std::string_view::npos) {
    throw std::invalid_argument("Override '" + std::string(assignment) + "' must look like key=value");
  }
  const std::string key = trim(assignment.substr(0, eq));
  const std::string value = trim(assignment.substr(eq + 1));
  json::value root(to_object(config));
  // Walk the dotted path (numbers index arrays).
  json::value *node = &root;
  std::string last;
  std::stringstream ks(key);
  std::string part;
  std::vector<std::string> parts;
  while (std::getline(ks, part, '.')) parts.push_back(part);
  if (parts.empty()) throw std::invalid_argument("Empty override key");
  for (std::size_t i = 0; i + 1 < parts.size(); ++i) {
    const auto &p = parts[i];
    if (node->is_object()) {
      auto &obj = node->as_object();
      if (!obj.contains(p)) throw std::invalid_argument("Unknown key '" + key + "'");
      node = &obj[p];
    } else if (node->is_array()) {
      std::size_t idx = 0;
      try {
        idx = std::stoul(p);
      } catch (const std::exception &) {
        throw std::invalid_argument("'" + p + "' in '" + key + "' must be an array index");
      }
      auto &arr = node->as_array();
      if (idx >= arr.size()) throw std::invalid_argument("Index " + p + " out of range in '" + key + "'");
      node = &arr[idx];
    } else {
      throw std::invalid_argument("Unknown key '" + key + "'");
    }
  }
  if (!node->is_object()) throw std::invalid_argument("Unknown key '" + key + "'");
  auto &obj = node->as_object();
  const auto &leaf = parts.back();
  json::value parsed;
  // Numbers and booleans as such, everything else as a string. A list is
  // kept as a list where the key holds one (or is new), e.g.
  // layers.0.crystal.orientation.surface=[1,1,0].
  try {
    parsed = json::parse(value);
    const bool list_for_list = parsed.is_array() && (!obj.contains(leaf) || obj[leaf].is_array());
    if ((parsed.is_object() || parsed.is_array()) && !list_for_list) parsed = json::value(value);
  } catch (const std::exception &) {
    parsed = json::value(value);
  }
  if (obj.contains(leaf) && obj[leaf].is_string()) parsed = json::value(value);
  obj[leaf] = parsed;
  config = parse_config(json::serialize(root));
}

// ---------------------------------------------------------------------------
// Examples and engine construction
// ---------------------------------------------------------------------------

RunConfig example_config(Mode mode) {
  RunConfig c;
  c.mode = mode;
  auto &t = c.trim;
  t.title = "Fe (1 MeV) into Fe";
  t.ion = trim::IonBeam::defaults(26, 1.0 * mega_electron_volt);
  trim::TargetLayer fe;
  fe.name = "Fe";
  fe.width = 1.2 * micrometer;
  fe.density = Element::by_z(26).solid_density();
  auto el = trim::TargetElement::defaults(26);
  el.displacement_energy = 40.0 * electron_volt; // ASTM E521 value for Fe
  el.lattice_binding_energy = 0.0 * electron_volt;
  fe.elements.push_back(el);
  t.layers.push_back(fe);
  t.ion_count = 1000;
  t.damage = trim::DamageCalculation::kinchin_pease;
  t.output_directory = "output";
  return c;
}

namespace {

unsigned isotope_mass_number(unsigned z, mass_quantity mass) {
  const auto &el = Element::by_z(z);
  if (std::abs(((mass - el.natural_mass()) / atomic_mass_unit).value()) < 0.02 ||
      std::abs(((mass - el.most_abundant_isotope_mass()) / atomic_mass_unit).value()) < 0.02) {
    return el.most_abundant_mass_number();
  }
  return static_cast<unsigned>(std::lround((mass / atomic_mass_unit).value()));
}

Simulation::Settings::ScreeningModel native_screening(const std::string &s) {
  using M = Simulation::Settings::ScreeningModel;
  if (s == "thomas-fermi" || s == "thomasfermi" || s == "tf") return M::ThomasFermi;
  if (s == "quickbasic" || s == "firsov") return M::QuickBasic;
  if (s == "zbl" || s == "universal") return M::ZBL;
  if (s == "moliere") return M::Moliere;
  if (s == "custom") return M::Custom;
  throw std::invalid_argument("native.screening_model: expected thomas-fermi, quickbasic, zbl, moliere or custom");
}

Simulation::Settings::ElectronicStoppingModel native_electronic(const std::string &s) {
  using M = Simulation::Settings::ElectronicStoppingModel;
  if (s == "zbl") return M::ZBL;
  if (s == "lss") return M::LSS;
  if (s == "custom" || s == "custom-power-law") return M::CustomPowerLaw;
  if (s == "bethe-bloch" || s == "bethe") return M::BetheBloch;
  if (s == "builtin") return M::Builtin;
  if (s == "srim") return M::SrimTable;
  throw std::invalid_argument(
      "native.electronic_stopping_model: expected zbl, lss, custom, bethe-bloch, builtin or srim");
}

} // namespace

std::unique_ptr<SimulationEngine> make_engine(const RunConfig &c) {
  if (c.mode == Mode::srim) {
    trim::TrimConfig t = c.trim;
    if (t.configuration_json.empty()) {
      RunConfig copy = c;
      copy.trim.particles.clear(); // (TRIM.DAT lists can be long)
      // The record says what the run does: the crystal transport takes none
      // of TRIM's bugs, whatever the scenario's switches say.
      copy.trim.bugs = c.trim.effective_bugs();
      copy.trim.single_precision_tables = c.trim.effective_single_precision_tables();
      t.configuration_json = to_json(copy);
    }
    return std::make_unique<trim::TrimSimulation>(std::move(t));
  }
  // Native engine.
  const auto &t = c.trim;
  const auto &nv = c.native;
  t.validate();
  std::vector<Layer> layers;
  length_quantity depth = 0.0 * angstrom;
  for (const auto &l : t.layers) {
    depth += l.width;
    Layer::material_vector material;
    std::vector<energy_quantity> eds;
    for (const auto &e : l.elements) {
      material.emplace_back(e.stoichiometry,
                            Particle::Properties(e.z, isotope_mass_number(e.z, e.mass)));
      eds.push_back(e.displacement_energy);
    }
    Layer layer(std::move(material), depth, l.density);
    layer.set_displacement_energies(std::move(eds));
    layers.push_back(std::move(layer));
  }
  const std::filesystem::path out = t.output_directory.empty() ? "output" : t.output_directory;
  const unsigned threads =
      t.threads != 0 ? t.threads : std::max(1u, std::thread::hardware_concurrency());
  const unsigned ion_a = nv.ion_mass_number != 0 ? nv.ion_mass_number
                                                 : isotope_mass_number(t.ion.z, t.ion.mass);
  Simulation::Settings settings(
      nv.electron_stop_energy, nv.electron ? 0 : t.ion.z, nv.electron ? 0 : ion_a,
      nv.damage_cascade, nv.ion_stop_energy, nv.displacement_energy, nv.angle_divisors,
      nv.flying_distances, nv.range > 0.0 * angstrom ? nv.range : depth, t.ion_count, nv.electron,
      t.ion.energy, threads, t.ion.angle, nv.bins, nv.bins, nv.bins, nv.bins, nv.bins,
      (out / "displacement_histogram.csv").string(), (out / "electronic_loss_histogram.csv").string(),
      (out / "nuclear_loss_histogram.csv").string(), (out / "total_loss_histogram.csv").string(),
      nv.per_bombardment_output, (out / "bombardments").string() + "/", nv.electron_screening_table,
      nv.mott_table, native_screening(nv.screening_model), (out / "range_histogram.csv").string(),
      native_electronic(nv.electronic_stopping_model), nv.custom_electronic_coefficient,
      nv.custom_electronic_exponent, nv.custom_recoil_electronic_coefficient,
      nv.custom_recoil_electronic_exponent, nv.custom_nuclear_screening_scale, t.stopping_table);
  return std::make_unique<Simulation>(settings, Volume(std::move(layers)));
}

} // namespace config
} // namespace amcset
