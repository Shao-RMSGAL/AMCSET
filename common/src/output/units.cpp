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

/*! \file units.cpp
 *  \brief Output unit registry built from the Boost.Units quantities.
 */

#include "output/units.h"

#include "amcset_units.h"

#include <algorithm>
#include <stdexcept>

namespace amcset {
namespace output {

namespace {

namespace bu = boost::units;
using namespace common;

/*
 * A unit is defined by a quantity of its dimension; the number stored is how
 * many of it make one SI base unit: (base / unit) is dimensionless, so the
 * compiler checks that both have the same dimension.
 */
template <class Q>
Unit make(std::string name, Dimension d, const Q &unit, const Q &base,
          std::vector<std::string> aliases = {}) {
  const bu::quantity<bu::si::dimensionless> ratio = base / unit;
  return Unit{std::move(name), d, ratio.value(), std::move(aliases)};
}

std::vector<Unit> build() {
  std::vector<Unit> u;
  const dimensionless_quantity one = 1.0 * bu::si::si_dimensionless;
  // Dimensionless: the first entry is the default.
  u.push_back(make("1", Dimension::dimensionless, one, one,
                   {"", "-", "none", "count", "atoms/ion", "/ion", "per ion"}));
  u.push_back(make("%", Dimension::dimensionless, 0.01 * one, one, {"percent"}));
  u.push_back(make("ppm", Dimension::dimensionless, 1e-6 * one, one));

  // Length (TRIM: Angstrom).
  const length_quantity m = meter;
  u.push_back(make("A", Dimension::length, angstrom, m, {"Angstrom", "angstrom", "\xC3\x85", "Ang", "ang"}));
  u.push_back(make("nm", Dimension::length, nanometer, m, {"nanometer"}));
  u.push_back(make("um", Dimension::length, micrometer, m, {"\xC2\xB5m", "micron", "micrometer"}));
  u.push_back(make("mm", Dimension::length, 1e-3 * meter, m, {"millimeter"}));
  u.push_back(make("cm", Dimension::length, centimeter, m, {"centimeter"}));
  u.push_back(make("m", Dimension::length, meter, m, {"meter"}));

  // Inverse length: densities per ion per unit depth.
  const inverse_length_quantity per_m = 1.0 / meter;
  u.push_back(make("1/A", Dimension::inverse_length, 1.0 / angstrom, per_m,
                   {"/A", "A^-1", "1/Angstrom", "/(Angstrom-Ion)", "atoms/(A ion)"}));
  u.push_back(make("1/nm", Dimension::inverse_length, 1.0 / nanometer, per_m, {"/nm", "nm^-1"}));
  u.push_back(make("1/um", Dimension::inverse_length, 1.0 / micrometer, per_m, {"/um", "um^-1"}));
  u.push_back(make("1/mm", Dimension::inverse_length, 1.0 / (1e-3 * meter), per_m, {"/mm", "mm^-1"}));
  u.push_back(make("1/cm", Dimension::inverse_length, 1.0 / centimeter, per_m,
                   {"/cm", "cm^-1", "(atoms/cm3)/(atoms/cm2)", "(atoms/cm^3)/(atoms/cm^2)"}));
  u.push_back(make("1/m", Dimension::inverse_length, per_m, per_m, {"/m", "m^-1"}));

  // Energy (TRIM: eV).
  const energy_quantity j = joule;
  u.push_back(make("eV", Dimension::energy, electron_volt, j, {"electronvolt"}));
  u.push_back(make("meV", Dimension::energy, 1e-3 * electron_volt, j));
  u.push_back(make("keV", Dimension::energy, kilo_electron_volt, j));
  u.push_back(make("MeV", Dimension::energy, mega_electron_volt, j));
  u.push_back(make("GeV", Dimension::energy, giga_electron_volt, j));
  u.push_back(make("J", Dimension::energy, joule, j, {"joule"}));

  // Energy per length (TRIM: eV/Angstrom).
  const stopping_power_quantity j_per_m = joule / meter;
  u.push_back(make("eV/A", Dimension::energy_per_length, electron_volt / angstrom, j_per_m,
                   {"eV/Angstrom", "eV/(Angstrom-Ion)", "eV/(A ion)", "eV/Ang"}));
  u.push_back(make("eV/nm", Dimension::energy_per_length, electron_volt / nanometer, j_per_m));
  u.push_back(make("keV/nm", Dimension::energy_per_length, kilo_electron_volt / nanometer, j_per_m));
  u.push_back(make("keV/um", Dimension::energy_per_length, kilo_electron_volt / micrometer, j_per_m,
                   {"keV/\xC2\xB5m", "keV/micron"}));
  u.push_back(make("MeV/mm", Dimension::energy_per_length, mega_electron_volt / (1e-3 * meter), j_per_m));
  u.push_back(make("MeV/cm", Dimension::energy_per_length, mega_electron_volt / centimeter, j_per_m));
  u.push_back(make("eV/cm", Dimension::energy_per_length, electron_volt / centimeter, j_per_m));
  u.push_back(make("J/m", Dimension::energy_per_length, j_per_m, j_per_m));

  // Inverse energy (spectra per unit energy).
  const inverse_energy_quantity per_j = 1.0 / joule;
  u.push_back(make("1/eV", Dimension::inverse_energy, 1.0 / electron_volt, per_j, {"/eV", "eV^-1"}));
  u.push_back(make("1/keV", Dimension::inverse_energy, 1.0 / kilo_electron_volt, per_j, {"/keV", "keV^-1"}));
  u.push_back(make("1/MeV", Dimension::inverse_energy, 1.0 / mega_electron_volt, per_j, {"/MeV", "MeV^-1"}));
  u.push_back(make("1/J", Dimension::inverse_energy, per_j, per_j, {"/J"}));

  // Angle.
  u.push_back(make("deg", Dimension::angle, degree, radian, {"degree", "degrees", "\xC2\xB0"}));
  u.push_back(make("rad", Dimension::angle, radian, radian, {"radian", "radians"}));
  u.push_back(make("mrad", Dimension::angle, 1e-3 * radian, radian));

  // Time.
  const time_quantity s = 1.0 * bu::si::second;
  u.push_back(make("s", Dimension::time, s, s, {"second", "seconds"}));
  u.push_back(make("ms", Dimension::time, 1e-3 * s, s));
  u.push_back(make("us", Dimension::time, 1e-6 * s, s, {"\xC2\xB5s"}));
  u.push_back(make("ns", Dimension::time, 1e-9 * s, s));
  u.push_back(make("ps", Dimension::time, 1e-12 * s, s));
  u.push_back(make("fs", Dimension::time, 1e-15 * s, s));
  return u;
}

const std::vector<Unit> &registry() {
  static const std::vector<Unit> r = build();
  return r;
}

std::string known_units(const Dimension *d) {
  std::string s;
  for (const auto &u : registry()) {
    if (d != nullptr && u.dimension != *d) continue;
    if (!s.empty()) s += ", ";
    s += u.name;
  }
  return s;
}

} // namespace

std::string to_string(Dimension d) {
  switch (d) {
  case Dimension::dimensionless: return "dimensionless";
  case Dimension::length: return "length";
  case Dimension::inverse_length: return "inverse_length";
  case Dimension::energy: return "energy";
  case Dimension::energy_per_length: return "energy_per_length";
  case Dimension::inverse_energy: return "inverse_energy";
  case Dimension::angle: return "angle";
  case Dimension::time: return "time";
  }
  return "dimensionless";
}

Dimension parse_dimension(std::string_view s) {
  for (auto d : {Dimension::dimensionless, Dimension::length, Dimension::inverse_length,
                 Dimension::energy, Dimension::energy_per_length, Dimension::inverse_energy,
                 Dimension::angle, Dimension::time}) {
    if (s == to_string(d)) return d;
  }
  if (s == "stopping" || s == "stopping_power") return Dimension::energy_per_length;
  throw std::invalid_argument("unknown dimension '" + std::string(s) +
                              "' (dimensionless, length, inverse_length, energy, "
                              "energy_per_length, inverse_energy, angle, time)");
}

const Unit &find_unit(std::string_view name) {
  for (const auto &u : registry()) {
    if (u.name == name) return u;
  }
  for (const auto &u : registry()) {
    if (std::find(u.aliases.begin(), u.aliases.end(), name) != u.aliases.end()) return u;
  }
  throw std::invalid_argument("unknown unit '" + std::string(name) + "'; known units: " +
                              known_units(nullptr));
}

const Unit &find_unit(std::string_view name, Dimension dimension) {
  const Unit *found = nullptr;
  try {
    found = &find_unit(name);
  } catch (const std::invalid_argument &) {
  }
  if (found == nullptr || found->dimension != dimension) {
    throw std::invalid_argument("'" + std::string(name) + "' is not a unit of " +
                                to_string(dimension) + "; use one of: " + known_units(&dimension));
  }
  return *found;
}

std::vector<const Unit *> units_of(Dimension d) {
  std::vector<const Unit *> v;
  for (const auto &u : registry()) {
    if (u.dimension == d) v.push_back(&u);
  }
  return v;
}

const Unit &trim_unit(Dimension d) { return *units_of(d).front(); }

const Unit &base_unit(Dimension d) {
  for (const auto &u : registry()) {
    if (u.dimension == d && u.per_base == 1.0) return u;
  }
  return trim_unit(d);
}

double conversion_factor(const Unit &from, const Unit &to) {
  if (from.dimension != to.dimension) {
    throw std::invalid_argument("cannot convert " + from.name + " (" + to_string(from.dimension) +
                                ") to " + to.name + " (" + to_string(to.dimension) + ")");
  }
  // value[to] = value[from] * (base per from)^-1 ... : x_from / per_base_from = x_base.
  return to.per_base / from.per_base;
}

double conversion_factor(std::string_view from, std::string_view to) {
  return conversion_factor(find_unit(from), find_unit(to));
}

} // namespace output
} // namespace amcset
