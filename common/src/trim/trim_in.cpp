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

/*! \file trim_in.cpp
 *  \brief TRIM.IN reader and writer.
 */

#include "trim/trim_in.h"

#include "element_data.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <limits>
#include <fstream>
#include <istream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace amcset {
namespace trim {

using namespace common;

namespace {

class LineReader {
public:
  explicit LineReader(std::istream &in) : in_(in) {}

  std::string next(const char *what) {
    std::string line;
    if (!std::getline(in_, line)) {
      throw std::runtime_error("TRIM.IN: unexpected end of file while reading " +
                               std::string(what));
    }
    ++number_;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return line;
  }
  std::vector<double> numbers(const char *what, std::size_t at_least) {
    const auto line = next(what);
    std::vector<double> v;
    std::istringstream ss(line);
    std::string tok;
    while (ss >> tok) {
      try {
        std::size_t used = 0;
        const double d = std::stod(tok, &used);
        if (used != tok.size()) break;
        v.push_back(d);
      } catch (const std::exception &) {
        break;
      }
    }
    if (v.size() < at_least) {
      throw std::runtime_error("TRIM.IN line " + std::to_string(number_) + ": expected " +
                               std::to_string(at_least) + " numbers for " + what + ", got '" +
                               line + "'");
    }
    return v;
  }
  std::size_t line_number() const noexcept { return number_; }

private:
  std::istream &in_;
  std::size_t number_ = 0;
};

/*!
 * A whole number of a TRIM.IN field within [lo, hi]. Converting a double
 * outside an integer type's range (a negative count, 1e30, NaN) is undefined
 * behaviour in C++, and a wrapped-around count once made the layer reader
 * index past its numbers; such values are input errors.
 */
template <class T>
T whole(double v, double lo, double hi, const std::string &what) {
  if (!(v >= lo && v <= hi) || v != std::floor(v)) {
    throw std::runtime_error("TRIM.IN: " + what + " must be a whole number from " +
                             std::to_string(static_cast<long long>(lo)) + " to " +
                             std::to_string(static_cast<long long>(hi)) + ", got " +
                             std::to_string(v));
  }
  return static_cast<T>(v);
}

//! Split off a quoted string at the start (after optional leading tokens).
std::string take_quoted(std::string &line) {
  const auto a = line.find('"');
  if (a == std::string::npos) return {};
  const auto b = line.find('"', a + 1);
  if (b == std::string::npos) return {};
  std::string name = line.substr(a + 1, b - a - 1);
  line = line.substr(0, a) + " " + line.substr(b + 1);
  return name;
}

std::vector<double> parse_numbers(const std::string &s) {
  std::vector<double> v;
  std::istringstream ss(s);
  double d = 0.0;
  while (ss >> d) v.push_back(d);
  return v;
}

} // namespace

TrimConfig parse_trim_in(std::istream &in, TrimInExtras *extras_out) {
  LineReader r(in);
  TrimConfig c;
  // TRIM computes its stopping from SRIM-2013's coefficients: use the
  // SRIM-2013 tables that come with AMCSET (empty table path).
  c.stopping_model = physics::ElectronicStoppingModel::srim_tables;
  TrimInExtras extras;

  r.next("the title line");
  r.next("the ion header");
  const auto ion = r.numbers("the ion data (Z1 M1 E angle N bragg autosave)", 5);
  c.ion.z = whole<unsigned>(ion[0], 1, Element::max_z, "the ion's Z");
  c.ion.mass = ion[1] * atomic_mass_unit;
  c.ion.energy = ion[2] * kilo_electron_volt;
  c.ion.angle = ion[3] * degree;
  c.ion_count = whole<std::size_t>(ion[4], 0, 1e12, "the number of ions");
  if (ion.size() > 6) extras.autosave = whole<int>(ion[6], -1e9, 1e9, "the autosave number");

  r.next("the cascade header");
  const auto casc = r.numbers("the cascade type, seed and reminders", 2);
  extras.cascade_code = whole<int>(casc[0], 1, 7, "the cascade type");
  c.seed = whole<std::uint32_t>(casc[1], 0, 4294967295.0, "the random number seed");
  if (casc.size() > 2) extras.reminders = whole<int>(casc[2], -1e9, 1e9, "the reminders flag");
  switch (extras.cascade_code) {
  case 1: c.damage = DamageCalculation::kinchin_pease; break;
  case 2: c.damage = DamageCalculation::full_cascade; break;
  case 3: c.damage = DamageCalculation::monolayer; break;
  // Codes 4-7 take their particles from TRIM.DAT (read_trim_in loads it from
  // the folder of TRIM.IN); TRIM runs every line and ignores the ion number.
  case 4:
    c.damage = DamageCalculation::kinchin_pease;
    c.source = TrimConfig::Source::particle_ions;
    break;
  case 5:
    c.damage = DamageCalculation::full_cascade;
    c.source = TrimConfig::Source::particle_ions;
    break;
  // TRIM.exe maps code 6 to its internal 7 and 7 to 8 (0x43acb0): 6 starts
  // a full cascade at each TRIM.DAT atom (0x449330); 7 runs each atom
  // through the ion loop, with monolayer steps.
  case 6:
    c.damage = DamageCalculation::full_cascade;
    c.source = TrimConfig::Source::particle_recoils;
    break;
  case 7:
    c.damage = DamageCalculation::monolayer;
    c.source = TrimConfig::Source::particle_recoils;
    break;
  default:
    throw std::runtime_error("TRIM.IN: unknown cascade type " +
                             std::to_string(extras.cascade_code) + " (expected 1-7)");
  }
  if (c.source != TrimConfig::Source::beam) c.ion_count = 0;

  r.next("the diskfiles header");
  const auto disk = r.numbers("the diskfile flags", 4);
  extras.ranges_file = disk[0] != 0;
  extras.backscatter_file = disk[1] != 0;
  extras.transmit_file = disk[2] != 0;
  extras.sputter_file = disk[3] != 0;
  if (disk.size() > 4) extras.collisions_file = whole<int>(disk[4], 0, 2, "the collisions-file flag");
  if (disk.size() > 5) extras.exyz_interval = whole<int>(disk[5], 0, 1e9, "the EXYZ energy step");
  c.record_ion_positions = extras.ranges_file;
  c.record_exiting_atoms = extras.backscatter_file || extras.transmit_file || extras.sputter_file;
  c.exyz_interval = static_cast<double>(std::max(extras.exyz_interval, 0)) * electron_volt;
  if (extras.collisions_file > 0) {
    c.collision_details_ions =
        c.source == TrimConfig::Source::beam ? c.ion_count : std::numeric_limits<std::size_t>::max();
    c.collision_recoil_details = extras.collisions_file >= 2;
  }

  r.next("the target header");
  std::string target_line = r.next("the target name and counts");
  c.title = take_quoted(target_line);
  const auto counts = parse_numbers(target_line);
  if (counts.size() < 2) {
    throw std::runtime_error("TRIM.IN line " + std::to_string(r.line_number()) +
                             ": expected the number of elements and layers");
  }
  const auto n_elements = whole<std::size_t>(counts[0], 1, 1000, "the number of target elements");
  const auto n_layers = whole<std::size_t>(counts[1], 1, 100000, "the number of layers");
  if (n_elements == 0 || n_layers == 0) throw std::runtime_error("TRIM.IN: empty target");

  r.next("the plot header");
  const auto plot = r.numbers("the plot type and depth window", 3);
  extras.plot_type = whole<int>(plot[0], -1e9, 1e9, "the plot type");
  c.window_start = plot[1] * angstrom;
  c.window_end = plot[2] * angstrom;

  r.next("the target elements header");
  struct Atom {
    unsigned z;
    double mass;
    std::string symbol;
  };
  std::vector<Atom> atoms;
  for (std::size_t i = 0; i < n_elements; ++i) {
    auto line = r.next("a target element");
    // "Atom 1 = Fe =      26  55.847"
    const auto eq = line.rfind('=');
    if (eq == std::string::npos) {
      throw std::runtime_error("TRIM.IN line " + std::to_string(r.line_number()) +
                               ": expected 'Atom n = X = Z mass'");
    }
    const auto v = parse_numbers(line.substr(eq + 1));
    if (v.size() < 2) {
      throw std::runtime_error("TRIM.IN line " + std::to_string(r.line_number()) +
                               ": expected Z and mass");
    }
    const auto z = whole<unsigned>(v[0], 1, Element::max_z, "a target element's Z");
    atoms.push_back({z, v[1], std::string(Element::by_z(z).symbol())});
  }
  // Elements appearing more than once get distinct names (TRIM columns).
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    for (std::size_t k = 0; k < i; ++k) {
      if (atoms[k].symbol == atoms[i].symbol) {
        atoms[i].symbol += "(" + std::to_string(i + 1) + ")";
        break;
      }
    }
  }

  r.next("the layer header (1)");
  r.next("the layer header (2)");
  struct LayerRow {
    std::string name;
    double width;
    double density;
    std::vector<double> stoich;
  };
  std::vector<LayerRow> rows;
  for (std::size_t i = 0; i < n_layers; ++i) {
    auto line = r.next("a layer");
    LayerRow row;
    row.name = take_quoted(line);
    const auto v = parse_numbers(line);
    if (v.size() < 3 + n_elements) {
      throw std::runtime_error("TRIM.IN line " + std::to_string(r.line_number()) +
                               ": expected layer number, width, density and " +
                               std::to_string(n_elements) + " stoichiometries");
    }
    row.width = v[1];
    row.density = v[2];
    row.stoich.assign(v.begin() + 3, v.begin() + 3 + static_cast<long>(n_elements));
    rows.push_back(row);
  }
  r.next("the layer phase header");
  const auto phases = r.numbers("the layer phases", n_layers);
  r.next("the compound correction header");
  const auto bragg = r.numbers("the compound corrections", n_layers);
  r.next("the displacement energy header");
  const auto ed = r.numbers("the displacement energies", n_elements);
  r.next("the lattice binding energy header");
  const auto eb = r.numbers("the lattice binding energies", n_elements);
  r.next("the surface binding energy header");
  const auto es = r.numbers("the surface binding energies", n_elements);
  try {
    r.next("the stopping version header");
    const auto sv = r.numbers("the stopping version", 1);
    extras.stopping_version = whole<int>(sv[0], -1e9, 1e9, "the stopping version");
  } catch (const std::runtime_error &) {
    // optional
  }

  for (std::size_t li = 0; li < n_layers; ++li) {
    TargetLayer layer;
    layer.name = rows[li].name;
    layer.width = rows[li].width * angstrom;
    layer.density = rows[li].density * gram_per_cubic_centimeter;
    layer.gas = phases[li] != 0.0;
    layer.compound_correction = bragg[li];
    for (std::size_t e = 0; e < n_elements; ++e) {
      if (!(rows[li].stoich[e] > 0.0)) continue;
      TargetElement el;
      el.z = atoms[e].z;
      el.mass = atoms[e].mass * atomic_mass_unit;
      el.stoichiometry = rows[li].stoich[e];
      el.displacement_energy = ed[e] * electron_volt;
      el.lattice_binding_energy = eb[e] * electron_volt;
      el.surface_binding_energy = es[e] * electron_volt;
      el.name = atoms[e].symbol;
      el.trim_atom = static_cast<unsigned>(e + 1);
      layer.elements.push_back(el);
    }
    c.layers.push_back(std::move(layer));
  }
  if (extras_out != nullptr) *extras_out = extras;
  return c;
}

std::filesystem::path find_trim_dat(const std::filesystem::path &directory) {
  std::error_code ec;
  for (const auto &entry : std::filesystem::directory_iterator(directory, ec)) {
    std::string name = entry.path().filename().string();
    for (char &ch : name) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
    if (name == "trim.dat" && entry.is_regular_file()) return entry.path();
  }
  return {};
}

TrimConfig read_trim_in(const std::filesystem::path &path, TrimInExtras *extras) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("Cannot open " + path.string());
  TrimConfig c;
  try {
    c = parse_trim_in(in, extras);
  } catch (const std::runtime_error &e) {
    throw std::runtime_error(path.string() + ": " + e.what());
  }
  if (c.source != TrimConfig::Source::beam) {
    auto dir = path.parent_path();
    if (dir.empty()) dir = ".";
    const auto dat = find_trim_dat(dir);
    if (dat.empty()) {
      throw std::runtime_error(path.string() + ": cascade type " +
                               std::to_string(c.trim_code()) +
                               " needs TRIM.DAT in the same folder");
    }
    auto list = read_trim_dat(dat);
    c.particles = std::move(list.particles);
    c.particles_title = list.title;
  }
  return c;
}

TrimInExtras extras_for(const TrimConfig &c) {
  TrimInExtras x;
  x.ranges_file = c.record_ion_positions;
  x.backscatter_file = x.transmit_file = x.sputter_file = c.record_exiting_atoms;
  x.collisions_file = c.collision_details_ions == 0 ? 0 : c.collision_recoil_details ? 2 : 1;
  x.exyz_interval = static_cast<int>(std::lround(in_units(c.exyz_interval, electron_volt)));
  return x;
}

void write_trim_in(std::ostream &o, const TrimConfig &c, const TrimInExtras &x) {
  // TRIM lists every distinct atom once; layers give a stoichiometry for each.
  struct Atom {
    std::string name;
    unsigned z;
    double mass;
    double ed, eb, es;
  };
  std::vector<Atom> atoms;
  std::vector<std::vector<double>> stoich(c.layers.size());
  for (const TargetElement *e : c.trim_atoms()) {
    atoms.push_back({e->name.empty() ? std::string(Element::by_z(e->z).symbol()) : e->name, e->z,
                     (e->mass / atomic_mass_unit).value(), in_units(e->displacement_energy, electron_volt),
                     in_units(e->lattice_binding_energy, electron_volt),
                     in_units(e->surface_binding_energy, electron_volt)});
  }
  for (std::size_t li = 0; li < c.layers.size(); ++li) {
    stoich[li].assign(atoms.size(), 0.0);
    double total = 0.0;
    for (const auto &e : c.layers[li].elements) total += e.stoichiometry;
    for (const auto &e : c.layers[li].elements) {
      const std::string name = e.name.empty() ? std::string(Element::by_z(e.z).symbol()) : e.name;
      for (std::size_t k = 0; k < atoms.size(); ++k) {
        if (atoms[k].name == name && atoms[k].z == e.z) stoich[li][k] += e.stoichiometry / total;
      }
    }
  }
  const int cascade = c.trim_code();
  char buf[512];
  o << "==> SRIM-2013.00 This file controls TRIM Calculations.\r\n";
  o << "Ion: Z1 ,  M1,  Energy (keV), Angle,Number,Bragg Corr,AutoSave Number.\r\n";
  std::snprintf(buf, sizeof buf, "     %u   %.3f   %g   %g   %zu   1   %zu\r\n", c.ion.z,
                (c.ion.mass / atomic_mass_unit).value(), in_units(c.ion.energy, kilo_electron_volt),
                in_units(c.ion.angle, degree), c.ion_count,
                x.autosave > 0 ? static_cast<std::size_t>(x.autosave) : c.ion_count);
  o << buf;
  o << "Cascades(1=No;2=Full;3=Sputt;4-5=Ions;6-7=Neutrons), Random Number Seed, Reminders\r\n";
  std::snprintf(buf, sizeof buf, "                      %d                                   %u       %d\r\n",
                cascade, c.seed, x.reminders);
  o << buf;
  o << "Diskfiles (0=no,1=yes): Ranges, Backscatt, Transmit, Sputtered, Collisions(1=Ion;2=Ion+Recoils), Special EXYZ.txt file\r\n";
  std::snprintf(buf, sizeof buf, "                          %d       %d       %d       %d       %d       %d\r\n",
                x.ranges_file ? 1 : 0, x.backscatter_file ? 1 : 0, x.transmit_file ? 1 : 0,
                x.sputter_file ? 1 : 0, x.collisions_file, x.exyz_interval);
  o << buf;
  o << "Target material : Number of Elements & Layers\r\n";
  o << "\"" << c.title << "\"       " << atoms.size() << "               " << c.layers.size() << "\r\n";
  o << "PlotType (0-5); Plot Depths: Xmin, Xmax(Ang.) [=0 0 for Viewing Full Target]\r\n";
  std::snprintf(buf, sizeof buf, "       %d                         %g            %g\r\n", x.plot_type,
                in_units(c.window_start, angstrom), in_units(c.window_end, angstrom));
  o << buf;
  o << "Target Elements:    Z   Mass(amu)\r\n";
  for (std::size_t k = 0; k < atoms.size(); ++k) {
    std::snprintf(buf, sizeof buf, "Atom %zu = %s =      %u  %.3f\r\n", k + 1,
                  std::string(Element::by_z(atoms[k].z).symbol()).c_str(), atoms[k].z, atoms[k].mass);
    o << buf;
  }
  o << "Layer   Layer Name /               Width Density";
  for (const auto &a : atoms) o << "    " << Element::by_z(a.z).symbol() << "(" << a.z << ")";
  o << "\r\nNumb.   Description                (Ang) (g/cm3)";
  for (std::size_t k = 0; k < atoms.size(); ++k) o << "    Stoich";
  o << "\r\n";
  for (std::size_t li = 0; li < c.layers.size(); ++li) {
    const auto &l = c.layers[li];
    std::snprintf(buf, sizeof buf, " %zu      \"%s\"           %g  %.4f", li + 1, l.name.c_str(),
                  in_units(l.width, angstrom), in_units(l.density, gram_per_cubic_centimeter));
    o << buf;
    for (double s : stoich[li]) {
      std::snprintf(buf, sizeof buf, "  %.6f", s);
      o << buf;
    }
    o << "\r\n";
  }
  o << "0  Target layer phases (0=Solid, 1=Gas)\r\n";
  for (const auto &l : c.layers) o << (l.gas ? "1 " : "0 ");
  o << "\r\nTarget Compound Corrections (Bragg)\r\n";
  for (const auto &l : c.layers) o << " " << l.compound_correction << "  ";
  o << "\r\nIndividual target atom displacement energies (eV)\r\n";
  for (const auto &a : atoms) {
    std::snprintf(buf, sizeof buf, "   %.2f", a.ed);
    o << buf;
  }
  o << "\r\nIndividual target atom lattice binding energies (eV)\r\n";
  for (const auto &a : atoms) {
    std::snprintf(buf, sizeof buf, "   %.2f", a.eb);
    o << buf;
  }
  o << "\r\nIndividual target atom surface binding energies (eV)\r\n";
  for (const auto &a : atoms) {
    std::snprintf(buf, sizeof buf, "   %.2f", a.es);
    o << buf;
  }
  o << "\r\nStopping Power Version (1=2011, 0=2011)\r\n " << x.stopping_version << " \r\n";
}

} // namespace trim
} // namespace amcset
