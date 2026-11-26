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

/*! \file cli.cpp
 *  \brief Command-line parsing, option application, progress bar and the
 *  stopping-power table.
 */

#include "data_files.h"
#include "cli.h"
#include "config.h"

#include "crystal/cell_io.h"
#include "crystal/crystal.h"
#include "element_data.h"
#include "output/table.h"
#include "trim/trim_dat.h"
#include "physics/scattering.h"
#include "physics/stopping.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#if defined(_WIN32)
#include <io.h>
#define AMCSET_ISATTY _isatty
#define AMCSET_FILENO _fileno
#else
#include <unistd.h>
#define AMCSET_ISATTY isatty
#define AMCSET_FILENO fileno
#endif

namespace amcset {
namespace cli {

using namespace common;

Options parse_arguments(int argc, const char *const *argv) {
  std::vector<std::string> args;
  for (int i = 1; i < argc; ++i) args.emplace_back(argv[i]);
  return parse_arguments(args);
}

Options parse_arguments(const std::vector<std::string> &args) {
  Options o;
  std::size_t i = 0;
  const auto value = [&](const std::string &flag) -> std::string {
    if (i + 1 >= args.size()) throw UsageError(flag + " needs a value");
    return args[++i];
  };
  const auto to_count = [](const std::string &flag, const std::string &v) -> std::size_t {
    try {
      if (v.empty() || !std::isdigit(static_cast<unsigned char>(v[0]))) throw std::invalid_argument(v);
      std::size_t used = 0;
      const auto n = std::stoull(v, &used);
      if (used != v.size()) throw std::invalid_argument(v);
      return static_cast<std::size_t>(n);
    } catch (const std::exception &) {
      throw UsageError(flag + ": expected a non-negative integer, got '" + v + "'");
    }
  };
  for (; i < args.size(); ++i) {
    const std::string &a = args[i];
    if (a == "-h" || a == "--help") o.help = true;
    else if (a == "-q" || a == "--quiet") o.quiet = true;
    else if (a == "-v" || a == "--verbose") o.verbose = true;
    else if (a == "--json") o.json = true;
    else if (a == "--csv") o.csv = true;
    else if (a == "--gas") o.gas = true;
    else if (a == "-o" || a == "--output") o.output = value(a);
    else if (a == "-n" || a == "--ions") o.ions = to_count(a, value(a));
    else if (a == "-j" || a == "--threads") o.threads = static_cast<unsigned>(to_count(a, value(a)));
    else if (a == "-s" || a == "--seed") {
      const std::string v = value(a);
      if (v == "random") {
        o.random_seed = true;
        o.seed.reset();
      } else {
        o.seed = static_cast<std::uint32_t>(to_count(a, v));
        o.random_seed = false;
      }
    }
    else if (a == "--mode") o.mode = value(a);
    else if (a == "--damage") o.damage = value(a);
    else if (a == "--stopping" || a == "--model") o.stopping = value(a);
    else if (a == "--stopping-table") o.stopping_table = value(a);
    else if (a == "--srim-data") o.srim_data = value(a);
    else if (a == "--set") o.overrides.push_back(value(a));
    else if (a == "--fix-trim-bug") o.fixed_bugs.push_back(value(a));
    else if (a == "--fix-trim-bugs") o.fixed_bugs.emplace_back("all");
    else if (a == "--fast" || a == "--fast-math") o.fast_math = true;
    else if (a == "--log") o.log_file = value(a);
    else if (a == "--auto-depth") o.auto_depth = true;
    else if (a == "--save-config") o.save_config = value(a);
    else if (a == "--ion") o.ion = value(a);
    else if (a == "--ion-mass") o.ion_mass = value(a);
    else if (a == "--target") o.target = value(a);
    else if (a == "--density") o.density = value(a);
    else if (a == "--from") o.from = value(a);
    else if (a == "--to") o.to = value(a);
    else if (a == "--per-decade") o.per_decade = static_cast<unsigned>(to_count(a, value(a)));
    else if (a == "--units") o.units = value(a);
    else if (a == "--formats" || a == "--format") o.formats = value(a);
    else if (a == "--output-spec") o.output_spec = value(a);
    else if (a == "--output-units") o.output_units = value(a);
    else if (a == "--number-format") o.number_format = value(a);
    else if (a == "--particles") o.particles = value(a);
    else if (a == "--source") o.source = value(a);
    else if (a == "--recoil") o.recoil = value(a);
    else if (a == "--exyz") o.exyz = value(a);
    else if (a == "--trajectory-file") o.trajectory_file = value(a);
    else if (a == "--trajectory-points") o.trajectory_points = value(a);
    else if (a == "--md-events") o.md_events = value(a);
    else if (a == "--md-ec") o.md_ec.push_back(value(a));
    else if (a == "--md-emin") o.md_emin = value(a);
    else if (a == "--md-handoff") o.md_handoff = true;
    else if (a == "--lammps") o.lammps = true;
    else if (a == "--correct-nickel") o.correct_nickel = true;
    else if (a == "--ntypes") o.ntypes = static_cast<unsigned>(to_count(a, value(a)));
    else if (a == "--list") o.list = true;
    else if (a == "--structure") o.structure = value(a);
    else if (a == "--sublattice") o.sublattices.push_back(value(a));
    else if (a == "--interstitial") o.interstitials.push_back(value(a));
    else if (a == "--a" || a == "--lattice-constant") o.lattice_a = value(a);
    else if (a == "--c-over-a") o.c_over_a = value(a);
    else if (a == "--from-format") o.cell_from_format = value(a);
    else if (a == "--to-format") o.cell_to_format = value(a);
    else if (a == "--elements") o.cell_elements = value(a);
    else if (a == "--surface") o.surface = value(a);
    else if (a == "--reference") o.reference = value(a);
    else if (a == "--twist") o.twist = value(a);
    else if (a == "--beam") o.beam = value(a);
    else if (a == "--tilt") o.tilt = value(a);
    else if (a.size() > 1 && a[0] == '-' && !(a.size() > 1 && std::isdigit(static_cast<unsigned char>(a[1])))) {
      throw UsageError("unknown option '" + a + "'");
    } else if (o.command.empty()) {
      o.command = a;
    } else {
      o.positional.push_back(a);
    }
  }
  if (o.per_decade == 0) throw UsageError("--per-decade must be positive");
  return o;
}

void apply_options(config::RunConfig &c, const Options &o) {
  // --fix-trim-bug(s) first, so that --set physics.trim_bugs.NAME=true can
  // keep one of TRIM's behaviours (e.g. zero_energy_defaults_skip_first_atom,
  // which lets the first target atom keep Ed/Eb/Es = 0).
  for (const auto &name : o.fixed_bugs) {
    if (name == "all") {
      c.trim.bugs.set_all(false);
      continue;
    }
    bool found = false;
    for (const auto &m : trim::TrimBugs::members()) {
      if (name == m.first) {
        c.trim.bugs.*(m.second) = false;
        found = true;
      }
    }
    if (name == "single_precision_tables") {
      c.trim.single_precision_tables = false;
      found = true;
    }
    if (!found) {
      std::string list;
      for (const auto &m : trim::TrimBugs::members()) list += std::string(" ") + m.first;
      throw UsageError("--fix-trim-bug: unknown bug '" + name + "'; known: all single_precision_tables" + list);
    }
  }
  for (const auto &ov : o.overrides) config::apply_override(c, ov);
  if (o.mode) c.mode = config::parse_mode(*o.mode);
  if (o.output) c.trim.output_directory = *o.output;
  if (o.ions) c.trim.ion_count = *o.ions;
  if (o.threads) c.trim.threads = *o.threads;
  if (o.seed) {
    c.trim.seed = *o.seed;
    c.trim.random_seed = false;
  }
  if (o.random_seed) c.trim.random_seed = true;
  if (o.damage) c.trim.damage = trim::parse_damage_calculation(*o.damage);
  if (o.stopping_table) {
    c.trim.stopping_table = *o.stopping_table;
    if (!o.stopping) c.trim.stopping_model = physics::ElectronicStoppingModel::srim_tables;
  }
  if (o.stopping) c.trim.stopping_model = physics::parse_electronic_stopping_model(*o.stopping);
  if (o.srim_data) c.trim.srim_data_directory = *o.srim_data;
  if (o.fast_math) c.trim.fast_math = true;
  if (o.correct_nickel) c.trim.stopping_corrections.nickel_ions = true;
  if (o.auto_depth) c.trim.auto_depth.enabled = true;
  if (o.trajectory_file) {
    c.trim.trajectory_export.file = *o.trajectory_file;
    // The file needs recorded trajectories: every ion unless a number is set.
    if (c.trim.trajectory_ions == 0) c.trim.trajectory_ions = trim::TrimConfig::all_ions;
  }
  if (o.trajectory_points) {
    if (*o.trajectory_points == "all") c.trim.trajectory_export.births_only = false;
    else if (*o.trajectory_points == "births") c.trim.trajectory_export.births_only = true;
    else throw UsageError("--trajectory-points: all or births, got '" + *o.trajectory_points + "'");
  }
  if (o.md_events) c.trim.md_events.file = *o.md_events;
  if (!o.md_ec.empty()) {
    c.trim.md_events.transition_energies.clear();
    for (const auto &list : o.md_ec) {
      for (const auto &e : config::parse_energy_list(list, "--md-ec")) c.trim.md_events.transition_energies.push_back(e);
    }
  }
  if (o.md_emin) {
    c.trim.md_events.min_energy = (*o.md_emin == "Ed" || *o.md_emin == "ed")
                                      ? 0.0 * electron_volt
                                      : config::parse_energy(*o.md_emin, "--md-emin");
  }
  if (o.md_handoff) c.trim.md_events.handoff = true;
  auto &spec = c.trim.outputs;
  if (o.output_spec) {
    std::ifstream in(*o.output_spec, std::ios::binary);
    if (!in) throw std::runtime_error("Cannot open " + *o.output_spec);
    std::stringstream ss;
    ss << in.rdbuf();
    try {
      spec = output::parse_output_spec(ss.str());
    } catch (const std::invalid_argument &e) {
      throw std::invalid_argument(*o.output_spec + ": " + e.what());
    }
  }
  const auto split = [](const std::string &text, char sep) {
    std::vector<std::string> parts;
    std::string cur;
    for (char ch : text) {
      if (ch == sep) {
        parts.push_back(cur);
        cur.clear();
      } else if (!std::isspace(static_cast<unsigned char>(ch))) {
        cur += ch;
      }
    }
    if (!cur.empty()) parts.push_back(cur);
    return parts;
  };
  if (o.formats) {
    spec.formats.clear();
    for (const auto &f : split(*o.formats, ',')) {
      if (!f.empty()) spec.formats.insert(output::parse_format(f));
    }
  }
  if (o.output_units) {
    for (const auto &kv : split(*o.output_units, ',')) {
      const auto eq = kv.find('=');
      if (eq == std::string::npos) {
        // A bare unit sets its own dimension, e.g. "nm" or "keV".
        const auto &u = output::find_unit(kv);
        spec.units[u.dimension] = u.name;
        continue;
      }
      const auto d = output::parse_dimension(kv.substr(0, eq));
      const auto &u = output::find_unit(kv.substr(eq + 1), d);
      spec.units[d] = u.name;
    }
  }
  if (o.number_format) {
    const auto colon = o.number_format->find(':');
    spec.format.notation = output::parse_notation(o.number_format->substr(0, colon));
    if (colon != std::string::npos) {
      try {
        spec.format.digits = std::stoi(o.number_format->substr(colon + 1));
      } catch (const std::exception &) {
        throw UsageError("--number-format: expected NOTATION[:DIGITS], e.g. scientific:4");
      }
    }
  }
  if (o.source) {
    const auto &v = *o.source;
    if (v == "beam") c.trim.source = trim::TrimConfig::Source::beam;
    else if (v == "ions") c.trim.source = trim::TrimConfig::Source::particle_ions;
    else if (v == "recoils") c.trim.source = trim::TrimConfig::Source::particle_recoils;
    else throw UsageError("--source: expected beam, ions or recoils");
  }
  if (o.particles) {
    const std::filesystem::path p(*o.particles);
    std::string name = p.filename().string();
    for (char &ch : name) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
    const bool exits = name == "BACKSCAT.TXT" || name == "TRANSMIT.TXT" || name == "SPUTTER.TXT" ||
                       name == "TRIMOUT.TXT";
    auto list = exits ? trim::read_exit_file_as_particles(p) : trim::read_trim_dat(p);
    c.trim.particles = std::move(list.particles);
    c.trim.particles_title = list.title;
    if (c.trim.source == trim::TrimConfig::Source::beam) {
      c.trim.source = exits ? trim::TrimConfig::Source::particle_recoils
                            : trim::TrimConfig::Source::particle_ions;
    }
    if (!o.ions) c.trim.ion_count = 0;
  }
  if (o.recoil) {
    // --recoil ELEMENT,ENERGY,DEPTH[,COSX,COSY,COSZ]: one target atom set in
    // motion at DEPTH (TRIM's x; lateral 0), run as many times as there are
    // ions (each with its own random stream).
    std::vector<std::string> parts;
    std::string part;
    for (const char ch : *o.recoil + ",") {
      if (ch == ',') {
        parts.push_back(part);
        part.clear();
      } else {
        part += ch;
      }
    }
    if (parts.size() != 3 && parts.size() != 6) {
      throw UsageError("--recoil: expected ELEMENT,ENERGY,DEPTH[,COSX,COSY,COSZ], e.g. Fe,10keV,50nm");
    }
    trim::SourceParticle sp;
    sp.name = "PKA";
    try {
      sp.z = Element::by_symbol(parts[0]).z();
    } catch (const std::exception &) {
      throw UsageError("--recoil: unknown element '" + parts[0] + "'");
    }
    sp.energy = config::parse_energy(parts[1], "--recoil energy");
    sp.x = config::parse_length(parts[2], "--recoil depth");
    sp.y = 0.0 * angstrom;
    sp.z_position = 0.0 * angstrom;
    if (parts.size() == 6) {
      double v[3];
      for (int k = 0; k < 3; ++k) {
        try {
          std::size_t used = 0;
          v[k] = std::stod(parts[3 + static_cast<std::size_t>(k)], &used);
          if (used != parts[3 + static_cast<std::size_t>(k)].size()) throw std::invalid_argument("trailing");
        } catch (const std::exception &) {
          throw UsageError("--recoil: the direction is three numbers (cosines along depth, y, z)");
        }
      }
      sp.cos_x = v[0];
      sp.cos_y = v[1];
      sp.cos_z = v[2];
    }
    c.trim.source = trim::TrimConfig::Source::particle_recoils;
    c.trim.particles = {sp};
    c.trim.particles_title = "PKA from --recoil";
    c.trim.particle_repeat = std::max<std::size_t>(c.trim.ion_count, 1);
    c.trim.ion_count = 0;
    // (TRIM's ion energy is the highest particle energy the stopping tables cover.)
    if (c.trim.ion.energy < sp.energy) c.trim.ion.energy = sp.energy;
  }
  if (o.exyz) c.trim.exyz_interval = config::parse_energy(*o.exyz, "--exyz");
}

bool looks_like_trim_in(const std::string &path) {
  const auto name = std::filesystem::path(path).filename().string();
  std::string upper;
  for (char ch : name) upper.push_back(static_cast<char>(std::toupper(static_cast<unsigned char>(ch))));
  if (upper == "TRIM.IN" || (upper.size() > 3 && upper.substr(upper.size() - 3) == ".IN")) return true;
  std::ifstream in(path);
  std::string first;
  if (in && std::getline(in, first)) {
    return first.find("TRIM") != std::string::npos && first.find('{') == std::string::npos;
  }
  return false;
}

std::vector<std::pair<unsigned, double>> parse_composition(const std::string &spec) {
  std::vector<std::pair<unsigned, double>> out;
  std::stringstream ss(spec);
  std::string part;
  double total = 0.0;
  while (std::getline(ss, part, ',')) {
    const auto colon = part.find(':');
    const std::string sym = part.substr(0, colon);
    double frac = 1.0;
    if (colon != std::string::npos) {
      try {
        frac = std::stod(part.substr(colon + 1));
      } catch (const std::exception &) {
        throw std::invalid_argument("Bad fraction in '" + part + "'");
      }
    }
    if (!(frac > 0.0)) throw std::invalid_argument("Fractions must be positive in '" + spec + "'");
    out.emplace_back(Element::by_symbol(sym).z(), frac);
    total += frac;
  }
  if (out.empty()) throw std::invalid_argument("Empty composition");
  for (auto &p : out) p.second /= total;
  return out;
}

namespace {

//! An ion given by symbol or Z.
unsigned ion_z(const std::string &s) {
  try {
    std::size_t used = 0;
    const auto z = std::stoul(s, &used);
    if (used == s.size()) return static_cast<unsigned>(z);
  } catch (const std::exception &) {
  }
  return Element::by_symbol(s).z();
}

} // namespace

int print_stopping_table(const Options &o, std::ostream &out) {
  if (o.ion.empty() || o.target.empty()) throw UsageError("stopping needs --ion and --target");
  // --lammps takes a list of ions (one column per LAMMPS atom type).
  std::vector<std::string> ion_names;
  {
    std::string cur;
    for (char ch : o.ion + ",") {
      if (ch == ',') {
        if (!cur.empty()) ion_names.push_back(cur);
        cur.clear();
      } else if (!std::isspace(static_cast<unsigned char>(ch))) {
        cur += ch;
      }
    }
  }
  if (ion_names.empty()) throw UsageError("--ion: no projectile given");
  if (ion_names.size() > 1 && !o.lammps) throw UsageError("--ion: one projectile (a list only with --lammps)");
  if (o.ntypes > 0 && (!o.lammps || ion_names.size() != 1)) {
    throw UsageError("--ntypes: only with --lammps and a single --ion");
  }
  if (ion_names.size() > 1 && !o.ion_mass.empty()) throw UsageError("--ion-mass: only with a single --ion");
  const unsigned z1 = ion_z(ion_names.front());
  const auto &ion = Element::by_z(z1);
  const mass_quantity m1 = o.ion_mass.empty() ? ion.most_abundant_isotope_mass()
                                              : config::parse_mass(o.ion_mass, "--ion-mass");
  const auto comp = parse_composition(o.target);
  mass_quantity mean_mass = 0.0 * kilogram;
  std::vector<physics::StoppingComponent> layer;
  for (const auto &[z, f] : comp) {
    layer.push_back({z, f});
    mean_mass += f * Element::by_z(z).natural_mass();
  }
  mass_density_quantity density = 0.0 * gram_per_cubic_centimeter;
  if (!o.density.empty()) {
    density = config::parse_density(o.density, "--density");
  } else if (comp.size() == 1) {
    const auto &e = Element::by_z(comp.front().first);
    density = o.gas ? e.gas_density() : e.solid_density();
  } else {
    throw UsageError("--density is required for compounds");
  }
  const number_density_quantity n = density / mean_mass;
  physics::StoppingCorrections corrections;
  corrections.nickel_ions = o.correct_nickel;

  // Default: the SRIM-2013 tables that come with AMCSET, if they are found.
  const std::string bundled = data::bundled_srim_stopping_table();
  physics::ElectronicStoppingModel model = (o.stopping_table || !bundled.empty())
                                               ? physics::ElectronicStoppingModel::srim_tables
                                               : physics::ElectronicStoppingModel::builtin;
  if (o.stopping) model = physics::parse_electronic_stopping_model(*o.stopping);
  std::shared_ptr<const physics::SrimStoppingDatabase> db;
  if (model == physics::ElectronicStoppingModel::srim_tables) {
    const std::string path = o.stopping_table ? *o.stopping_table : bundled;
    if (path.empty()) throw UsageError("--model srim: the bundled SRIM tables were not found; give --stopping-table");
    db = physics::SrimStoppingDatabase::load(path);
  }
  const energy_quantity lo = config::parse_energy(o.from, "--from");
  const energy_quantity hi = config::parse_energy(o.to, "--to");
  if (!(hi >= lo) || !(lo > 0.0 * electron_volt)) throw UsageError("--from must be positive and not above --to");
  if (o.lammps) {
    // LAMMPS fix electron/stopping (metal units): the kinetic energy in eV,
    // then the electronic stopping N S_e in eV/Angstrom of each atom type,
    // in increasing energy. Electronic stopping only; Bragg's rule with the
    // chosen model, as SRIM mode uses it.
    struct Type {
      unsigned z;
      mass_quantity m;
    };
    std::vector<Type> types;
    for (const auto &name : ion_names) {
      const unsigned z = ion_z(name);
      types.push_back({z, types.empty() ? m1 : Element::by_z(z).most_abundant_isotope_mass()});
    }
    if (o.ntypes > 1) types.assign(o.ntypes, types.front());
    out << "# AMCSET " << AMCSET_VERSION << ": electronic stopping for LAMMPS fix electron/stopping (metal units)\n";
    out << "# target " << o.target << ", " << in_units(density, gram_per_cubic_centimeter) << " g/cm3, "
        << in_units(n, per_cubic_angstrom) << " atoms/A3, model " << physics::to_string(model) << "\n";
    out << "# column 1: kinetic energy (eV); then N*S_e (eV/Angstrom) of atom type";
    for (std::size_t k = 0; k < types.size(); ++k) {
      out << " " << k + 1 << "=" << Element::by_z(types[k].z).symbol() << " (" << (types[k].m / atomic_mass_unit).value()
          << " amu)";
    }
    out << "\n";
    const double decades = std::log10((hi / lo).value());
    const auto steps = static_cast<std::size_t>(std::ceil(decades * o.per_decade));
    for (std::size_t i = 0; i <= steps; ++i) {
      const energy_quantity e = lo * std::pow(10.0, std::min(decades, static_cast<double>(i) / o.per_decade));
      char buf[64];
      std::snprintf(buf, sizeof buf, "%.8g", in_units(e, electron_volt));
      out << buf;
      for (const auto &t : types) {
        const auto se = n * physics::layer_stopping(model, db.get(), t.z, t.m, layer, o.gas, 1.0, e, n, {}, corrections);
        std::snprintf(buf, sizeof buf, " %.8g", in_units(se, electron_volt_per_angstrom));
        out << buf;
      }
      out << "\n";
    }
    return 0;
  }
  std::shared_ptr<const physics::PairPotentialTable> pots;
  if (o.srim_data) pots = physics::PairPotentialTable::load_srim(*o.srim_data);
  std::vector<physics::CollisionPair> pairs;
  for (const auto &[z, f] : comp) {
    pairs.emplace_back(z1, m1, z, Element::by_z(z).natural_mass(), pots.get());
  }

  // Output unit.
  stopping_power_quantity unit = electron_volt_per_angstrom;
  bool per_areal = false;
  bool per_mass = false;
  std::string unit_name = o.units;
  if (o.units == "eV/A") unit = electron_volt_per_angstrom;
  else if (o.units == "keV/um") unit = kilo_electron_volt_per_micrometer;
  else if (o.units == "MeV/(mg/cm2)") per_mass = true;
  else if (o.units == "eV/(1e15 atoms/cm2)") per_areal = true;
  else throw UsageError("--units: expected eV/A, keV/um, MeV/(mg/cm2) or eV/(1e15 atoms/cm2)");

  const auto convert = [&](physics::stopping_cross_section_quantity s) {
    if (per_areal) return in_units(s, srim_stopping_unit);
    if (per_mass) {
      const auto per_mass_unit = mega_electron_volt / (1e-3 * gram / square_centimeter);
      return ((s * n / density) / per_mass_unit).value();
    }
    return in_units(s * n, unit);
  };

  const energy_quantity e0 = config::parse_energy(o.from, "--from");
  const energy_quantity e1 = config::parse_energy(o.to, "--to");
  if (!(e1 >= e0) || !(e0 > 0.0 * electron_volt)) {
    throw UsageError("--from must be positive and not above --to");
  }
  const double decades = std::log10((e1 / e0).value());
  const auto steps = static_cast<std::size_t>(std::ceil(decades * o.per_decade));
  const char sep = o.csv ? ',' : ' ';
  out << "# AMCSET stopping table: " << ion.symbol() << " (" << (m1 / atomic_mass_unit).value()
      << " amu) in " << o.target << ", density "
      << in_units(density, gram_per_cubic_centimeter) << " g/cm3, model "
      << physics::to_string(model) << (pots ? ", SRIM pair nuclear stopping" : "") << "\n";
  out << "# energy_keV" << sep << "S_electronic" << sep << "S_nuclear" << sep << "S_total" << sep
      << "range_um  (stopping in " << unit_name << "; range = continuous-slowing-down path length)\n";
  // CSDA range by integrating 1/S from 10 eV (velocity-proportional below).
  double range_um = 0.0;
  energy_quantity prev_e = 10.0 * electron_volt;
  auto total_power = [&](energy_quantity e) {
    stopping_power_quantity s =
        n * physics::layer_stopping(model, db.get(), z1, m1, layer, o.gas, 1.0, e, n, {}, corrections);
    for (std::size_t k = 0; k < pairs.size(); ++k) s += n * comp[k].second * pairs[k].nuclear_stopping(e);
    return s;
  };
  // Integrate up to the first tabulated energy.
  {
    const int sub = 200;
    for (int k = 1; k <= sub; ++k) {
      const energy_quantity e = prev_e * std::pow((e0 / prev_e).value(), static_cast<double>(k) / sub);
      const energy_quantity em = prev_e * std::pow((e0 / prev_e).value(), (k - 0.5) / sub);
      const energy_quantity ep = prev_e * std::pow((e0 / prev_e).value(), (k - 1.0) / sub);
      range_um += in_units((e - ep) / total_power(em), micrometer);
    }
    prev_e = e0;
  }
  for (std::size_t i = 0; i <= steps; ++i) {
    const energy_quantity e = e0 * std::pow(10.0, std::min(decades, static_cast<double>(i) / o.per_decade));
    if (i > 0) {
      const int sub = 20;
      for (int k = 1; k <= sub; ++k) {
        const energy_quantity ep = prev_e * std::pow((e / prev_e).value(), (k - 1.0) / sub);
        const energy_quantity eq = prev_e * std::pow((e / prev_e).value(), static_cast<double>(k) / sub);
        const energy_quantity em = prev_e * std::pow((e / prev_e).value(), (k - 0.5) / sub);
        range_um += in_units((eq - ep) / total_power(em), micrometer);
      }
    }
    prev_e = e;
    const auto se = physics::layer_stopping(model, db.get(), z1, m1, layer, o.gas, 1.0, e, n, {}, corrections);
    physics::stopping_cross_section_quantity sn = 0.0 * srim_stopping_unit;
    for (std::size_t k = 0; k < pairs.size(); ++k) sn += comp[k].second * pairs[k].nuclear_stopping(e);
    char buf[256];
    if (o.csv) {
      std::snprintf(buf, sizeof buf, "%.6g,%.6g,%.6g,%.6g,%.6g\n", in_units(e, kilo_electron_volt),
                    convert(se), convert(sn), convert(se + sn), range_um);
    } else {
      std::snprintf(buf, sizeof buf, "%12.5E %12.5E %12.5E %12.5E %12.5E\n",
                    in_units(e, kilo_electron_volt), convert(se), convert(sn), convert(se + sn), range_um);
    }
    out << buf;
  }
  return 0;
}

namespace {
double now_seconds() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
} // namespace

ProgressBar::ProgressBar(std::size_t total, bool enabled)
    : total_(total), enabled_(enabled && AMCSET_ISATTY(AMCSET_FILENO(stderr))),
      start_(now_seconds()) {}

void ProgressBar::update(std::size_t done) {
  if (!enabled_ || total_ == 0) return;
  const int permille = static_cast<int>(1000.0 * static_cast<double>(done) / static_cast<double>(total_));
  if (permille == last_permille_) return;
  last_permille_ = permille;
  const double elapsed = now_seconds() - start_;
  const double eta = done > 0 ? elapsed * (static_cast<double>(total_) / static_cast<double>(done) - 1.0) : 0.0;
  const int width = 40;
  const int filled = permille * width / 1000;
  std::string bar(static_cast<std::size_t>(filled), '#');
  bar.resize(static_cast<std::size_t>(width), '-');
  std::fprintf(stderr, "\r[%s] %5.1f%%  %zu/%zu  elapsed %.0fs  eta %.0fs   ", bar.c_str(),
               permille / 10.0, done, total_, elapsed, eta);
  std::fflush(stderr);
}

void ProgressBar::clear_line() {
  if (!enabled_) return;
  std::fprintf(stderr, "\r%90s\r", "");
  last_permille_ = -1;
}

void ProgressBar::finish(std::size_t done) {
  if (!enabled_) return;
  last_permille_ = -1;
  update(done);
  std::fprintf(stderr, "\n");
}

// ---------------------------------------------------------------------------
// cell: unit cells
// ---------------------------------------------------------------------------

namespace {

std::string lower_text(std::string s) {
  for (auto &c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

//! "Fe" or "Fe:0.9,Cr:0.1": elements with their site fractions (not normalized).
std::vector<crystal::Occupant> parse_occupant_spec(const std::string &spec, const std::string &flag) {
  std::vector<crystal::Occupant> out;
  std::stringstream ss(spec);
  std::string part;
  while (std::getline(ss, part, ',')) {
    const auto colon = part.find(':');
    std::string sym = part.substr(0, colon);
    sym.erase(std::remove_if(sym.begin(), sym.end(), [](unsigned char c) { return std::isspace(c); }), sym.end());
    const auto z = Element::find_z(sym);
    if (!z) throw std::invalid_argument(flag + ": unknown element '" + sym + "' in '" + spec + "'");
    double fraction = 1.0;
    if (colon != std::string::npos) {
      try {
        std::size_t used = 0;
        const std::string f = part.substr(colon + 1);
        fraction = std::stod(f, &used);
        if (used != f.size()) throw std::invalid_argument(f);
      } catch (const std::exception &) {
        throw std::invalid_argument(flag + ": bad fraction in '" + part + "'");
      }
    }
    out.push_back({*z, fraction});
  }
  if (out.empty()) throw std::invalid_argument(flag + ": empty");
  return out;
}

//! "1 1 0", "1,-1,0" or "[1 -1 0]": three or four numbers.
std::vector<double> parse_indices(const std::string &text, const std::string &flag) {
  std::string t = text;
  for (auto &c : t) {
    if (c == ',' || c == '[' || c == ']' || c == '(' || c == ')' || c == '<' || c == '>') c = ' ';
  }
  std::vector<double> out;
  std::stringstream ss(t);
  std::string tok;
  while (ss >> tok) {
    try {
      std::size_t used = 0;
      out.push_back(std::stod(tok, &used));
      if (used != tok.size()) throw std::invalid_argument(tok);
    } catch (const std::exception &) {
      throw std::invalid_argument(flag + ": '" + tok + "' is not a number (expected e.g. \"1 1 0\")");
    }
  }
  if (out.size() != 3 && out.size() != 4) {
    throw std::invalid_argument(flag + ": expected three indices (four for hexagonal cells), e.g. \"1 1 0\"");
  }
  return out;
}

std::string fixed(double v, int decimals) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.*f", decimals, std::abs(v) < 0.5 * std::pow(10.0, -decimals) ? 0.0 : v);
  return buf;
}

std::string general(double v, int digits = 6) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.*g", digits, v);
  return buf;
}

std::string occupants_text(const crystal::Site &s) {
  if (s.occupants.empty()) return "(empty)";
  std::string t;
  for (const auto &o : s.occupants) {
    t += (t.empty() ? "" : " + ") + std::string(Element::by_z(o.z).symbol()) + " " + general(o.fraction, 6);
  }
  return t;
}

//! The low-index crystal direction of a vector, e.g. "[1 1 0]" or "[1 1 4] (0.3 deg off)".
std::string direction_text(const crystal::Lattice &lattice, const crystal::Vec3 &v) {
  const auto n = crystal::nearest_direction(lattice, v);
  std::string t = crystal::format_indices({double(n.uvw[0]), double(n.uvw[1]), double(n.uvw[2])});
  const double off = in_units(n.off, degree);
  if (off > 0.005) t += " (" + fixed(off, 2) + " deg off)";
  return t;
}

} // namespace

int print_cell(const Options &o, std::ostream &out) {
  if (o.list) {
    out << "Built-in crystal structures (amcset cell --structure NAME --sublattice ...):\n";
    for (const auto &s : crystal::builtin_structures()) {
      std::string subs;
      for (auto n : s.sublattices) subs += (subs.empty() ? "" : "; ") + std::string(n);
      std::string ins;
      for (auto n : s.interstitials) ins += (ins.empty() ? "" : ", ") + std::string(n);
      char buf[512];
      std::snprintf(buf, sizeof buf, "  %-16s %-28s %s\n", std::string(s.name).c_str(), std::string(s.title).c_str(),
                    std::string(s.examples).c_str());
      out << buf;
      out << "                   sublattices: " << subs << (s.hexagonal || s.tetragonal ? "; needs c/a" : "")
          << (s.has_u ? "; internal parameter u" : "") << "\n";
      if (!ins.empty()) out << "                   interstitial sites: " << ins << "\n";
    }
    return 0;
  }

  crystal::Crystal c;
  bool from_density = false;
  if (!o.structure.empty()) {
    if (!o.positional.empty()) throw UsageError("cell: give a cell file or --structure, not both");
    crystal::StructureSpec spec;
    spec.name = o.structure;
    for (const auto &s : o.sublattices) spec.sublattices.push_back(parse_occupant_spec(s, "--sublattice"));
    for (const auto &s : o.interstitials) {
      const auto eq = s.find('=');
      if (eq == std::string::npos) {
        throw UsageError("--interstitial: expected SITES=ELEMENTS, e.g. octahedral=C:0.01");
      }
      spec.interstitials.push_back({lower_text(s.substr(0, eq)), parse_occupant_spec(s.substr(eq + 1), "--interstitial")});
    }
    if (!o.lattice_a.empty()) {
      spec.a_A = in_units(config::parse_length(o.lattice_a, "--a"), angstrom);
      if (!(spec.a_A > 0.0)) throw std::invalid_argument("--a: must be positive");
    } else {
      from_density = true;
    }
    if (!o.c_over_a.empty()) {
      try {
        spec.c_over_a = std::stod(o.c_over_a);
      } catch (const std::exception &) {
        throw UsageError("--c-over-a: expected a number");
      }
    }
    c.cell = crystal::build_structure(spec);
  } else if (o.positional.size() == 1) {
    for (const auto *given : {&o.lattice_a, &o.c_over_a}) {
      if (!given->empty()) throw UsageError("cell: --a and --c-over-a belong to --structure");
    }
    if (!o.sublattices.empty() || !o.interstitials.empty()) {
      throw UsageError("cell: --sublattice and --interstitial belong to --structure");
    }
    std::optional<crystal::CellFormat> format;
    if (!o.cell_from_format.empty()) format = crystal::parse_cell_format(o.cell_from_format);
    crystal::CellReadOptions opt;
    if (!o.cell_elements.empty()) {
      for (const auto &occ : parse_occupant_spec(o.cell_elements, "--elements")) opt.elements.push_back(occ.z);
    }
    c.cell = crystal::read_cell(o.positional.front(), format, opt);
  } else {
    throw UsageError("cell needs one cell file or --structure NAME (amcset cell --list shows the structures)");
  }
  if (!o.surface.empty()) c.orientation.surface = parse_indices(o.surface, "--surface");
  if (!o.reference.empty()) c.orientation.reference = parse_indices(o.reference, "--reference");
  if (!o.twist.empty()) c.orientation.twist = config::parse_angle(o.twist, "--twist");

  // A layer makes the composition and density the cell's and checks it all.
  trim::TargetLayer layer;
  layer.width = 1.0 * micrometer;
  layer.crystal = c;
  if (from_density) {
    const auto composition = c.cell.composition();
    mass_density_quantity rho = 0.0 * gram_per_cubic_centimeter;
    if (!o.density.empty()) {
      rho = config::parse_density(o.density, "--density");
    } else if (composition.size() == 1) {
      rho = Element::by_z(composition.begin()->first).solid_density();
    }
    if (!(rho > 0.0 * gram_per_cubic_centimeter)) {
      throw std::invalid_argument("cell: give the lattice constant (--a) or the density (--density)");
    }
    layer.crystal->cell.scale_to_density(rho, {});
  } else if (!o.density.empty()) {
    throw UsageError("cell: --density sizes a built-in structure that has no --a");
  }
  trim::apply_crystal(layer);
  trim::check_layer_crystal(layer, "cell");
  const crystal::Crystal &cr = *layer.crystal;
  const crystal::UnitCell &cell = cr.cell;

  if (o.output) {
    std::optional<crystal::CellFormat> to;
    if (!o.cell_to_format.empty()) to = crystal::parse_cell_format(o.cell_to_format);
    crystal::save_cell(cell, *o.output, to);
  }
  if (o.json) {
    out << config::crystal_to_json(cr);
    return 0;
  }

  const auto &lat = cell.lattice;
  out << "Unit cell: " << (cell.name.empty() ? std::string("(unnamed)") : cell.name) << "\n";
  out << "  a = " << fixed(in_units(lat.a(), angstrom), 5) << " A   b = " << fixed(in_units(lat.b(), angstrom), 5)
      << " A   c = " << fixed(in_units(lat.c(), angstrom), 5) << " A\n";
  out << "  alpha = " << fixed(in_units(lat.alpha(), degree), 4) << "   beta = "
      << fixed(in_units(lat.beta(), degree), 4) << "   gamma = " << fixed(in_units(lat.gamma(), degree), 4) << " deg\n";
  out << "  volume " << general(lat.volume_A3(), 8) << " A^3   " << general(cell.atoms_per_cell(), 8)
      << " atoms per cell   " << general(in_units(cell.atomic_density(), per_cubic_angstrom) * 1e24, 6)
      << " atoms/cm3\n";
  out << "  density " << general(in_units(layer.density, gram_per_cubic_centimeter), 6)
      << " g/cm3 (standard atomic weights)" << (from_density ? "; lattice constant from the density" : "") << "\n";
  out << "  composition:";
  for (const auto &e : layer.elements) {
    out << "  " << Element::by_z(e.z).symbol() << " " << general(e.stoichiometry, 6) << " ("
        << fixed(100.0 * e.stoichiometry / cell.atoms_per_cell(), 3) << " at.%)";
  }
  out << "\n\nSites (fractional coordinates): " << cell.sites.size() << "\n";
  out << "  label     kind           x         y         z        occupants\n";
  const std::size_t shown = std::min<std::size_t>(cell.sites.size(), 60);
  for (std::size_t i = 0; i < shown; ++i) {
    const auto &s = cell.sites[i];
    char buf[256];
    std::snprintf(buf, sizeof buf, "  %-9s %-12s %9s %9s %9s   %s\n", s.label.c_str(), crystal::to_string(s.kind).c_str(),
                  fixed(s.position.x, 5).c_str(), fixed(s.position.y, 5).c_str(), fixed(s.position.z, 5).c_str(),
                  occupants_text(s).c_str());
    out << buf;
  }
  if (shown < cell.sites.size()) out << "  ... " << cell.sites.size() - shown << " more\n";

  // Neighbour shells of the first lattice site (how to recognise a structure).
  for (const auto &s : cell.sites) {
    if (s.kind != crystal::SiteKind::lattice) continue;
    const double reach = 1.8 * std::cbrt(lat.volume_A3() / std::max(1.0, cell.atoms_per_cell()));
    const auto shells = cell.neighbour_shells(s.position, reach);
    out << "\nLattice sites around " << s.label << ":";
    for (std::size_t k = 0; k < shells.size() && k < 4; ++k) {
      out << (k ? "," : "") << " " << general(shells[k].count, 4) << " at " << fixed(shells[k].distance_A, 4) << " A";
    }
    out << "\n";
    break;
  }

  const crystal::Mat3 rot = cr.rotation();
  const crystal::Mat3 to_crystal = rot.transposed();
  out << "\nOrientation: ";
  if (cr.orientation.kind == crystal::OrientationKind::euler) {
    out << "Bunge Euler angles " << general(cr.orientation.euler_deg[0]) << ", " << general(cr.orientation.euler_deg[1])
        << ", " << general(cr.orientation.euler_deg[2]) << " deg";
  } else {
    out << "surface " << crystal::format_indices(cr.orientation.surface, '(', ')');
    if (!cr.orientation.reference.empty()) {
      out << ", reference " << crystal::format_indices(cr.orientation.reference);
    }
  }
  out << ", twist " << general(in_units(cr.orientation.twist, degree)) << " deg\n";
  out << "  depth axis (lab x) = " << direction_text(lat, to_crystal * crystal::Vec3{1, 0, 0})
      << "   lab y = " << direction_text(lat, to_crystal * crystal::Vec3{0, 1, 0})
      << "   lab z = " << direction_text(lat, to_crystal * crystal::Vec3{0, 0, 1}) << "\n";
  if (!o.beam.empty()) {
    const auto uvw = parse_indices(o.beam, "--beam");
    const crystal::Vec3 dir = crystal::direction_vector(lat, crystal::direction_indices(lat, uvw, "--beam"));
    const auto tt = crystal::tilt_twist_for_beam(lat, cr.orientation, dir, "--beam");
    out << "Beam along " << crystal::format_indices(uvw) << ": ion.angle = " << general(in_units(tt.tilt, degree), 8)
        << " deg, orientation.twist = " << general(in_units(tt.twist, degree), 8) << " deg\n";
  }
  if (!o.tilt.empty()) {
    const angle_quantity tilt = config::parse_angle(o.tilt, "--tilt");
    if (!(in_units(tilt, degree) > -90.0 && in_units(tilt, degree) < 90.0)) {
      throw std::invalid_argument("--tilt: between -90 and 90 deg");
    }
    out << "Beam with ion.angle = " << general(in_units(tilt, degree)) << " deg travels along "
        << direction_text(lat, crystal::beam_in_crystal(lat, cr.orientation, tilt)) << "\n";
  }
  if (o.output) out << "\nWritten to " << *o.output << "\n";
  return 0;
}

} // namespace cli
} // namespace amcset
