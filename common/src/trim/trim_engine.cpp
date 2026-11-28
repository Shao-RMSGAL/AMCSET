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

/*! \file trim_engine.cpp
 *  \brief SRIM-mode transport: target precomputation, the flight/collision
 *  loop, damage bookkeeping and the parallel driver.
 */

#include "trim/trim.h"
#include "trim_export.h"
#include "target_model.h"
#include "../compute/device.h"

#include "data_files.h"
#include "diagnostics.h"

#include "amcset_common.h"
#include "config.h"
#include "crystal/grain_map.h"
#include "crystal/lattice_search.h"
#include "crystal/crystal_database.h"
#include "crystal/materials.h"
#include "element_data.h"
#include "physics/fast_math.h"
#include "trim/random.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <thread>
#include <unordered_set>

namespace amcset {
namespace trim {

using namespace common;
using physics::CollisionPair;

// ---------------------------------------------------------------------------
// Configuration helpers
// ---------------------------------------------------------------------------

std::string to_string(DamageCalculation d) {
  switch (d) {
  case DamageCalculation::kinchin_pease: return "kinchin-pease";
  case DamageCalculation::full_cascade: return "full-cascade";
  case DamageCalculation::monolayer: return "monolayer";
  }
  return "kinchin-pease";
}

std::uint32_t draw_random_seed() {
  // random_device may be a fixed sequence (some MinGW builds), so the clock
  // and a call counter are mixed in, then a splitmix64 finaliser spreads the
  // bits. 31 bits: the seed fits the GUI's field and TRIM.IN.
  static std::atomic<std::uint64_t> calls{0};
  std::uint64_t x = 0;
  try {
    std::random_device rd;
    x = (static_cast<std::uint64_t>(rd()) << 32) ^ static_cast<std::uint64_t>(rd());
  } catch (...) {
    // no entropy source: the clock alone
  }
  x ^= static_cast<std::uint64_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  x += 0x9E3779B97F4A7C15ull * (calls.fetch_add(1) + 1);
  x ^= x >> 30;
  x *= 0xBF58476D1CE4E5B9ull;
  x ^= x >> 27;
  x *= 0x94D049BB133111EBull;
  x ^= x >> 31;
  return static_cast<std::uint32_t>(x >> 33);
}

std::string to_string(TrimConfig::Transport t) {
  switch (t) {
  case TrimConfig::Transport::trim2013: return "trim2013";
  case TrimConfig::Transport::refined: return "refined";
  case TrimConfig::Transport::crystal: return "crystal";
  }
  return "trim2013";
}

TrimConfig::Transport parse_transport(std::string_view s) {
  std::string t;
  for (char c : s) t.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  if (t == "trim2013" || t == "trim") return TrimConfig::Transport::trim2013;
  if (t == "refined") return TrimConfig::Transport::refined;
  if (t == "crystal") return TrimConfig::Transport::crystal;
  throw std::invalid_argument("expected trim2013, refined or crystal, got '" + std::string(s) + "'");
}

std::string to_string(BeamSpot::Shape s) {
  switch (s) {
  case BeamSpot::Shape::point: return "point";
  case BeamSpot::Shape::rectangle: return "rectangle";
  case BeamSpot::Shape::circle: return "circle";
  case BeamSpot::Shape::gaussian: return "gaussian";
  }
  return "rectangle";
}

BeamSpot::Shape parse_beam_spot_shape(std::string_view s) {
  std::string t;
  for (char c : s) t.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  if (t == "point") return BeamSpot::Shape::point;
  if (t == "rectangle" || t == "uniform" || t == "square") return BeamSpot::Shape::rectangle;
  if (t == "circle" || t == "disc" || t == "disk") return BeamSpot::Shape::circle;
  if (t == "gaussian" || t == "gauss") return BeamSpot::Shape::gaussian;
  throw std::invalid_argument("expected point, rectangle, circle or gaussian, got '" + std::string(s) + "'");
}

DamageCalculation parse_damage_calculation(const std::string &s) {
  std::string t;
  for (char c : s) t.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
  if (t == "1" || t == "kp" || t == "quick" || t == "kinchin-pease" || t == "kinchin_pease") {
    return DamageCalculation::kinchin_pease;
  }
  if (t == "2" || t == "full" || t == "full-cascade" || t == "full_cascade" || t == "cascade") {
    return DamageCalculation::full_cascade;
  }
  if (t == "3" || t == "monolayer" || t == "sputtering" || t == "surface-sputtering") {
    return DamageCalculation::monolayer;
  }
  throw std::invalid_argument("Unknown damage calculation '" + s +
                              "' (expected kp, full or monolayer)");
}

TargetElement TargetElement::defaults(unsigned z, double stoichiometry) {
  const auto &el = Element::by_z(z);
  TargetElement e;
  e.z = z;
  e.mass = el.natural_mass();
  e.stoichiometry = stoichiometry;
  e.displacement_energy = default_displacement_energy;
  e.lattice_binding_energy = default_lattice_binding_energy;
  e.surface_binding_energy = el.sublimation_energy();
  e.name = std::string(el.symbol());
  return e;
}

IonBeam IonBeam::defaults(unsigned z, energy_quantity energy, angle_quantity angle) {
  IonBeam b;
  b.z = z;
  b.mass = Element::by_z(z).most_abundant_isotope_mass();
  b.energy = energy;
  b.angle = angle;
  return b;
}

std::vector<const TargetElement *> TrimConfig::trim_atoms() const {
  auto name_of = [](const TargetElement &e) {
    return e.name.empty() ? std::string(Element::by_z(e.z).symbol()) : e.name;
  };
  std::vector<const TargetElement *> atoms;
  for (const auto &l : layers) {
    for (const auto &e : l.elements) {
      bool known = false;
      for (const auto *a : atoms) known = known || (a->z == e.z && name_of(*a) == name_of(e));
      if (!known) atoms.push_back(&e);
    }
  }
  // Numbered atoms first, by number; the others keep their order.
  std::stable_sort(atoms.begin(), atoms.end(), [](const TargetElement *a, const TargetElement *b) {
    const unsigned ka = a->trim_atom == 0 ? std::numeric_limits<unsigned>::max() : a->trim_atom;
    const unsigned kb = b->trim_atom == 0 ? std::numeric_limits<unsigned>::max() : b->trim_atom;
    return ka < kb;
  });
  return atoms;
}

int TrimConfig::trim_code() const {
  const bool full = damage != DamageCalculation::kinchin_pease;
  switch (source) {
  case Source::particle_ions:
    return full ? 5 : 4;
  case Source::particle_recoils:
    // TRIM's code 6 runs full cascades from the atoms; code 7 runs the atoms
    // as ions with monolayer steps. (TRIM has no quick calculation for
    // recoils; AMCSET's is written as code 6.)
    return damage == DamageCalculation::monolayer ? 7 : 6;
  case Source::beam:
  default:
    return static_cast<int>(damage);
  }
}

std::size_t TrimConfig::particle_count() const {
  if (source == Source::beam) return ion_count;
  const std::size_t all = particles.size() * std::max<std::size_t>(particle_repeat, 1);
  return ion_count == 0 ? all : std::min(ion_count, all);
}

length_quantity TrimConfig::target_thickness() const {
  length_quantity t = 0.0 * angstrom;
  for (const auto &l : layers) t += l.width;
  return t;
}

std::vector<std::pair<const char *, bool TrimBugs::*>> TrimBugs::members() {
  return {{"stale_primary_energy", &TrimBugs::stale_primary_energy},
          {"recoil_energy_created", &TrimBugs::recoil_energy_created},
          {"struck_atom_direction", &TrimBugs::struck_atom_direction},
          {"replacement_interstitial", &TrimBugs::replacement_interstitial},
          {"subthreshold_sputter_counted", &TrimBugs::subthreshold_sputter_counted},
          {"ion_phonons_counted_twice", &TrimBugs::ion_phonons_counted_twice},
          {"cascade_random_replay", &TrimBugs::cascade_random_replay},
          {"zero_energy_defaults_skip_first_atom", &TrimBugs::zero_energy_defaults_skip_first_atom},
          {"integer_layer_z", &TrimBugs::integer_layer_z},
          {"exyz_skips_first_step", &TrimBugs::exyz_skips_first_step},
          {"recoil_event_counted_twice", &TrimBugs::recoil_event_counted_twice}};
}

void TrimBugs::set_all(bool reproduce) {
  for (const auto &m : members()) this->*(m.second) = reproduce;
}

bool TrimBugs::all_reproduced() const {
  for (const auto &m : members()) {
    if (!(this->*(m.second))) return false;
  }
  return true;
}

TrimBugs TrimConfig::effective_bugs() const {
  TrimBugs b = bugs;
  if (transport == Transport::crystal) b.set_all(false);
  return b;
}

bool TrimConfig::effective_single_precision_tables() const {
  return single_precision_tables && transport != Transport::crystal;
}

namespace {

//! "%g"-style text of a number for messages.
std::string num_text(double v) {
  char buf[48];
  std::snprintf(buf, sizeof buf, "%.6g", v);
  return buf;
}

} // namespace

void TrimConfig::validate(bool complete) const {
  // Messages name the field, the value, the allowed range and what to use
  // instead; the GUI shows them as they are.
  const auto fail = [](const std::string &m) { throw std::invalid_argument(m); };
  // Every value must be finite: "> 0" alone lets +inf through (an infinite
  // width or energy has no meaningful transport and gave NaN tallies).
  const auto finite = [](const auto &q) { return std::isfinite(static_cast<double>(q.value())); };
  const auto eV = [](energy_quantity e) { return num_text(in_units(e, electron_volt)) + " eV"; };
  const auto check_mass = [&](unsigned z, mass_quantity m, const std::string &where) {
    const double amu = (m / atomic_mass_unit).value();
    const auto [lo, hi] = mass_range_amu(z);
    if (!(amu >= lo && amu <= hi) || !std::isfinite(amu)) {
      fail(where + ": " + num_text(amu) + " amu is not a mass of " + std::string(Element::by_z(z).symbol()) +
           " (Z = " + std::to_string(z) + "); expected " + num_text(lo) + " to " + num_text(hi) +
           " amu, e.g. the natural " + num_text((Element::by_z(z).natural_mass() / atomic_mass_unit).value()) +
           " amu");
    }
  };

  // --- Ion --------------------------------------------------------------
  if (ion.z < 1 || ion.z > Element::max_z) {
    fail("ion.z: " + std::to_string(ion.z) + " is not an element AMCSET knows (1 to " +
         std::to_string(Element::max_z) + ")");
  }
  check_mass(ion.z, ion.mass, "ion.mass");
  if (!(ion.energy > 0.0 * electron_volt) || !finite(ion.energy)) {
    fail("ion.energy: must be positive and finite, got " + eV(ion.energy));
  }
  if (in_units(ion.energy, electron_volt) > max_ion_energy_eV) {
    fail("ion.energy: " + num_text(in_units(ion.energy, giga_electron_volt)) +
         " GeV is above the highest energy SRIM mode simulates (100 GeV: its stopping tables and "
         "Bethe-Bloch continuation end there); use at most 100 GeV");
  }
  {
    const energy_quantity stop = ion_stop_energy > 0.0 * electron_volt ? ion_stop_energy : final_energy;
    if (source == Source::beam && !(ion.energy > stop)) {
      fail("ion.energy: " + eV(ion.energy) + " is not above the energy at which ions stop (" + eV(stop) +
           "); give a higher ion energy or a lower final energy");
    }
  }
  if (!(ion.angle >= 0.0 * degree) || !(ion.angle < 90.0 * degree)) {
    fail("ion.angle: must be 0 to below 90 degrees (to the surface normal), got " +
         num_text(in_units(ion.angle, degree)) + " deg");
  }

  // --- Target -----------------------------------------------------------
  if (layers.empty()) fail("the target needs at least one layer");
  for (std::size_t i = 0; i < layers.size(); ++i) {
    const auto &l = layers[i];
    const std::string where = "layer " + std::to_string(i + 1) + " ('" + l.name + "')";
    if (!(l.width > 0.0 * angstrom) || !finite(l.width)) {
      fail(where + ": the thickness must be positive and finite, got " + num_text(in_units(l.width, angstrom)) +
           " A");
    }
    if (!(l.density > 0.0 * gram_per_cubic_centimeter) || !finite(l.density)) {
      fail(where + ": the density must be positive and finite, got " +
           num_text(in_units(l.density, gram_per_cubic_centimeter)) + " g/cm3");
    }
    if (in_units(l.density, gram_per_cubic_centimeter) > max_density_g_cm3) {
      fail(where + ": " + num_text(in_units(l.density, gram_per_cubic_centimeter)) +
           " g/cm3 is denser than any material (osmium: 22.6 g/cm3); at most " + num_text(max_density_g_cm3) +
           " g/cm3 is accepted");
    }
    if (!(l.compound_correction > 0.0) || !(l.compound_correction <= max_compound_correction)) {
      fail(where + ": the compound correction must be above 0 and at most " + num_text(max_compound_correction) +
           " (1 = Bragg's rule), got " + num_text(l.compound_correction));
    }
    if (l.elements.empty()) fail(where + ": needs at least one element");
    double sum = 0.0;
    mass_quantity mean_mass = 0.0 * kilogram;
    for (std::size_t k = 0; k < l.elements.size(); ++k) {
      const auto &e = l.elements[k];
      const std::string ew = where + ", element " + std::to_string(k + 1);
      if (e.z < 1 || e.z > Element::max_z) {
        fail(ew + ": Z = " + std::to_string(e.z) + " is not an element AMCSET knows (1 to " +
             std::to_string(Element::max_z) + ")");
      }
      check_mass(e.z, e.mass, ew + " mass");
      if (!(e.stoichiometry >= 0.0) || !std::isfinite(e.stoichiometry)) {
        fail(ew + ": the stoichiometry must be >= 0 and finite, got " + num_text(e.stoichiometry));
      }
      const auto check_energy = [&](energy_quantity v, double max_eV, const char *what) {
        if (!(v >= 0.0 * electron_volt) || !finite(v)) {
          fail(ew + ": " + what + " must be >= 0 and finite, got " + eV(v));
        }
        if (in_units(v, electron_volt) > max_eV) {
          fail(ew + ": " + what + " of " + eV(v) + " is not physical (at most " + num_text(max_eV) +
               " eV is accepted; TRIM's defaults are Ed 25 eV, Eb 3 eV, Es = heat of sublimation)");
        }
      };
      check_energy(e.displacement_energy, max_displacement_energy_eV, "the displacement energy Ed");
      if (const auto &de = e.damage_efficiency) {
        if (!(de->arc_b <= 0.0) || !std::isfinite(de->arc_b) || !(de->arc_c > 0.0 && de->arc_c <= 1.0) ||
            !(de->rpa_b > 0.0 * electron_volt) || !finite(de->rpa_b) || !(de->rpa_c > 0.0) ||
            !std::isfinite(de->rpa_c)) {
          fail(ew + ": damage_efficiency needs arc_b <= 0, 0 < arc_c <= 1, rpa_b > 0 and rpa_c > 0 "
                    "(Nordlund et al. 2018, Table 1)");
        }
      }
      check_energy(e.lattice_binding_energy, max_binding_energy_eV, "the lattice binding energy Eb");
      check_energy(e.surface_binding_energy, max_binding_energy_eV, "the surface binding energy Es");
      sum += e.stoichiometry;
    }
    if (!(sum > 0.0)) fail(where + ": the stoichiometries must sum to a positive value");
    for (const auto &e : l.elements) mean_mass += (e.stoichiometry / sum) * e.mass;
    check_layer_crystal(l, where); // not TRIM: the layer's optional crystal
    // A layer thinner than the distance between its atoms (N^-1/3) holds no
    // atomic layer: the transport needs at least one monolayer per layer.
    const number_density_quantity n = l.density / mean_mass;
    const double monolayer_A = in_units(1.0 / boost::units::root<3>(n), angstrom);
    if (in_units(l.width, angstrom) < monolayer_A) {
      fail(where + ": " + num_text(in_units(l.width, angstrom)) +
           " A is thinner than one atomic layer of this material (" + num_text(monolayer_A) +
           " A at this density); use at least " + num_text(std::ceil(monolayer_A * 100.0) / 100.0) + " A");
    }
  }

  // --- Particles, calculation ---------------------------------------------
  if (source == Source::beam) {
    if (ion_count == 0) fail("ions: the number of ions must be positive");
  } else {
    if (complete && particles.empty()) fail("the particle source (TRIM.DAT) has no particles");
    if (particle_repeat < 1 || particle_repeat > max_particle_repeat) {
      fail("source.repeat: 1 to " + std::to_string(max_particle_repeat) + " runs of every particle, got " +
           std::to_string(particle_repeat));
    }
    if (source == Source::particle_ions && damage == DamageCalculation::monolayer) {
      fail("TRIM.DAT ions run with the quick (Kinchin-Pease) or full-cascade damage calculation");
    }
    // (The refined transport starts every ion at the surface: it would run
    // the beam in place of the list.)
    if (transport == Transport::refined) {
      fail("physics.transport: the refined transport runs a beam only; particle sources (TRIM.DAT ions or "
           "recoils) run with trim2013 or crystal");
    }
  }
  if (threads > max_threads) {
    fail("threads: at most " + std::to_string(max_threads) + " (0 = all cores), got " + std::to_string(threads));
  }
  if (!(dpa_fluence_per_cm2 > 0.0) || !std::isfinite(dpa_fluence_per_cm2) || dpa_fluence_per_cm2 > 1e30) {
    fail("output.dpa_fluence: ions per cm2, above 0 and at most 1e30, got " + num_text(dpa_fluence_per_cm2));
  }
  {
    // The crystal transport and the beam spot (not TRIM).
    const double pm = in_units(crystal.max_impact_parameter, angstrom);
    if (!(pm >= 0.0) || !(pm <= 10.0)) {
      fail("physics.crystal.max_impact_parameter: 0 (automatic) to 10 A, got " + num_text(pm) + " A");
    }
    const double sd = in_units(crystal.simultaneous_distance, angstrom);
    if (!(sd >= 0.0) || !(sd <= 5.0)) {
      fail("physics.crystal.simultaneous_distance: 0 to 5 A, got " + num_text(sd) + " A");
    }
    if (!(crystal.local_stopping_fraction >= 0.0 && crystal.local_stopping_fraction <= 1.0)) {
      fail("physics.crystal.local_stopping_fraction: 0 to 1, got " + num_text(crystal.local_stopping_fraction));
    }
    if (!(crystal.local_stopping_decay >= 0.05 && crystal.local_stopping_decay <= 5.0)) {
      fail("physics.crystal.local_stopping_decay: 0.05 to 5 (0.3: the slowest Moliere exponential), got " +
           num_text(crystal.local_stopping_decay));
    }
    if (ion.spot.shape != BeamSpot::Shape::point) {
      for (const auto &[size, key] : {std::pair{ion.spot.width, "ion.spot.width"},
                                      std::pair{ion.spot.height, "ion.spot.height"}}) {
        const double a = in_units(size, angstrom);
        if (!(a > 0.0) || !(a <= 1e10)) {
          fail(std::string(key) + ": above 0 and at most 1 m, got " + num_text(a) + " A");
        }
      }
    }
  }
  {
    // Export files (not TRIM outputs).
    const auto check_file = [&](const std::string &f, const char *key) {
      const auto ext = std::filesystem::path(f).extension().string();
      if (ext != ".csv" && ext != ".parquet") {
        fail(std::string(key) + ": '" + f + "' must end in .csv or .parquet (the format)");
      }
    };
    if (!trajectory_export.file.empty()) {
      check_file(trajectory_export.file, "output.trajectory_file.file");
      if (trajectory_ions == 0) {
        fail("output.trajectory_file: no trajectories are recorded; set output.trajectories to a number of "
             "ions or \"all\"");
      }
      if (!(trajectory_export.min_energy >= 0.0 * electron_volt) || !finite(trajectory_export.min_energy)) {
        fail("output.trajectory_file.min_energy: must be >= 0 and finite, got " + eV(trajectory_export.min_energy));
      }
    }
    if (md_events.enabled()) {
      check_file(md_events.file, "output.md_events.file");
      if (damage == DamageCalculation::kinchin_pease) {
        fail("output.md_events: needs recoils that are followed: the full-cascade or monolayer damage "
             "calculation, not the quick (Kinchin-Pease) one");
      }
      const auto &ec = md_events.transition_energies;
      if (ec.empty() || ec.size() > 32) {
        fail("output.md_events.E_c: 1 to 32 transition energies, got " + std::to_string(ec.size()));
      }
      for (std::size_t i = 0; i < ec.size(); ++i) {
        if (!(ec[i] > 0.0 * electron_volt) || !finite(ec[i])) {
          fail("output.md_events.E_c: must be positive and finite, got " + eV(ec[i]));
        }
        if (i > 0 && !(ec[i] > ec[i - 1])) fail("output.md_events.E_c: give the values in increasing order, without repeats");
      }
      if (!(md_events.min_energy >= 0.0 * electron_volt) || !finite(md_events.min_energy)) {
        fail("output.md_events.E_min: must be >= 0 and finite (or \"Ed\"), got " + eV(md_events.min_energy));
      }
      if (md_events.min_energy >= ec.front()) {
        fail("output.md_events.E_min: " + eV(md_events.min_energy) + " must be below the smallest E_c (" +
             eV(ec.front()) + ")");
      }
      if (md_events.handoff && ec.size() != 1) {
        fail("output.md_events: the handoff mode takes exactly one E_c (a list of E_c needs mode \"record\")");
      }
    }
  }
  if (complete) {
    // Files the run will open (checked here so that `validate` and the GUI
    // report them before a run, not after the set-up).
    std::error_code ec;
    if (stopping_model == physics::ElectronicStoppingModel::srim_tables && !stopping_table.empty() &&
        !std::filesystem::exists(std::filesystem::path(stopping_table), ec)) {
      fail("stopping.table: '" + stopping_table +
           "' does not exist (leave it empty for the SRIM-2013 tables that come with AMCSET)");
    }
    if (!srim_data_directory.empty() && !std::filesystem::is_directory(std::filesystem::path(srim_data_directory), ec)) {
      fail("srim_data: '" + srim_data_directory +
           "' is not a directory (the Data folder of a SRIM-2013 installation, or empty)");
    }
  }
  if (!(auto_depth.margin >= 0.01 && auto_depth.margin <= 0.5)) {
    fail("auto_depth.margin: 0.01 to 0.5 (the deepest ion must stop in this fraction at the back), got " +
         num_text(auto_depth.margin));
  }
  if (auto_depth.test_ions < 1 || auto_depth.test_ions > 10000000) {
    fail("auto_depth.test_ions: 1 to 10000000 ions per test run, got " + std::to_string(auto_depth.test_ions));
  }
  if (!(exyz_interval >= 0.0 * electron_volt) || !finite(exyz_interval)) {
    fail("output.exyz_interval: must be >= 0 and finite, got " + eV(exyz_interval));
  }
  if (depth_bins == 0 || depth_bins > max_depth_bins) {
    fail("output.bins: 1 to " + std::to_string(max_depth_bins) + " depth bins (TRIM uses 100), got " +
         std::to_string(depth_bins));
  }
  const auto check_stop = [&](energy_quantity v, bool zero_ok, const char *what) {
    if (!(zero_ok ? v >= 0.0 * electron_volt : v > 0.0 * electron_volt) || !finite(v)) {
      fail(std::string(what) + (zero_ok ? ": must be >= 0" : ": must be positive") + " and finite, got " + eV(v));
    }
    if (in_units(v, electron_volt) > max_final_energy_eV) {
      fail(std::string(what) + ": " + eV(v) + " is too high (at most " + num_text(max_final_energy_eV) +
           " eV; TRIM uses 1 eV)");
    }
  };
  check_stop(final_energy, false, "physics.final_energy");
  check_stop(ion_stop_energy, true, "physics.ion_stop_energy");
  check_stop(recoil_stop_energy, true, "physics.recoil_stop_energy");
  check_stop(min_energy_transfer, false, "physics.min_energy_transfer");
  if (!(max_impact_parameter_ws >= 0.0) || !std::isfinite(max_impact_parameter_ws)) {
    fail("physics.max_impact_parameter_ws: must be >= 0 and finite");
  }
  if (!(surface_atom_depth >= 0.0 * angstrom) || !finite(surface_atom_depth)) {
    fail("physics.surface_atom_depth: must be >= 0 and finite");
  }
  const auto thickness = target_thickness();
  const auto end = window_end > 0.0 * angstrom ? window_end : thickness;
  if (window_start >= thickness) {
    fail("output.depth_window: starts at " + num_text(in_units(window_start, angstrom)) +
         " A, beyond the target (" + num_text(in_units(thickness, angstrom)) + " A thick)");
  }
  if (!(window_start >= 0.0 * angstrom) || !(end > window_start) || !finite(end)) {
    fail("output.depth_window: must satisfy 0 <= start < end, got " + num_text(in_units(window_start, angstrom)) +
         " to " + num_text(in_units(end, angstrom)) + " A");
  }
}

std::string TrimConfig::describe() const {
  std::ostringstream s;
  const auto &el = Element::by_z(ion.z);
  s << "SRIM mode: " << title << "\n";
  s << "  Ion: " << el.symbol() << " (Z=" << ion.z << ", M=" << (ion.mass / atomic_mass_unit).value()
    << " amu), E=" << in_units(ion.energy, kilo_electron_volt) << " keV, angle="
    << in_units(ion.angle, degree) << " deg\n";
  s << "  Damage: " << to_string(damage) << ", ions=" << ion_count << ", seed=" << seed << "\n";
  s << "  Stopping: " << physics::to_string(stopping_model);
  if (!stopping_table.empty()) s << " (" << stopping_table << ")";
  s << "; potential: ";
  if (srim_data_directory.empty() || pair_potentials == PairPotentialUse::none) {
    s << "ZBL universal";
  } else if (pair_potentials == PairPotentialUse::recoils) {
    s << "ZBL universal (ion), SRIM pair potentials (recoils)";
  } else {
    s << "SRIM pair potentials";
  }
  s << "\n";
  for (std::size_t i = 0; i < layers.size(); ++i) {
    const auto &l = layers[i];
    s << "  Layer " << i + 1 << " '" << l.name << "': " << in_units(l.width, angstrom)
      << " A, " << in_units(l.density, gram_per_cubic_centimeter) << " g/cm3"
      << (l.gas ? " (gas)" : "") << "\n";
    for (const auto &e : l.elements) {
      s << "    " << Element::by_z(e.z).symbol() << " x" << e.stoichiometry
        << "  M=" << (e.mass / atomic_mass_unit).value()
        << " Ed=" << in_units(e.displacement_energy, electron_volt)
        << " Eb=" << in_units(e.lattice_binding_energy, electron_volt)
        << " Es=" << in_units(e.surface_binding_energy, electron_volt) << " eV\n";
    }
  }
  return s.str();
}

// ---------------------------------------------------------------------------
// Results helpers
// ---------------------------------------------------------------------------

void DepthHistogram::add(const DepthHistogram &o) {
  if (values.size() != o.values.size()) {
    if (values.empty()) {
      values = o.values;
      return;
    }
    throw std::logic_error("DepthHistogram size mismatch");
  }
  for (std::size_t i = 0; i < values.size(); ++i) values[i] += o.values[i];
}

double Moments::stddev() const {
  if (n <= 0.0) return 0.0;
  const double m = mean();
  return std::sqrt(std::max(0.0, s2 / n - m * m));
}

double Moments::skewness() const {
  const double sd = stddev();
  if (n <= 0.0 || sd <= 0.0) return 0.0;
  const double m = mean();
  const double m3 = s3 / n - 3.0 * m * s2 / n + 2.0 * m * m * m;
  return m3 / (sd * sd * sd);
}

double Moments::kurtosis() const {
  const double sd = stddev();
  if (n <= 0.0 || sd <= 0.0) return 0.0;
  const double m = mean();
  const double m4 =
      s4 / n - 4.0 * m * s3 / n + 6.0 * m * m * s2 / n - 3.0 * m * m * m * m;
  return m4 / (sd * sd * sd * sd);
}

void TrimTallies::resize(std::size_t bins, std::size_t elements) {
  for (auto *h : {&ion_stopped, &knock_ons, &vacancies_by_ion, &vacancies_by_recoils,
                  &replacements, &displacements, &interstitials, &ionization_by_ion,
                  &ionization_by_recoils, &phonons_by_ion, &phonons_by_recoils,
                  &energy_to_recoils, &pka_damage_energy, &nrt_displacements, &kp_displacements,
                  &arc_displacements, &rpa_replacements}) {
    h->resize(bins);
  }
  lateral_bins.assign(bins, {0.0, 0.0, 0.0, 0.0, 0.0});
  recoil_stopped.assign(elements, DepthHistogram{});
  vacancies.assign(elements, DepthHistogram{});
  vacancies_removed.assign(elements, DepthHistogram{});
  energy_absorbed.assign(elements, DepthHistogram{});
  for (std::size_t e = 0; e < elements; ++e) {
    recoil_stopped[e].resize(bins);
    vacancies[e].resize(bins);
    vacancies_removed[e].resize(bins);
    energy_absorbed[e].resize(bins);
  }
  pka_spectrum.assign(spectrum_bins_per_decade * spectrum_decades, 0.0);
  recoil_spectrum.assign(spectrum_bins_per_decade * spectrum_decades, 0.0);
  sputtered.assign(elements, 0.0);
  sputtered_energy.assign(elements, 0.0);
  transmitted_recoils.assign(elements, 0.0);
  surface_energy.assign(elements, std::vector<double>(surface_energy_bins, 0.0));
}

TrimTallies TrimTallies::without_lists() const {
  TrimTallies t;
  t.ions = ions;
  t.ion_stopped = ion_stopped;
  t.range = range;
  t.lateral_y = lateral_y;
  t.lateral_abs_y = lateral_abs_y;
  t.radial = radial;
  t.lateral_bins = lateral_bins;
  t.recoil_stopped = recoil_stopped;
  t.knock_ons = knock_ons;
  t.vacancies = vacancies;
  t.vacancies_removed = vacancies_removed;
  t.vacancies_by_ion = vacancies_by_ion;
  t.vacancies_by_recoils = vacancies_by_recoils;
  t.replacements = replacements;
  t.displacements = displacements;
  t.interstitials = interstitials;
  t.ionization_by_ion = ionization_by_ion;
  t.ionization_by_recoils = ionization_by_recoils;
  t.phonons_by_ion = phonons_by_ion;
  t.phonons_by_recoils = phonons_by_recoils;
  t.energy_to_recoils = energy_to_recoils;
  t.energy_absorbed = energy_absorbed;
  t.pka_damage_energy = pka_damage_energy;
  t.nrt_displacements = nrt_displacements;
  t.kp_displacements = kp_displacements;
  t.arc_displacements = arc_displacements;
  t.rpa_replacements = rpa_replacements;
  t.total_pka_damage_energy = total_pka_damage_energy;
  t.total_nrt = total_nrt;
  t.total_kp = total_kp;
  t.total_arc = total_arc;
  t.total_rpa = total_rpa;
  t.pka_without_efficiency = pka_without_efficiency;
  t.md = md;
  t.grains = grains;
  t.total_displacements = total_displacements;
  t.total_vacancies = total_vacancies;
  t.total_replacements = total_replacements;
  t.total_interstitials = total_interstitials;
  t.total_ionization_ion = total_ionization_ion;
  t.total_ionization_recoils = total_ionization_recoils;
  t.total_phonons_ion = total_phonons_ion;
  t.total_phonons_recoils = total_phonons_recoils;
  t.total_energy_to_recoils = total_energy_to_recoils;
  t.total_binding_energy_loss = total_binding_energy_loss;
  t.backscattered = backscattered;
  t.backscattered_energy = backscattered_energy;
  t.transmitted = transmitted;
  t.transmitted_energy = transmitted_energy;
  t.sputtered = sputtered;
  t.sputtered_energy = sputtered_energy;
  t.transmitted_recoils = transmitted_recoils;
  t.recoils_followed = recoils_followed;
  t.recoil_limit_hits = recoil_limit_hits;
  t.pka_spectrum = pka_spectrum;
  t.recoil_spectrum = recoil_spectrum;
  t.surface_energy = surface_energy;
  return t;
}

void TrimTallies::add(const TrimTallies &o, std::size_t max_exit_records) {
  if (md.size() < o.md.size()) md.resize(o.md.size());
  for (std::size_t i = 0; i < o.md.size(); ++i) {
    md[i].transition_energy_eV = o.md[i].transition_energy_eV;
    md[i].add(o.md[i]);
  }
  for (const auto &[key, g] : o.grains) {
    auto &mine = grains[key];
    mine.entered += g.entered;
    mine.entered_stopped += g.entered_stopped;
    mine.depth_sum_A += g.depth_sum_A;
    mine.depth2_sum_A2 += g.depth2_sum_A2;
    mine.stopped_here += g.stopped_here;
  }
  // The per-element arrays are indexed by this tally's sizes below: they
  // must match (an empty tally takes the other's layout).
  if (recoil_stopped.empty() && vacancies.empty() && sputtered.empty() && !o.recoil_stopped.empty()) {
    resize(o.ion_stopped.values.size(), o.recoil_stopped.size());
  }
  const std::size_t ne = recoil_stopped.size();
  if (o.recoil_stopped.size() != ne || o.vacancies.size() != ne || vacancies.size() != ne ||
      o.energy_absorbed.size() != ne || energy_absorbed.size() != ne || o.sputtered.size() != ne ||
      sputtered.size() != ne || o.sputtered_energy.size() != ne || sputtered_energy.size() != ne ||
      o.transmitted_recoils.size() != ne || transmitted_recoils.size() != ne) {
    throw std::logic_error("TrimTallies::add: the tallies are for different targets (element count)");
  }
  ions += o.ions;
  ion_stopped.add(o.ion_stopped);
  range.add(o.range);
  lateral_y.add(o.lateral_y);
  lateral_abs_y.add(o.lateral_abs_y);
  radial.add(o.radial);
  for (std::size_t i = 0; i < lateral_bins.size() && i < o.lateral_bins.size(); ++i) {
    for (std::size_t k = 0; k < 5; ++k) lateral_bins[i][k] += o.lateral_bins[i][k];
  }
  for (std::size_t e = 0; e < recoil_stopped.size(); ++e) {
    recoil_stopped[e].add(o.recoil_stopped[e]);
    vacancies[e].add(o.vacancies[e]);
    if (e < vacancies_removed.size() && e < o.vacancies_removed.size()) vacancies_removed[e].add(o.vacancies_removed[e]);
    energy_absorbed[e].add(o.energy_absorbed[e]);
    sputtered[e] += o.sputtered[e];
    sputtered_energy[e] += o.sputtered_energy[e];
    transmitted_recoils[e] += o.transmitted_recoils[e];
    if (e < surface_energy.size() && e < o.surface_energy.size()) {
      for (std::size_t k = 0; k < surface_energy[e].size() && k < o.surface_energy[e].size(); ++k) {
        surface_energy[e][k] += o.surface_energy[e][k];
      }
    }
  }
  for (std::size_t i = 0; i < pka_spectrum.size() && i < o.pka_spectrum.size(); ++i) {
    pka_spectrum[i] += o.pka_spectrum[i];
    recoil_spectrum[i] += o.recoil_spectrum[i];
  }
  knock_ons.add(o.knock_ons);
  vacancies_by_ion.add(o.vacancies_by_ion);
  vacancies_by_recoils.add(o.vacancies_by_recoils);
  replacements.add(o.replacements);
  displacements.add(o.displacements);
  interstitials.add(o.interstitials);
  ionization_by_ion.add(o.ionization_by_ion);
  ionization_by_recoils.add(o.ionization_by_recoils);
  phonons_by_ion.add(o.phonons_by_ion);
  phonons_by_recoils.add(o.phonons_by_recoils);
  energy_to_recoils.add(o.energy_to_recoils);
  pka_damage_energy.add(o.pka_damage_energy);
  nrt_displacements.add(o.nrt_displacements);
  kp_displacements.add(o.kp_displacements);
  arc_displacements.add(o.arc_displacements);
  rpa_replacements.add(o.rpa_replacements);
  total_pka_damage_energy += o.total_pka_damage_energy;
  total_nrt += o.total_nrt;
  total_kp += o.total_kp;
  total_arc += o.total_arc;
  total_rpa += o.total_rpa;
  pka_without_efficiency += o.pka_without_efficiency;
  total_displacements += o.total_displacements;
  total_vacancies += o.total_vacancies;
  total_replacements += o.total_replacements;
  total_interstitials += o.total_interstitials;
  total_ionization_ion += o.total_ionization_ion;
  total_ionization_recoils += o.total_ionization_recoils;
  total_phonons_ion += o.total_phonons_ion;
  total_phonons_recoils += o.total_phonons_recoils;
  total_energy_to_recoils += o.total_energy_to_recoils;
  total_binding_energy_loss += o.total_binding_energy_loss;
  backscattered += o.backscattered;
  backscattered_energy += o.backscattered_energy;
  transmitted += o.transmitted;
  transmitted_energy += o.transmitted_energy;
  recoils_followed += o.recoils_followed;
  recoil_limit_hits += o.recoil_limit_hits;
  for (const auto &r : o.exits) {
    if (exits.size() >= max_exit_records) break;
    exits.push_back(r);
  }
  positions.insert(positions.end(), o.positions.begin(), o.positions.end());
  cascades.insert(cascades.end(), o.cascades.begin(), o.cascades.end());
  ion_summaries.insert(ion_summaries.end(), o.ion_summaries.begin(), o.ion_summaries.end());
  for (const auto &r : o.exyz) {
    if (exyz.size() >= max_exit_records) break;
    exyz.push_back(r);
  }
}

// ---------------------------------------------------------------------------
// Precomputed target model
// ---------------------------------------------------------------------------

namespace {

// Robinson partition with precomputed k_d and eps_d / T[eV] (see
// robinson_damage_energy()); energies in eV.
double robinson_damage(double kd, double eps_per_eV, double T_eV) {
  const double eps = eps_per_eV * T_eV;
  const double g = eps + 0.40244 * std::pow(eps, 0.75) + 3.4008 * std::pow(eps, 1.0 / 6.0);
  return T_eV / (1.0 + kd * g);
}

double robinson_kd(unsigned z, mass_quantity mass) {
  return 0.1334 * std::pow(static_cast<double>(z), 2.0 / 3.0) /
         std::sqrt((mass / atomic_mass_unit).value());
}

double robinson_eps_per_eV(unsigned z) {
  return 0.01014 * std::pow(static_cast<double>(z), -7.0 / 3.0);
}

} // namespace

energy_quantity robinson_damage_energy(unsigned z, mass_quantity mass,
                                       energy_quantity recoil_energy) {
  if (recoil_energy <= 0.0 * electron_volt) return 0.0 * electron_volt;
  return robinson_damage(robinson_kd(z, mass), robinson_eps_per_eV(z),
                         in_units(recoil_energy, electron_volt)) *
         electron_volt;
}

double kinchin_pease_vacancies(unsigned z, mass_quantity mass,
                               energy_quantity displacement_energy,
                               energy_quantity recoil_energy) {
  const double ev = in_units(robinson_damage_energy(z, mass, recoil_energy), electron_volt);
  const double ed = in_units(displacement_energy, electron_volt);
  if (ev < ed) return 0.0;
  if (ev > 2.5 * ed) return 0.8 * ev / (2.0 * ed);
  return 1.0;
}


namespace {

std::size_t species_index(std::vector<SpeciesInfo> &species, unsigned z,
                          mass_quantity m) {
  for (std::size_t i = 0; i < species.size(); ++i) {
    if (species[i].z == z &&
        std::abs((species[i].mass / m).value() - 1.0) < 1e-9) {
      return i;
    }
  }
  species.push_back({z, m});
  return species.size() - 1;
}

constexpr double ln10 = 2.302585092994046;

} // namespace

// ---------------------------------------------------------------------------
// TrimSimulation construction
// ---------------------------------------------------------------------------

TrimSimulation::TrimSimulation(TrimConfig config)
    : config_(std::move(config)), model_(std::make_unique<TargetModel>()) {
  const auto t_setup = std::chrono::steady_clock::now();
  // physics.crystal.auto_crystal: layers without a crystal whose composition
  // is a database crystal become that crystal (the crystal transport only).
  if (config_.transport == TrimConfig::Transport::crystal && config_.crystal.auto_crystal) {
    for (auto &layer : config_.layers) {
      if (layer.crystal || layer.gas) continue;
      std::vector<std::pair<unsigned, double>> composition;
      for (const auto &e : layer.elements) composition.emplace_back(e.z, e.stoichiometry);
      const auto *row = crystal::match_composition(composition);
      if (row == nullptr) continue;
      const double before = in_units(layer.density, gram_per_cubic_centimeter);
      layer.crystal = crystal::make_database_crystal(*row);
      apply_crystal(layer);
      char buf[64];
      std::snprintf(buf, sizeof buf, "%.5g A, density %.4g g/cm3 (was %.4g)", row->a_A,
                    in_units(layer.density, gram_per_cubic_centimeter), before);
      warnings_.push_back("layer '" + layer.name + "': followed as the crystal " + std::string(row->name) + " (" +
                          std::string(row->structure) + ", a = " + buf + ") of the crystal database, (0 0 1) " +
                          "surface, 293 K (physics.crystal.auto_crystal)");
    }
  }
  config_.validate();
  // A random seed is drawn once, here: the run is then deterministic, and
  // config() and the outputs carry the seed it used (random_seed stays set
  // so that the outputs say the seed was drawn).
  if (config_.random_seed) config_.seed = draw_random_seed();
  // The crystal transport is not TRIM: it takes none of TRIM.exe's bugs,
  // whatever the switches say (TrimConfig::effective_bugs).
  config_.bugs = config_.effective_bugs();
  config_.single_precision_tables = config_.effective_single_precision_tables();
  auto &m = *model_;

  // TRIM.exe (0x44efc0) replaces Ed <= 0 by 20 eV, Es <= 0 by 4 eV and
  // Eb <= 0 by 2 eV for its target atoms 2..n; its loop starts at the second
  // atom, so the first keeps zeros. TDATA.txt shows the replaced values.
  layers_ = config_.layers;
  if (config_.transport == TrimConfig::Transport::trim2013) {
    const auto atoms = config_.trim_atoms();
    auto name_of = [](const TargetElement &e) {
      return e.name.empty() ? std::string(Element::by_z(e.z).symbol()) : e.name;
    };
    for (std::size_t k = 0; k < atoms.size(); ++k) {
      const unsigned position = atoms[k]->trim_atom > 0 ? atoms[k]->trim_atom : static_cast<unsigned>(k + 1);
      if (position < 2 && config_.bugs.zero_energy_defaults_skip_first_atom) continue;
      std::vector<std::string> changed;
      for (auto &l : layers_) {
        for (auto &e : l.elements) {
          if (e.z != atoms[k]->z || name_of(e) != name_of(*atoms[k])) continue;
          auto fix = [&](energy_quantity &v, double value, const char *what) {
            if (v > 0.0 * electron_volt) return;
            v = value * electron_volt;
            const std::string text = std::string(what) + " " + std::to_string(static_cast<int>(value)) + " eV";
            if (std::find(changed.begin(), changed.end(), text) == changed.end()) changed.push_back(text);
          };
          fix(e.displacement_energy, 20.0, "displacement energy");
          fix(e.surface_binding_energy, 4.0, "surface binding energy");
          fix(e.lattice_binding_energy, 2.0, "lattice binding energy");
        }
      }
      if (!changed.empty()) {
        std::string list;
        for (std::size_t i = 0; i < changed.size(); ++i) list += (i ? ", " : "") + changed[i];
        warnings_.push_back("Target atom " + std::to_string(position) + " (" + name_of(*atoms[k]) +
                            (config_.bugs.zero_energy_defaults_skip_first_atom
                                 ? "): a zero energy is replaced as TRIM does for every target atom after "
                                   "the first: "
                                 : "): a zero energy is replaced by TRIM's default: ") +
                            list + ".");
      }
    }
  }

  if (config_.stopping_model == physics::ElectronicStoppingModel::srim_tables) {
    if (config_.stopping_table.empty()) {
      // The SRIM-2013 tables that come with AMCSET: TRIM.exe computes its
      // stopping from the same coefficients (Data/SCOEF03.DAT).
      config_.stopping_table = data::bundled_srim_stopping_table();
      if (config_.stopping_table.empty()) {
        warnings_.push_back("The SRIM-2013 stopping tables that come with AMCSET (data/" +
                            std::string(data::srim_stopping_table) +
                            ") were not found: the built-in stopping is used instead, which differs "
                            "from TRIM's by up to 25 %. Set the stopping table or AMCSET_DATA_DIR.");
        config_.stopping_model = physics::ElectronicStoppingModel::builtin;
      }
    }
    if (config_.stopping_model == physics::ElectronicStoppingModel::srim_tables) {
      m.stopping_db = physics::SrimStoppingDatabase::load(config_.stopping_table);
    }
  }
  if (!config_.srim_data_directory.empty()) {
    m.potentials = physics::PairPotentialTable::load_srim(config_.srim_data_directory);
  }

  // Species: the ion first, then every distinct target atom.
  m.species.push_back({config_.ion.z, config_.ion.mass});

  // Tally element index: one entry per (layer, element) as TRIM numbers its
  // "target atoms" (elements repeated in several layers get one column each).
  length_quantity depth = 0.0 * angstrom;
  m.lowest_ed = std::numeric_limits<double>::max() * electron_volt;
  // Tally columns in TRIM's atom order.
  for (const auto *a : config_.trim_atoms()) {
    element_names_.push_back(a->name.empty() ? std::string(Element::by_z(a->z).symbol()) : a->name);
    element_z_.push_back(a->z);
  }
  for (const auto &l : layers_) {
    LayerModel lm;
    lm.top = depth;
    depth += l.width;
    lm.bottom = depth;
    lm.gas = l.gas;
    double total = 0.0;
    for (const auto &e : l.elements) total += e.stoichiometry;
    mass_quantity mean_mass = 0.0 * kilogram;
    double cumulative = 0.0;
    for (const auto &e : l.elements) {
      ElementModel em;
      em.z = e.z;
      em.mass = e.mass;
      em.fraction = e.stoichiometry / total;
      cumulative += em.fraction;
      em.cumulative = cumulative;
      em.ed = e.displacement_energy;
      em.eb = e.lattice_binding_energy;
      em.es = e.surface_binding_energy;
      // Elements with the same name and Z share one tally column (as one
      // TRIM "target atom" that appears in several layers).
      const std::string name =
          e.name.empty() ? std::string(Element::by_z(e.z).symbol()) : e.name;
      em.tally_index = element_names_.size();
      for (std::size_t k = 0; k < element_names_.size(); ++k) {
        if (element_names_[k] == name && element_z_[k] == e.z) em.tally_index = k;
      }
      if (em.tally_index == element_names_.size()) {
        element_names_.push_back(name);
        element_z_.push_back(e.z);
      }
      em.species = species_index(m.species, e.z, e.mass);
      // Robinson's partition, as used by TRIM (SRIM book Eqs. 7-50..7-53):
      // k_d = 0.1334 Z^(2/3) M^(-1/2) (M in amu) and
      // eps_d = 0.01014 Z^(-7/3) T (T in eV).
      em.kp_kd = robinson_kd(e.z, e.mass);
      em.kp_eps_per_eV = robinson_eps_per_eV(e.z);
      if (const auto de = e.damage_efficiency ? e.damage_efficiency : published_damage_efficiency(e.z)) {
        em.has_efficiency = true;
        em.arc_b = de->arc_b;
        em.arc_c = de->arc_c;
        em.rpa_c = de->rpa_c;
        const double knee = 2.0 * in_units(em.ed, electron_volt) / 0.8;
        em.arc_scale = (1.0 - em.arc_c) / std::pow(knee, em.arc_b);
        em.rpa_bc = std::pow(in_units(de->rpa_b, electron_volt), em.rpa_c);
        em.rpa_scale = em.rpa_bc / std::pow(knee, em.rpa_c) + 1.0;
      }
      mean_mass += em.fraction * e.mass;
      if (e.displacement_energy > 0.0 * electron_volt) {
        m.lowest_ed = std::min(m.lowest_ed, e.displacement_energy);
      }
      lm.elements.push_back(em);
    }
    lm.elements.back().cumulative = 1.0;
    lm.density = l.density / mean_mass;
    lm.monolayer = 1.0 / boost::units::root<3>(lm.density);
    if (config_.transport == TrimConfig::Transport::crystal && l.crystal) {
      // The lattice the crystal transport follows: thermal amplitudes per
      // element and the largest impact parameter (see CrystalPhysics).
      std::map<unsigned, double> rms_A;
      if (config_.crystal.thermal_vibrations) {
        for (const auto &e : l.elements) {
          if (const auto u = crystal::thermal_amplitude(*l.crystal, e.z, e.mass)) rms_A[e.z] = in_units(*u, angstrom);
        }
      }
      double pmax_A = in_units(config_.crystal.max_impact_parameter, angstrom);
      if (!(pmax_A > 0.0)) pmax_A = crystal::largest_hole_radius_A(l.crystal->cell);
      lm.lattice = std::make_shared<crystal::LatticeSearch>(l.crystal->cell, l.crystal->rotation(), rms_A, pmax_A,
                                                            pmax_A);
      if (l.crystal->polycrystal) {
        // A polycrystal: the grains of the layer, each with its own lattice
        // orientation (the site data is the layer's, shared).
        lm.grains = crystal::make_grain_map(*l.crystal, in_units(l.width, angstrom));
        lm.boundary_half_A = 0.5 * in_units(l.crystal->polycrystal->boundary_width, angstrom);
        lm.grain_lattices.reserve(lm.grains->grain_count());
        for (const auto &g : lm.grains->grains()) {
          lm.grain_lattices.emplace_back(*lm.lattice, l.crystal->cell, g.rotation);
        }
      }
    }
    lm.domain = m.layers.empty() ? 0u
                                 : m.layers.back().domain +
                                       static_cast<std::uint32_t>(std::max<std::size_t>(m.layers.back().grain_lattices.size(), 1));
    {
      // TRIM-2013 layer averages (TRIMSTART). Constants as in TRIM.exe:
      // a0 = 0.5292 A, 0.8853, e^2 = 14.4 eV A, E_min = 5 eV, Bohr
      // straggling 177 = sqrt(4 pi) e^2 sqrt(12).
      auto &t = lm.t13;
      double zacc = 0.0;
      for (const auto &em : lm.elements) {
        if (config_.bugs.integer_layer_z) {
          zacc = std::nearbyint(em.fraction * static_cast<double>(em.z) + zacc); // VB CInt
        } else {
          zacc += em.fraction * static_cast<double>(em.z);
        }
      }
      t.z_int = std::max(zacc, 1.0);
      t.mass = mean_mass;
      t.pmax0 = lm.monolayer / std::sqrt(pi_v);
      const double z1 = static_cast<double>(config_.ion.z);
      const double m1 = (config_.ion.mass / atomic_mass_unit).value();
      const double m2 = (mean_mass / atomic_mass_unit).value();
      t.screening = 0.5292 * 0.8853 / (std::pow(z1, 0.23) + std::pow(t.z_int, 0.23)) * angstrom;
      const auto e2 = 14.4 * electron_volt * angstrom;
      t.fi = t.screening * m2 / (z1 * t.z_int * e2 * (m1 + m2));
      const double mu = m1 / m2;
      t.epsdg = (5.0 * electron_volt * t.fi).value() * (1.0 + mu) * (1.0 + mu) / (4.0 * mu);
      t.kp_eps_per_eV = 0.01 * std::pow(t.z_int, -7.0 / 3.0);
      t.kp_k = 0.1334 * std::pow(t.z_int, 2.0 / 3.0) / std::sqrt(m2);
      const double n_per_A3 = in_units(lm.density, per_cubic_angstrom);
      t.k_low = std::sqrt(3.2 * 12.0 * n_per_A3 / (1000.0 * m1) * std::pow(z1 + t.z_int, 8.0 / 3.0));
      t.k_bohr = 177.0 * z1 * std::sqrt(t.z_int * n_per_A3);
    }
    if (config_.max_impact_parameter_ws > 0.0) {
      // Wigner-Seitz radius r_ws = (3 / (4 pi N))^(1/3).
      const length_quantity r_ws =
          boost::units::root<3>(3.0 / (4.0 * pi_v * lm.density));
      const length_quantity p_cap = config_.max_impact_parameter_ws * r_ws;
      lm.min_flight = 1.0 / (pi_v * p_cap * p_cap * lm.density);
    } else {
      lm.min_flight = lm.monolayer;
    }
    m.layers.push_back(std::move(lm));
  }
  if (m.lowest_ed == std::numeric_limits<double>::max() * electron_volt) {
    m.lowest_ed = 1.0 * electron_volt;
  }
  m.thickness = depth;

  // Per (species, layer) stopping, straggling and collision pairs.
  m.data.resize(m.species.size());
  for (std::size_t s = 0; s < m.species.size(); ++s) {
    const auto &sp = m.species[s];
    for (std::size_t li = 0; li < m.layers.size(); ++li) {
      const auto &lm = m.layers[li];
      const auto &lc = config_.layers[li];
      SpeciesLayer d;
      std::vector<physics::StoppingComponent> comps;
      for (const auto &e : lm.elements) comps.push_back({e.z, e.fraction});
      d.stopping = physics::StoppingCurve(config_.stopping_model, m.stopping_db.get(), sp.z,
                                          sp.mass, comps, lc.gas, lc.compound_correction,
                                          lm.density, config_.bethe_bloch, 1.0 * electron_volt,
                                          100.0 * giga_electron_volt, 50, config_.stopping_corrections);
      if (config_.fast_math) d.stopping.enable_fast_lookup(); // no log/exp per lookup
      for (const auto &e : lm.elements) {
        // The ion (species 0) and the recoils may use different potentials.
        const bool use_pair =
            config_.pair_potentials == TrimConfig::PairPotentialUse::all ||
            (config_.pair_potentials == TrimConfig::PairPotentialUse::recoils && s != 0);
        d.pairs.emplace_back(sp.z, sp.mass, e.z, e.mass,
                             use_pair ? m.potentials.get() : nullptr);
        d.eps_min.push_back(d.pairs.back().reduced_energy(config_.min_energy_transfer));
      }
      // Straggling table (1 eV .. 100 GeV, 20 points per decade).
      constexpr unsigned ppd = 20;
      d.strag_log_e0 = 0.0;
      d.strag_inv_step = ppd / ln10;
      const std::size_t n = 11 * ppd + 1;
      d.strag.resize(n);
      for (std::size_t i = 0; i < n; ++i) {
        const energy_quantity e = std::exp(static_cast<double>(i) / d.strag_inv_step) * electron_volt;
        double v = 0.0;
        for (const auto &el : lm.elements) {
          v += el.fraction *
               in_units(physics::trim_electronic_straggling(sp.z, sp.mass, el.z, e, lm.density),
                        electron_volt * electron_volt / angstrom);
        }
        d.strag[i] = v;
      }
      m.data[s].push_back(std::move(d));
    }
  }

  // TRIM.DAT particles: TRIM skips (with a message) lines whose energy exceeds
  // the TRIM.IN energy (the stopping tables end there) or whose depth is
  // outside the target, and for recoils atoms that are not target elements.
  if (config_.source != TrimConfig::Source::beam) {
    const bool recoils = config_.source == TrimConfig::Source::particle_recoils;
    const double emax = in_units(config_.ion.energy, electron_volt);
    std::size_t line = 0;
    for (const auto &p : config_.particles) {
      ++line;
      const auto skip = [&](const std::string &why) {
        warnings_.push_back("particle " + std::to_string(line) + " ('" + p.name + "'): " + why +
                            "; skipped");
      };
      const double e = in_units(p.energy, electron_volt);
      if (!(e >= 0.0)) {
        skip("negative energy");
        continue;
      }
      if (e > emax * (1.0 + 1e-12)) {
        skip("energy " + std::to_string(e) + " eV above the maximum (ion energy) " +
             std::to_string(emax) + " eV");
        continue;
      }
      if (!(p.x >= 0.0 * angstrom) || p.x > m.thickness) {
        skip("depth " + std::to_string(in_units(p.x, angstrom)) + " A outside the target");
        continue;
      }
      if (!std::isfinite(in_units(p.y, angstrom)) || !std::isfinite(in_units(p.z_position, angstrom))) {
        skip("lateral position is not a finite number");
        continue;
      }
      {
        // A NaN or zero direction cannot be followed (NaN would never stop).
        const double n2 = p.cos_x * p.cos_x + p.cos_y * p.cos_y + p.cos_z * p.cos_z;
        if (!(n2 > 0.0) || !std::isfinite(n2)) {
          skip("direction is zero or not finite");
          continue;
        }
      }
      ResolvedParticle rp;
      rp.particle = p;
      rp.layer = m.layer_index(p.x);
      if (recoils) {
        // The atom is the first target element with this Z (TRIM searches its
        // element list); its constants come from the layer at the start
        // depth, or from the first layer that contains it.
        bool found = false;
        for (std::size_t pass = 0; pass < 2 && !found; ++pass) {
          for (std::size_t li = 0; li < m.layers.size() && !found; ++li) {
            if (pass == 0 && li != rp.layer) continue;
            const auto &els = m.layers[li].elements;
            for (std::size_t j = 0; j < els.size(); ++j) {
              if (els[j].z == p.z) {
                rp.element_layer = li;
                rp.element = j;
                rp.species = els[j].species;
                found = true;
                break;
              }
            }
          }
        }
        if (!found) {
          skip("Z = " + std::to_string(p.z) + " is not a target element");
          continue;
        }
        if (rp.element_layer != rp.layer) {
          warnings_.push_back("particle " + std::to_string(line) + " ('" + p.name + "'): Z = " +
                              std::to_string(p.z) + " is not an element of layer " +
                              std::to_string(rp.layer + 1) + " where it starts");
        }
      } else if (p.z != config_.ion.z) {
        warnings_.push_back("particle " + std::to_string(line) + " ('" + p.name + "'): Z = " +
                            std::to_string(p.z) + " differs from the ion; run as Z = " +
                            std::to_string(config_.ion.z) + " (as TRIM does)");
      }
      m.sources.push_back(std::move(rp));
    }
    if (m.sources.empty()) {
      throw std::invalid_argument("no usable particle in the particle source (TRIM.DAT)");
    }
    // Every particle runs source.repeat times, one after the other (run i is
    // particle i / repeat, with its own random-number stream).
    const std::size_t runs = m.sources.size() * config_.particle_repeat;
    total_ = config_.ion_count == 0 ? runs : std::min(config_.ion_count, runs);
    if (config_.first_ion >= runs) {
      throw std::invalid_argument("first_ion is beyond the end of the particle list");
    }
    total_ = std::min(total_, runs - config_.first_ion);
  } else {
    total_ = config_.ion_count;
  }

  tallies_.resize(config_.depth_bins, element_names_.size());
  if (config_.fast_math) {
    // The tables of the fast approximations are built now, not in the first collision.
    physics::MagicTable::universal();
    physics::fast::kp_partition(1.0);
    physics::fast::flight_path_denominator(1.0);
    physics::fast::log10_at_least_one(2.0);
  }
  AMCSET_LOG(info) << "SRIM-mode setup: " << m.species.size() << " species, " << m.layers.size()
                   << " layers, stopping " << physics::to_string(config_.stopping_model)
                   << (config_.stopping_table.empty() ? "" : " (" + config_.stopping_table + ")")
                   << ", " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t_setup).count()
                   << " s";
  for (const auto &w : config_.plausibility_warnings()) warnings_.push_back(w);
  if (config_.compute.backend != compute::Backend::cpu) {
    // GPU backends: a run they cannot take goes to the CPU (with a warning).
    const std::string why = compute::unsupported_reason(config_);
    if (!why.empty()) {
      warnings_.push_back("compute." + compute::to_string(config_.compute.backend) + ": " + why +
                          "; this calculation runs on the CPU");
    } else {
      for (const auto &n : compute::device_notes(config_)) warnings_.push_back(n);
    }
  }
  for (const auto &w : warnings_) AMCSET_LOG(warning) << "SRIM-mode input: " << w;
}

TrimSimulation::ComputeReport TrimSimulation::compute_report() const {
  std::lock_guard<std::mutex> lock(compute_mutex_);
  return compute_report_;
}

TrimTallies TrimSimulation::snapshot(bool with_lists) const {
  std::lock_guard<std::mutex> lock(tallies_mutex_);
  return with_lists ? tallies_ : tallies_.without_lists();
}

void TrimSimulation::continue_from(const TrimTallies &previous) {
  const std::size_t ne = element_names_.size();
  if (previous.ion_stopped.values.size() != config_.depth_bins || previous.vacancies.size() != ne ||
      previous.recoil_stopped.size() != ne || previous.energy_absorbed.size() != ne ||
      previous.sputtered.size() != ne || previous.sputtered_energy.size() != ne ||
      previous.transmitted_recoils.size() != ne ||
      previous.lateral_bins.size() != config_.depth_bins) {
    throw std::invalid_argument(
        "the tallies to continue from do not match this target (depth bins or elements)");
  }
  std::lock_guard<std::mutex> lock(tallies_mutex_);
  tallies_ = previous;
  transmitted_.store(tallies_.transmitted, std::memory_order_relaxed);
  if (tallies_.surface_energy.size() != element_names_.size()) {
    tallies_.surface_energy.assign(element_names_.size(),
                                   std::vector<double>(TrimTallies::surface_energy_bins, 0.0));
  }
  if (tallies_.vacancies_removed.size() != element_names_.size()) {
    tallies_.vacancies_removed.assign(element_names_.size(), DepthHistogram{});
    for (auto &h : tallies_.vacancies_removed) h.resize(config_.depth_bins);
  }
}

void TrimSimulation::request_cancel() noexcept {
  cancel_.store(true);
  try {
    std::lock_guard<std::mutex> lock(pause_mutex_);
    paused_.store(false);
  } catch (...) {
    paused_.store(false);
  }
  pause_cv_.notify_all();
}

void TrimSimulation::set_paused(bool paused) {
  {
    std::lock_guard<std::mutex> lock(pause_mutex_);
    paused_.store(paused);
  }
  pause_cv_.notify_all();
}

std::optional<DamageEfficiency> published_damage_efficiency(unsigned z) {
  // Nordlund et al., Nat. Commun. 9, 1084 (2018), Table 1 (fits to MD; the
  // Ed each was fitted with).
  struct Row {
    unsigned z;
    double ed, arc_b, arc_c, rpa_b, rpa_c;
  };
  static constexpr Row rows[] = {
      {26, 40.0, -0.568, 0.286, 1018.0, 0.95},  // Fe
      {29, 33.0, -0.68, 0.16, 3319.0, 0.97},    // Cu
      {28, 39.0, -1.01, 0.23, 3325.0, 0.92},    // Ni
      {46, 41.0, -0.88, 0.15, 2065.0, 1.08},    // Pd
      {78, 42.0, -1.12, 0.11, 5531.0, 0.87},    // Pt
      {74, 70.0, -0.56, 0.12, 12332.0, 0.73},   // W
  };
  for (const auto &r : rows) {
    if (r.z == z) {
      return DamageEfficiency{r.arc_b, r.arc_c, r.rpa_b * electron_volt, r.rpa_c, r.ed * electron_volt};
    }
  }
  return std::nullopt;
}

TrimSimulation::~TrimSimulation() = default;

length_quantity TrimSimulation::depth_extent() const { return model_->thickness; }

std::shared_ptr<const crystal::GrainMap> TrimSimulation::grain_map(std::size_t layer) const {
  return layer < model_->layers.size() ? model_->layers[layer].grains : nullptr;
}

energy_quantity TrimSimulation::lowest_displacement_energy() const { return model_->lowest_ed; }

double TrimSimulation::bin_width() const {
  const auto end = config_.window_end > 0.0 * angstrom ? config_.window_end : model_->thickness;
  return in_units(end - config_.window_start, angstrom) / static_cast<double>(config_.depth_bins);
}

double TrimSimulation::bin_start(std::size_t i) const {
  return in_units(config_.window_start, angstrom) + bin_width() * static_cast<double>(i);
}

// ---------------------------------------------------------------------------
// Transport
// ---------------------------------------------------------------------------

namespace {

//! A moving atom (the ion or a recoil). Positions in Angstrom, energies in eV
//! are kept as Boost.Units quantities.
struct Mover {
  std::size_t species = 0;
  unsigned z = 0;
  length_quantity x, y, zz; // TRIM axes: x = depth
  double cx = 1.0, cy = 0.0, cz = 0.0;
  energy_quantity energy;
  bool is_ion = false;
  unsigned generation = 0;
  std::size_t tally_element = 0; // for recoils: tally index of its element
  bool surface_only = false;      // set in motion below Ed: no defect if it stays
  double origin_x = 0.0;          // depth where the atom was set in motion (Angstrom)
  std::size_t layer_hint = 0;     // trim2013: layer of the atom's site
  std::size_t element_hint = 0;   // trim2013: element index within origin_layer
  // trim2013: the layer whose element entry (element_hint) describes the
  // atom. layer_hint follows the atom from layer to layer; the element index
  // belongs to the layer the atom was struck in.
  std::size_t origin_layer = 0;
  std::uint32_t track = 0;        // path number within the ion's trajectory
  // Trajectory lineage and clock (trajectories, MD events; no effect on the
  // transport): the collision that set the atom in motion, and the time
  // since the ion entered the target (advanced only when recorded).
  std::uint32_t parent_track = CollisionPoint::no_parent;
  std::uint32_t parent_point = CollisionPoint::no_parent;
  common::time_quantity time{};
  std::uint32_t md_mask = 0;      // MD events: bit i = the branch is no longer the BCA's for E_c[i]
  bool handed = false;            // MD handoff: MD takes this atom; the BCA does not follow it
  bool untracked = false;         // never moved (TRIM's stale-energy primary): no trajectory points
};

/*!
 * Gamma(k, 1) / k deviate (mean 1, variance 1/k) by Marsaglia & Tsang's
 * method (ACM TOMS 26, 363 (2000)), with the k < 1 boost.
 */
inline double gamma_deviate(double k, LecuyerRandom &rng) {
  if (!(k > 0.0)) return 1.0;
  if (k > 1e6) return std::max(0.0, 1.0 + rng.gaussian() / std::sqrt(k));
  double boost = 1.0;
  double a = k;
  if (a < 1.0) {
    boost = std::pow(rng(), 1.0 / a);
    a += 1.0;
  }
  const double dd = a - 1.0 / 3.0;
  const double c = 1.0 / std::sqrt(9.0 * dd);
  while (true) {
    double x = 0.0;
    double v = 0.0;
    do {
      x = rng.gaussian();
      v = 1.0 + c * x;
    } while (v <= 0.0);
    v = v * v * v;
    const double u = rng();
    if (u < 1.0 - 0.0331 * x * x * x * x || std::log(u) < 0.5 * x * x + dd * (1.0 - v + std::log(v))) {
      return boost * dd * v / k;
    }
  }
}

//! Unit vector perpendicular to (cx, cy, cz) at azimuth phi.
inline void perpendicular(double cx, double cy, double cz, double cos_phi,
                          double sin_phi, double &px, double &py, double &pz) {
  const double s = std::sqrt(std::max(0.0, 1.0 - cx * cx));
  if (s < 1e-10) {
    px = 0.0;
    py = cos_phi;
    pz = sin_phi;
    return;
  }
  px = -s * cos_phi;
  py = (cx * cy * cos_phi - cz * sin_phi) / s;
  pz = (cx * cz * cos_phi + cy * sin_phi) / s;
}

//! Rotate (cx, cy, cz) by polar angle (cos_t, sin_t) at azimuth phi.
inline void rotate(double &cx, double &cy, double &cz, double cos_t, double sin_t,
                   double cos_phi, double sin_phi) {
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  perpendicular(cx, cy, cz, cos_phi, sin_phi, px, py, pz);
  double nx = cx * cos_t + px * sin_t;
  double ny = cy * cos_t + py * sin_t;
  double nz = cz * cos_t + pz * sin_t;
  const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  cx = nx / norm;
  cy = ny / norm;
  cz = nz / norm;
}

//! rotate() with one division instead of five (fast approximations: the
//! results differ from rotate() in the last bits).
inline void rotate_fast(double &cx, double &cy, double &cz, double cos_t, double sin_t,
                        double cos_phi, double sin_phi) {
  const double s = std::sqrt(std::max(0.0, 1.0 - cx * cx));
  double px = 0.0;
  double py = cos_phi;
  double pz = sin_phi;
  if (!(s < 1e-10)) {
    const double inv_s = 1.0 / s;
    px = -s * cos_phi;
    py = (cx * cy * cos_phi - cz * sin_phi) * inv_s;
    pz = (cx * cz * cos_phi + cy * sin_phi) * inv_s;
  }
  const double nx = cx * cos_t + px * sin_t;
  const double ny = cy * cos_t + py * sin_t;
  const double nz = cz * cos_t + pz * sin_t;
  // c and p are unit vectors, p perpendicular to c, and cos^2 + sin^2 = 1,
  // so |n| = 1 up to rounding: one Newton step of 1/sqrt about 1 corrects
  // it to ~1e-22 without sqrt or division (the exact form otherwise).
  const double n2 = nx * nx + ny * ny + nz * nz;
  const double inv_norm = std::abs(n2 - 1.0) < 1e-6 ? 1.5 - 0.5 * n2 : 1.0 / std::sqrt(n2);
  cx = nx * inv_norm;
  cy = ny * inv_norm;
  cz = nz * inv_norm;
}

// DIRCOS's branch for a direction along the x axis (sin x = 0): the new
// direction is (cx cos t, -sin t cos phi, -sin t sin phi), normalised. TRIM.exe
// uses it for every atom struck by a recoil, whatever the mover's direction.
template <bool Fast = false>
inline void rotate_about_x_2013(double &cx, double &cy, double &cz, double cos_t, double sin_t,
                                double cos_phi, double sin_phi) {
  const double nx = cx * cos_t;
  const double ny = -sin_t * cos_phi;
  const double nz = -sin_t * sin_phi;
  const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  if (norm == 0.0) return;
  if constexpr (Fast) {
    const double inv = 1.0 / norm;
    cx = nx * inv;
    cy = ny * inv;
    cz = nz * inv;
  } else {
    cx = nx / norm;
    cy = ny / norm;
    cz = nz / norm;
  }
}

class Transporter {
public:
  Transporter(const TrimConfig &cfg, const TargetModel &m, TrimTallies &t, double bin0,
              double bin_w)
      : cfg_(cfg), m_(m), t_(t), bin0_(bin0), inv_w_(1.0 / bin_w),
        nbins_(cfg.depth_bins),
        ion_stop_(cfg.ion_stop_energy > 0.0 * electron_volt ? cfg.ion_stop_energy
                                                            : cfg.final_energy),
        recoil_stop_(cfg.recoil_stop_energy > 0.0 * electron_volt ? cfg.recoil_stop_energy
                                                                  : cfg.final_energy),
        replacement_(cfg.final_energy),
        thickness_(m.thickness),
        mode_(cfg.damage) {
    w0_ = in_units(cfg.window_start, angstrom);
    w1_ = cfg.window_end > 0.0 * angstrom ? in_units(cfg.window_end, angstrom)
                                          : in_units(m.thickness, angstrom);
    // MD events (TrimConfig::md_events): the transition energies, and one
    // totals entry per E_c in this block's tallies.
    md_on_ = cfg.md_events.enabled();
    if (md_on_) {
      md_handoff_ = cfg.md_events.handoff;
      for (const auto &e : cfg.md_events.transition_energies) md_ec_eV_.push_back(in_units(e, electron_volt));
      md_emin_eV_ = in_units(cfg.md_events.min_energy, electron_volt);
      t_.md.resize(md_ec_eV_.size());
      for (std::size_t i = 0; i < md_ec_eV_.size(); ++i) t_.md[i].transition_energy_eV = md_ec_eV_[i];
    }
    // Bin edges bin0 + k / inv_w for k = 0 .. nbins: the expression deposit()
    // evaluated twice per bin, computed once (the same doubles).
    edges_.resize(nbins_ + 1);
    for (std::size_t k = 0; k <= nbins_; ++k) edges_[k] = bin0_ + static_cast<double>(k) / inv_w_;
    exyz_interval_eV_ = in_units(cfg.exyz_interval, electron_volt);
    slots_.resize(2);
  }

  void run_ion(std::size_t ion_index, std::vector<CollisionPoint> *trajectory);
  //! Where the MD events of the following ions go (null: only their totals).
  void set_md_output(std::vector<MdEvent> *events) noexcept { md_out_ = events; }
  //! Crystal transport: the last ion's entry point on the beam spot, in the
  //! output axes (TRIM's y, z): the lateral offset of its lattice.
  std::array<double, 2> lattice_offset_A() const noexcept { return {spot_y_, spot_z_}; }

private:
  // Bin index of a depth (Angstrom), or -1 outside the window.
  long bin_of(double x) const {
    const double f = (x - bin0_) * inv_w_;
    // Checked before the conversion (NaN or a huge f has no long value);
    // for finite f the same bins as (long)f < nbins.
    if (!(f >= 0.0) || f >= static_cast<double>(nbins_)) return -1;
    return static_cast<long>(f);
  }
  void add(DepthHistogram &h, double x, double v) const {
    const long b = bin_of(x);
    if (b >= 0) h.values[static_cast<std::size_t>(b)] += v;
  }
  //! Bin of a recoil energy er >= 1 eV in the log-spaced spectra
  //! (spectrum_bins_per_decade from 1 eV), or -1 (NaN, beyond 10^6 bins).
  long spectrum_bin(double er) const {
    const double l = cfg_.fast_math ? physics::fast::log10_at_least_one(er) : std::log10(er);
    const double f = std::floor(l * TrimTallies::spectrum_bins_per_decade);
    return f >= 0.0 && f < 1e6 ? static_cast<long>(f) : -1;
  }
  // Distribute \p amount uniformly over the depth interval [x0, x1].
  void deposit(DepthHistogram &h, double x0, double x1, double amount) const {
    if (amount == 0.0) return;
    if (x1 < x0) std::swap(x0, x1);
    const double len = x1 - x0;
    if (len <= 1e-12) {
      add(h, x0, amount);
      return;
    }
    const double lo = std::max(x0, bin0_);
    const double hi = std::min(x1, edges_[nbins_]);
    if (!(hi > lo)) return; // (also NaN)
    long b0 = static_cast<long>((lo - bin0_) * inv_w_);
    long b1 = static_cast<long>((hi - bin0_) * inv_w_);
    b0 = std::clamp(b0, 0L, static_cast<long>(nbins_) - 1);
    b1 = std::clamp(b1, 0L, static_cast<long>(nbins_) - 1);
    if (cfg_.fast_math) {
      // One division per flight instead of one per bin.
      const double per_length = amount / len;
      for (long b = b0; b <= b1; ++b) {
        const double e0 = std::max(lo, edges_[static_cast<std::size_t>(b)]);
        const double e1 = std::min(hi, edges_[static_cast<std::size_t>(b) + 1]);
        if (e1 > e0) h.values[static_cast<std::size_t>(b)] += per_length * (e1 - e0);
      }
      return;
    }
    for (long b = b0; b <= b1; ++b) {
      const double e0 = std::max(lo, edges_[static_cast<std::size_t>(b)]);
      const double e1 = std::min(hi, edges_[static_cast<std::size_t>(b) + 1]);
      if (e1 > e0) h.values[static_cast<std::size_t>(b)] += amount * (e1 - e0) / len;
    }
  }

  void follow(Mover mv, LecuyerRandom &rng, std::size_t ion_index,
              std::vector<CollisionPoint> *trajectory);
  void stop_mover(const Mover &mv, bool interstitial, std::size_t ion_index);
  void leave_front(Mover &mv, std::size_t ion_index);
  void leave_back(const Mover &mv, std::size_t ion_index);
  void record_point(std::vector<CollisionPoint> *trajectory, const Mover &mv,
                    CollisionPoint::Kind kind = CollisionPoint::Kind::collision);

  double kp_damage_energy(const ElementModel &e, double T_eV) const {
    return robinson_damage(e.kp_kd, e.kp_eps_per_eV, T_eV);
  }

  /*!
   * Kinematics of a collision at lab energy \p e, reduced energy \p eps and
   * reduced impact parameter \p b: TRIM's magic formula (or the scattering
   * integral), or with TrimConfig::fast_math its table (MagicTable).
   */
  physics::CollisionKinematics scatter(const CollisionPair &pair, energy_quantity e, double eps,
                                       double b) const {
    if (cfg_.fast_math) {
      double sin_half = 0.0;
      double cos_half = 1.0;
      physics::fast_half_angle(pair, eps, b, cfg_.scattering, sin_half, cos_half);
      return physics::collision_kinematics_half_angle(pair, e, sin_half, cos_half);
    }
    return physics::collision_kinematics(
        pair, e, physics::trim_sin2_half_theta(pair.screening, eps, b, cfg_.scattering));
  }
  //! A collision's azimuth phi = 2 pi R: its cosine and sine.
  void azimuth(LecuyerRandom &rng, double &cos_phi, double &sin_phi) const {
    const double u = rng();
    if (cfg_.fast_math) {
      physics::fast::sincos_turn(u, cos_phi, sin_phi);
      return;
    }
    const double phi = 2.0 * pi_v * u;
    cos_phi = std::cos(phi);
    sin_phi = std::sin(phi);
  }
  //! Electronic stopping N S_e(E): the curve's table with fast_math.
  common::stopping_power_quantity stopping_at(const physics::StoppingCurve &c, energy_quantity e) const {
    return cfg_.fast_math ? c.stopping_power_fast(e) : c.stopping_power(e);
  }
  //! rotate(), or rotate_fast() with fast_math.
  void turn(double &cx, double &cy, double &cz, double cos_t, double sin_t, double cos_phi,
            double sin_phi) const {
    if (cfg_.fast_math) rotate_fast(cx, cy, cz, cos_t, sin_t, cos_phi, sin_phi);
    else rotate(cx, cy, cz, cos_t, sin_t, cos_phi, sin_phi);
  }

  const TrimConfig &cfg_;
  const TargetModel &m_;
  TrimTallies &t_;
  double bin0_;
  double inv_w_;
  std::vector<double> edges_; // depth-bin edges (nbins_ + 1), see the constructor
  std::size_t nbins_;
  energy_quantity ion_stop_;
  energy_quantity recoil_stop_;
  energy_quantity replacement_;
  length_quantity thickness_;
  DamageCalculation mode_;
  void cascade(const Mover &primary, std::size_t ion_index,
               std::vector<CollisionPoint> *trajectory);

  std::vector<Mover> stack_;
  std::size_t recoils_this_ion_ = 0;
  std::size_t cascades_this_ion_ = 0; // fixed cascade-random-replay bug: stream number
  LecuyerRandom *recoil_rng_ = nullptr;
  bool details_ = false; // record CascadeRecords for the current ion
  // The current cascade's displacements and replacements (COLLISON.txt
  // with recoils), or null.
  std::vector<CascadeAtom> *cascade_atoms_ = nullptr;
  void note_atom(unsigned z, double energy_eV, const Mover &where, bool vacancy) {
    if (cascade_atoms_ == nullptr) return;
    cascade_atoms_->push_back({static_cast<unsigned>(cascade_atoms_->size() + 1), z, energy_eV,
                               in_units(where.x, angstrom), in_units(where.y, angstrom),
                               in_units(where.zz, angstrom), vacancy, !vacancy});
  }
  // --- Damage models of the primary knock-on atoms (TrimTallies) ---------
  //! TRIM's damage energy of a recoil of energy T in \p layer: the
  //! layer-averaged Robinson partition of the quick calculation (kp_2013,
  //! 0x44bae0), T / (1 + k g(eps)).
  double layer_damage_energy(const LayerModel &layer, double T_eV) const {
    const auto &t = layer.t13;
    const double eps = t.kp_eps_per_eV * T_eV;
    const double g = cfg_.fast_math ? physics::fast::kp_partition(eps)
                                    : eps + 0.4 * std::pow(eps, 0.75) + 3.4 * std::pow(eps, 1.0 / 6.0);
    return T_eV / (1.0 + t.kp_k * g);
  }
  //! A primary knock-on atom of element \p el with damage energy \p td_eV
  //! at depth \p x: NRT, original Kinchin-Pease, arc-dpa and rpa (Nordlund
  //! et al. 2018, Eqs. 1, 6, 7, 9). Only tallies; no effect on the transport.
  void book_pka(const ElementModel &el, double x, double td_eV) {
    const double ed = in_units(el.ed, electron_volt);
    if (!(td_eV > 0.0) || !(ed > 0.0)) return;
    const double knee = 2.0 * ed / 0.8; // NRT's 2 Ed / 0.8
    const double nrt = td_eV < ed ? 0.0 : td_eV < knee ? 1.0 : 0.8 * td_eV / (2.0 * ed);
    const double kp = td_eV < ed ? 0.0 : td_eV < 2.0 * ed ? 1.0 : td_eV / (2.0 * ed);
    add(t_.pka_damage_energy, x, td_eV);
    add(t_.nrt_displacements, x, nrt);
    add(t_.kp_displacements, x, kp);
    t_.total_pka_damage_energy += td_eV;
    t_.total_nrt += nrt;
    t_.total_kp += kp;
    if (!el.has_efficiency) {
      t_.pka_without_efficiency += 1.0;
      return;
    }
    double arc = nrt;
    double rpa = nrt;
    if (td_eV >= knee) {
      const double xi_arc = el.arc_scale * std::pow(td_eV, el.arc_b) + el.arc_c;
      const double tc = std::pow(td_eV, el.rpa_c);
      const double xi_rpa = el.rpa_scale * tc / (el.rpa_bc + tc);
      arc = nrt * xi_arc;
      rpa = nrt * xi_rpa;
    }
    add(t_.arc_displacements, x, arc);
    add(t_.rpa_replacements, x, rpa);
    t_.total_arc += arc;
    t_.total_rpa += rpa;
  }

  // --- Trajectory bookkeeping and MD events (opt-in; not TRIM) ------------
  // None of this draws random numbers or changes a tally, except that the
  // MD handoff mode stops following the atoms it hands to MD.
  bool timing_ = false;            // advance Mover::time (trajectories or MD events recorded)
  std::uint64_t ion_number_ = 0;   // 1-based, for CollisionPoint::ion and MdEvent::ion
  bool md_on_ = false;
  bool md_handoff_ = false;
  std::vector<double> md_ec_eV_;   // E_c, ascending
  double md_emin_eV_ = 0.0;        // E_min (0: the struck element's Ed)
  std::vector<MdEvent> *md_out_ = nullptr;

  //! Advance \p m's clock over a flight of \p L at energy \p e.
  void advance_time(Mover &m, length_quantity L, energy_quantity e) const {
    if (!timing_) return;
    const common::velocity_quantity v = physics::relativistic_speed(e, m_.species[m.species].mass);
    if (v > 0.0 * meters_per_second) m.time += L / v;
  }
  //! Index of \p m's last recorded point (it must be the trajectory's last).
  std::uint32_t last_point_of(const std::vector<CollisionPoint> *trajectory, const Mover &m) const {
    if (trajectory == nullptr || trajectory->empty() || m.untracked) return CollisionPoint::no_parent;
    const auto &pt = trajectory->back();
    return pt.track == m.track ? pt.point_index : CollisionPoint::no_parent;
  }
  //! The mover's collision point \p point (refined and crystal transports)
  //! is where it ends: a stop, at the mover's position (a replacing atom
  //! sits on the site it took).
  void end_point(std::vector<CollisionPoint> *trajectory, std::uint32_t point, const Mover &m) const {
    if (trajectory == nullptr || point == CollisionPoint::no_parent || point >= trajectory->size()) return;
    auto &pt = (*trajectory)[point];
    if (pt.track != m.track) return;
    pt.coordinate = Coordinate(m.y, m.zz, m.x);
    pt.kind = CollisionPoint::Kind::stop;
  }
  //! Mark \p m's last recorded point as a stop or exit.
  void mark_last(std::vector<CollisionPoint> *trajectory, const Mover &m, CollisionPoint::Kind kind) const {
    if (trajectory == nullptr || trajectory->empty() || m.untracked) return;
    auto &pt = trajectory->back();
    if (pt.track == m.track) pt.kind = kind;
  }
  //! A recoil's collision point (recorded where it arrived): the energies
  //! of the collision, and the energy and direction it leaves with.
  void annotate_collision(std::vector<CollisionPoint> *trajectory, const Mover &before, const Mover &after,
                          energy_quantity T, energy_quantity de) const {
    if (trajectory == nullptr || trajectory->empty() || before.untracked) return;
    auto &pt = trajectory->back();
    if (pt.track != before.track) return;
    pt.energy_before_eV = in_units(before.energy, electron_volt);
    pt.energy_transferred_eV = in_units(T, electron_volt);
    pt.electronic_loss_eV = in_units(de, electron_volt);
    pt.energy_eV = in_units(after.energy, electron_volt);
    const double speed = in_units(physics::relativistic_speed(after.energy, m_.species[after.species].mass),
                                  meters_per_second);
    pt.vx = after.cy * speed;
    pt.vy = after.cz * speed;
    pt.vz = after.cx * speed;
  }
  void md_record(std::size_t i, MdEvent::Kind kind, const Mover &m, double energy_eV) {
    auto &st = t_.md[i];
    st.events += 1.0;
    (kind == MdEvent::Kind::birth ? st.birth_events : st.crossing_events) += 1.0;
    st.energy_eV += energy_eV;
    if (md_out_ == nullptr) return;
    MdEvent ev;
    ev.ion = static_cast<std::size_t>(ion_number_);
    ev.transition_energy_eV = md_ec_eV_[i];
    ev.kind = kind;
    ev.time_s = m.time.value(); // SI: seconds
    // The trajectory axes: z = depth (TRIM's x), x and y lateral.
    ev.x_A = in_units(m.y, angstrom);
    ev.y_A = in_units(m.zz, angstrom);
    ev.z_A = in_units(m.x, angstrom);
    ev.ux = m.cy;
    ev.uy = m.cz;
    ev.uz = m.cx;
    ev.energy_eV = energy_eV;
    ev.generation = m.generation;
    ev.z = m.z;
    ev.mass_amu = (m_.species[m.species].mass / atomic_mass_unit).value();
    ev.track = m.track;
    ev.parent_track = m.parent_track;
    md_out_->push_back(ev);
  }
  //! A recoil \p r set in motion with energy \p T_eV (its position,
  //! direction, time and lineage set; md_mask = its parent's): events for
  //! E_min <= T < E_c, below-E_min counts; returns true if the handoff mode
  //! must not follow it.
  bool md_birth(Mover &r, double T_eV, const ElementModel &el) {
    const double emin = md_emin_eV_ > 0.0 ? md_emin_eV_ : in_units(el.ed, electron_volt);
    bool event = false;
    for (std::size_t i = 0; i < md_ec_eV_.size(); ++i) {
      const std::uint32_t bit = 1u << i;
      if ((r.md_mask & bit) != 0 || T_eV >= md_ec_eV_[i]) continue; // already MD's, or followed
      r.md_mask |= bit;
      if (T_eV < emin) {
        t_.md[i].below_min += 1.0;
        t_.md[i].below_min_energy_eV += T_eV;
        continue;
      }
      md_record(i, MdEvent::Kind::birth, r, T_eV);
      event = true;
    }
    return md_handoff_ && event;
  }
  //! A transfer TRIM treats as sub-threshold (T <= Ed from the ion, whose
  //! cascade TRIM still runs with a stale energy): never an event; counted
  //! below E_min. The handoff mode does not follow it.
  bool md_below_min(Mover &r, double T_eV) {
    for (std::size_t i = 0; i < md_ec_eV_.size(); ++i) {
      const std::uint32_t bit = 1u << i;
      if ((r.md_mask & bit) != 0) continue;
      r.md_mask |= bit;
      t_.md[i].below_min += 1.0;
      t_.md[i].below_min_energy_eV += T_eV;
    }
    return md_handoff_;
  }
  //! \p m's energy fell from \p before_eV to m.energy: events where it
  //! crossed an E_c; returns true if the handoff mode stops following it.
  bool md_crossing(Mover &m, double before_eV) {
    const double now = in_units(m.energy, electron_volt);
    bool event = false;
    for (std::size_t i = 0; i < md_ec_eV_.size(); ++i) {
      const std::uint32_t bit = 1u << i;
      if ((m.md_mask & bit) != 0 || !(before_eV >= md_ec_eV_[i] && now < md_ec_eV_[i])) continue;
      m.md_mask |= bit;
      md_record(i, MdEvent::Kind::crossing, m, now);
      if (m.is_ion) t_.md[i].ions_handed += 1.0;
      event = true;
    }
    return md_handoff_ && event;
  }
  //! An ion that enters below an E_c is an event at its start.
  bool md_ion_entry(Mover &m) {
    const double e = in_units(m.energy, electron_volt);
    bool event = false;
    for (std::size_t i = 0; i < md_ec_eV_.size(); ++i) {
      if (e >= md_ec_eV_[i]) continue;
      m.md_mask |= 1u << i;
      md_record(i, MdEvent::Kind::birth, m, e);
      t_.md[i].ions_handed += 1.0;
      event = true;
    }
    return md_handoff_ && event;
  }

  // Totals at the start of the current ion (COLLISON.txt's ion summary).
  struct IonStart {
    double d, r, v, s, t;
  };
  IonStart ion_start() const {
    double s = 0.0, t = 0.0;
    for (double v : t_.sputtered) s += v;
    for (double v : t_.transmitted_recoils) t += v;
    return {t_.total_displacements, t_.total_replacements, t_.total_vacancies, s, t};
  }
  void ion_summary(std::size_t ion_index, const IonStart &a) {
    if (!details_) return;
    const IonStart b = ion_start();
    // TRIM counts every displaced atom that stays in the target as an
    // interstitial (replacing atoms included): displacements - sputtered -
    // transmitted atoms.
    const double d = b.d - a.d, sp = b.s - a.s, tr = b.t - a.t;
    t_.ion_summaries.push_back({ion_index + 1, d, b.r - a.r, b.v - a.v, d - sp - tr, sp, tr});
  }

  // --- Crystal transport (TrimConfig::Transport::crystal; not TRIM) -------
  // After Robinson and Torrens, Phys. Rev. B 9, 5008 (1974) (MARLOWE);
  // docs/CRYSTAL_DESIGN.md.
  //! One atom struck in a collision of the crystal transport.
  //! A site of one of the target's lattices (LayerModel::domain, plus the grain).
  struct SiteKey {
    std::uint32_t domain = 0;
    crystal::SiteId id;
    bool operator==(const SiteKey &) const = default;
  };
  struct SiteKeyHash {
    std::size_t operator()(const SiteKey &k) const noexcept {
      return static_cast<std::size_t>(crystal::site_hash(k.id, k.domain));
    }
  };
  struct CrystalTarget {
    std::size_t layer = 0;       // layer of the atom
    std::size_t element = 0;     // its element in that layer
    double pos[3] = {0, 0, 0};   // the atom (laboratory, Angstrom)
    double xi = 0.0;             // path length to the closest approach (Angstrom)
    double impact = 0.0;         // impact parameter (Angstrom)
    bool on_site = false;        // a lattice site (id) or a random atom of an amorphous layer
    SiteKey id;
  };
  void run_ion_crystal(std::size_t ion_index, std::vector<CollisionPoint> *trajectory);
  void run_source_crystal(std::size_t ion_index, LecuyerRandom &ion_rng, std::vector<CollisionPoint> *trajectory);
  void crystal_follow(Mover mv, LecuyerRandom &rng, std::size_t ion_index,
                      std::vector<CollisionPoint> *trajectory);
  bool crystal_strike(const Mover &mv, const CrystalTarget &target, const double recoil_dir[3], energy_quantity T,
                      energy_quantity e1, bool may_replace, std::uint32_t parent_point);
  energy_quantity crystal_path_loss(const Mover &mv, const SpeciesLayer &d, double length_A, double share,
                                    LecuyerRandom &rng);
  double crystal_path_share(std::size_t species, std::size_t li);
  std::uint64_t history_seed_ = 0;  // the frozen lattice configuration of the current ion
  double spot_y_ = 0.0;             // the ion's entry point on the surface (Angstrom)
  double spot_z_ = 0.0;
  // Per-grain results (TrimTallies::grains): the ion's first grain, the
  // grain it was in at its last step (-1: none), and whether and how deep it
  // stopped (set by stop_mover). Bookkeeping only: no random number, no
  // other tally.
  int ion_entry_grain_ = -1;
  std::uint32_t ion_entry_layer_ = 0;
  int ion_last_grain_ = -1;
  std::uint32_t ion_last_layer_ = 0;
  bool ion_stopped_ = false;
  double ion_stop_depth_A_ = 0.0;
  void book_grains();
  std::unordered_set<SiteKey, SiteKeyHash> vacated_; // sites emptied in this history
  std::vector<double> path_share_;  // per species x layer: electronic loss along the path (cache)

  // --- TRIM-2013 transport (TrimConfig::Transport::trim2013) ---------------
  // See docs/TRIM2013_ALGORITHM.md for the TRIM.exe routines each step
  // reproduces.
  void run_ion_2013(std::size_t ion_index, std::vector<CollisionPoint> *trajectory);
  length_quantity ion_flight_path_2013(const Mover &mv, std::size_t li, LecuyerRandom &rng);
  void edge_2013(const Mover &mv, std::size_t li, length_quantity &L, LecuyerRandom &rng);
  void interface_2013(const Mover &mv, std::size_t li, length_quantity &L, LecuyerRandom &rng);
  //! \p known_se: the stopping power at mv.energy if the caller has it
  //! (the same value is then not computed twice).
  energy_quantity electronic_loss_2013(const Mover &mv, std::size_t li, length_quantity L,
                                       LecuyerRandom &rng,
                                       const common::stopping_power_quantity *known_se = nullptr) const;
  std::size_t pick_element(const LayerModel &layer, LecuyerRandom &rng) const;
  void kp_2013(const ElementModel &el, const LayerModel &layer, double T_eV, double x,
               std::size_t ion_index, const Mover &ion, common::stopping_power_quantity se);
  void displace_2013(const ElementModel &el, double x, bool by_ion, double T_eV);
  void cascade_2013(std::size_t ion_index, std::vector<CollisionPoint> *trajectory);
  void recoil_stop_2013(std::size_t k);
  bool in_window(double x) const { return x >= w0_ && x <= w1_; }

  void run_recoil_event_2013(std::size_t ion_index, std::vector<CollisionPoint> *trajectory);
  void start_cascade_2013(const ElementModel &el, std::size_t li, std::size_t j,
                          const Mover &site, energy_quantity T, std::size_t ion_index,
                          std::vector<CollisionPoint> *trajectory, LecuyerRandom &rng,
                          bool restore_rng, double ion_energy_keV,
                          common::stopping_power_quantity se);
  void record_exyz(std::size_t ion_index, const Mover &mv, double se_eV_per_A);
  void surface_arrival(const Mover &r);

  double w0_ = 0.0;            // TRIM's depth window [Xmin, Xmax] (Angstrom)
  double w1_ = 0.0;
  std::uint32_t next_track_ = 0; // last path number given out for this ion
  double exyz_interval_eV_ = 0.0;
  long exyz_next_ = 0;         // TRIM's G360: energy step of the last EXYZ record
  double last_T_eV_ = 0.0;     // TRIM's G5DC: energy given to the last struck atom
  int flight_counter_ = 0;     // TRIM's G518: randomised flights left
  std::vector<Mover> slots_;   // TRIM's recoil arrays (persist between ions)
  bool primary_displaced_ = false; // TRIM's G674
  // TRIM.DAT recoils: the depth where the atom of the event stops.
  bool track_primary_stop_ = false;
  double primary_stop_x_ = 0.0;
  bool sub_threshold_ = false;     // TRIM's G60C: current mover never displaced
};

void Transporter::record_point(std::vector<CollisionPoint> *trajectory, const Mover &mv,
                               CollisionPoint::Kind kind) {
  if (trajectory == nullptr || mv.untracked) return;
  // AMCSET convention: +z is depth; x, y lateral.
  CollisionPoint p{Coordinate(mv.y, mv.zz, mv.x)};
  const double speed = in_units(physics::relativistic_speed(mv.energy, m_.species[mv.species].mass),
                                meters_per_second);
  p.vx = mv.cy * speed;
  p.vy = mv.cz * speed;
  p.vz = mv.cx * speed;
  p.generation = mv.generation;
  p.z_number = mv.z;
  p.atomic_mass_number = (m_.species[mv.species].mass / atomic_mass_unit).value();
  p.energy_eV = in_units(mv.energy, electron_volt);
  p.track = mv.track;
  p.kind = kind;
  p.ion = ion_number_;
  p.point_index = static_cast<std::uint32_t>(trajectory->size());
  p.parent_track = mv.parent_track;
  p.parent_point = mv.parent_point;
  p.time = mv.time.value(); // SI: seconds (0 unless timing_)
  // The entry point on the beam spot (crystal transport; 0 otherwise), in
  // the same axes as the coordinate: x = TRIM's y, y = TRIM's z.
  p.entry_x_A = spot_y_;
  p.entry_y_A = spot_z_;
  trajectory->push_back(p);
}

void Transporter::stop_mover(const Mover &mv, bool interstitial, std::size_t ion_index) {
  const double x = in_units(mv.x, angstrom);
  const double e = in_units(mv.energy, electron_volt);
  if (mv.is_ion) {
    const double y = in_units(mv.y, angstrom);
    const double z = in_units(mv.zz, angstrom);
    const double r = std::sqrt(y * y + z * z);
    add(t_.ion_stopped, x, 1.0);
    t_.range.add(x);
    t_.lateral_y.add(y);
    t_.lateral_abs_y.add(std::abs(y));
    t_.radial.add(r);
    const long b = bin_of(x);
    if (b >= 0) {
      auto &lb = t_.lateral_bins[static_cast<std::size_t>(b)];
      lb[0] += 1.0;
      lb[1] += std::abs(y);
      lb[2] += y * y;
      lb[3] += r;
      lb[4] += r * r;
    }
    add(t_.phonons_by_ion, x, e);
    t_.total_phonons_ion += e;
    if (cfg_.record_ion_positions) t_.positions.push_back({ion_index + 1, x, y, z});
    ion_stopped_ = true;
    ion_stop_depth_A_ = x;
  } else {
    if (interstitial && !mv.surface_only) {
      add(t_.recoil_stopped[mv.tally_element], x, 1.0);
      add(t_.interstitials, x, 1.0);
      t_.total_interstitials += 1.0;
    }
    add(t_.phonons_by_recoils, x, e);
    t_.total_phonons_recoils += e;
  }
}

void Transporter::leave_front(Mover &mv, std::size_t ion_index) {
  // Move back to the surface plane along the trajectory.
  if (mv.cx < 0.0) {
    const double s = in_units(mv.x, angstrom) / mv.cx; // negative distance
    mv.y -= s * mv.cy * angstrom;
    mv.zz -= s * mv.cz * angstrom;
  }
  mv.x = 0.0 * angstrom;
  const double e = in_units(mv.energy, electron_volt);
  if (mv.is_ion) {
    t_.backscattered += 1;
    t_.backscattered_energy += e;
    if (cfg_.record_exiting_atoms && t_.exits.size() < cfg_.max_exit_records) {
      t_.exits.push_back({'B', ion_index + 1, mv.z, e, 0.0, in_units(mv.y, angstrom),
                          in_units(mv.zz, angstrom), mv.cx, mv.cy, mv.cz});
    }
    return;
  }
  // Planar surface potential: the normal energy must exceed Es.
  const auto &layer = m_.layers.front();
  energy_quantity es = 0.0 * electron_volt;
  for (const auto &el : layer.elements) {
    if (el.tally_index == mv.tally_element) es = el.es;
  }
  if (es == 0.0 * electron_volt) {
    // The recoil's element is not in the top layer: use its own Es from
    // whichever layer defines it.
    for (const auto &l : m_.layers) {
      for (const auto &el : l.elements) {
        if (el.tally_index == mv.tally_element) es = el.es;
      }
    }
  }
  const energy_quantity normal = mv.energy * mv.cx * mv.cx;
  if (normal > es) {
    const energy_quantity out = mv.energy - es;
    const double scale = std::sqrt((mv.energy / out).value());
    const double ncx = -std::sqrt(std::max(0.0, ((normal - es) / out).value()));
    const double ncy = mv.cy * scale;
    const double ncz = mv.cz * scale;
    const double eo = in_units(out, electron_volt);
    t_.sputtered[mv.tally_element] += 1.0;
    t_.sputtered_energy[mv.tally_element] += eo;
    if (mv.surface_only) {
      // The surface atom left its site for good: one vacancy.
      add(t_.vacancies[mv.tally_element], mv.origin_x, 1.0);
      t_.total_vacancies += 1.0;
    }
    // The binding energy is left in the surface as phonons.
    add(t_.phonons_by_recoils, 0.0, in_units(es, electron_volt));
    t_.total_phonons_recoils += in_units(es, electron_volt);
    if (cfg_.record_exiting_atoms && t_.exits.size() < cfg_.max_exit_records) {
      t_.exits.push_back({'S', ion_index + 1, mv.z, eo, 0.0, in_units(mv.y, angstrom),
                          in_units(mv.zz, angstrom), ncx, ncy, ncz});
    }
  } else {
    // Reflected by the surface barrier and trapped at the surface (TRIM
    // does not count these as interstitials).
    add(t_.phonons_by_recoils, 0.0, e);
    t_.total_phonons_recoils += e;
  }
}

void Transporter::leave_back(const Mover &mv, std::size_t ion_index) {
  const double e = in_units(mv.energy, electron_volt);
  if (mv.is_ion) {
    t_.transmitted += 1;
    t_.transmitted_energy += e;
  } else {
    t_.transmitted_recoils[mv.tally_element] += 1.0;
  }
  if (cfg_.record_exiting_atoms && t_.exits.size() < cfg_.max_exit_records) {
    t_.exits.push_back({'T', ion_index + 1, mv.z, e, in_units(thickness_, angstrom),
                        in_units(mv.y, angstrom), in_units(mv.zz, angstrom), mv.cx, mv.cy, mv.cz});
  }
}

void Transporter::follow(Mover mv, LecuyerRandom &rng, std::size_t ion_index,
                         std::vector<CollisionPoint> *trajectory) {
  const energy_quantity stop_energy = mv.is_ion ? ion_stop_ : recoil_stop_;
  const bool full = mode_ != DamageCalculation::kinchin_pease;
  const bool monolayer_only = !cfg_.free_flight_path || mode_ == DamageCalculation::monolayer;
  record_point(trajectory, mv, CollisionPoint::Kind::birth);

  while (true) {
    if (mv.energy <= stop_energy) {
      stop_mover(mv, true, ion_index);
      record_point(trajectory, mv, CollisionPoint::Kind::stop);
      return;
    }
    if (!(mv.energy >= 0.0 * electron_volt) || !std::isfinite(in_units(mv.x, angstrom))) {
      throw std::runtime_error("a moving atom's energy or depth is not a number");
    }
    const std::size_t li = m_.layer_index(mv.x);
    const LayerModel &layer = m_.layers[li];
    const SpeciesLayer &d = m_.data[mv.species][li];
    const auto stopping_power = stopping_at(d.stopping, mv.energy);

    // --- Free flight path ------------------------------------------------
    length_quantity L = layer.monolayer;
    if (!monolayer_only) {
      // Largest maximum impact parameter over the layer's elements
      // (energy-transfer criterion, SRIM book Eqs. 7-28, 7-40).
      length_quantity pmax = 0.0 * angstrom;
      for (std::size_t j = 0; j < d.pairs.size(); ++j) {
        const auto &pair = d.pairs[j];
        const double eps = pair.reduced_energy(mv.energy);
        const double b = physics::trim_max_reduced_impact_parameter(eps, d.eps_min[j], pair.gamma);
        pmax = std::max(pmax, b * pair.screening_length);
      }
      L = 1.0 / (pi_v * pmax * pmax * layer.density);
      // Electronic energy loss over L limited to 5 % of E.
      if (stopping_power * L > 0.05 * mv.energy) L = 0.05 * mv.energy / stopping_power;
      if (L < layer.min_flight) L = layer.min_flight;
      // Do not jump over interfaces or surfaces: truncate to a random
      // fraction of the distance to the nearest edge (never below the
      // shortest flight), so that the final approach is in short steps.
      length_quantity edge = std::numeric_limits<double>::max() * angstrom;
      if (mv.cx > 1e-12) edge = (layer.bottom - mv.x) / mv.cx;
      else if (mv.cx < -1e-12) edge = (mv.x - layer.top) / (-mv.cx);
      if (L > edge) L = std::max(layer.min_flight, rng() * edge);
    }
    // One target atom per cylinder of length L and radius pmax.
    const length_quantity pmax = 1.0 / boost::units::sqrt(pi_v * L * layer.density);

    // --- Electronic energy loss over the flight ---------------------------
    // As in TRIM, the loss is evaluated at the energy at the start of the
    // flight and subtracted together with the nuclear loss of the collision
    // at the end of the flight; the collision itself sees the start energy.
    const energy_quantity e_start = mv.energy;
    energy_quantity de = stopping_power * L;
    if (cfg_.electronic_straggling && de > 0.0 * electron_volt) {
      const double le = std::log(std::max(in_units(mv.energy, electron_volt), 1.0));
      double fi = (le - d.strag_log_e0) * d.strag_inv_step;
      fi = std::clamp(fi, 0.0, static_cast<double>(d.strag.size() - 1) - 1e-9);
      const auto i = static_cast<std::size_t>(fi);
      const double w = fi - static_cast<double>(i);
      const double omega2_per_A = d.strag[i] + w * (d.strag[i + 1] - d.strag[i]);
      const double variance = omega2_per_A * in_units(L, angstrom);
      // Gamma-distributed loss with the right mean and variance (a Gaussian
      // truncated at zero would bias the mean upwards at low energies).
      const double mean = in_units(de, electron_volt);
      if (variance > 0.0) {
        if (cfg_.straggling_distribution == TrimConfig::StragglingDistribution::gamma) {
          de = mean * gamma_deviate(mean * mean / variance, rng) * electron_volt;
        } else {
          de = std::max(0.0, mean + std::sqrt(variance) * rng.gaussian()) * electron_volt;
        }
      }
    }
    if (de > mv.energy) de = mv.energy;
    const double x0 = in_units(mv.x, angstrom);
    mv.x += L * mv.cx;
    mv.y += L * mv.cy;
    mv.zz += L * mv.cz;
    mv.energy = e_start - de;
    advance_time(mv, L, e_start);
    {
      const double x1 = in_units(mv.x, angstrom);
      const double ev = in_units(de, electron_volt);
      // Only the part of the flight inside the target deposits energy.
      double a = x0;
      double b = x1;
      const double th = in_units(thickness_, angstrom);
      double frac = 1.0;
      if (std::abs(b - a) > 1e-12) {
        const double lo = std::max(std::min(a, b), 0.0);
        const double hi = std::min(std::max(a, b), th);
        frac = hi > lo ? (hi - lo) / std::abs(b - a) : 0.0;
        a = std::clamp(a, 0.0, th);
        b = std::clamp(b, 0.0, th);
      }
      const double inside = ev * frac;
      if (mv.is_ion) {
        deposit(t_.ionization_by_ion, a, b, inside);
        t_.total_ionization_ion += inside;
      } else {
        deposit(t_.ionization_by_recoils, a, b, inside);
        t_.total_ionization_recoils += inside;
      }
    }

    // --- Leaving the target? ---------------------------------------------
    if (mv.x < 0.0 * angstrom) {
      leave_front(mv, ion_index);
      record_point(trajectory, mv, CollisionPoint::Kind::exit);
      return;
    }
    if (mv.x > thickness_) {
      leave_back(mv, ion_index);
      record_point(trajectory, mv, CollisionPoint::Kind::exit);
      return;
    }
    if (mv.energy <= stop_energy) {
      stop_mover(mv, true, ion_index);
      record_point(trajectory, mv, CollisionPoint::Kind::stop);
      return;
    }

    // --- Collision partner -----------------------------------------------
    const length_quantity p = pmax * std::sqrt(rng());
    double cphi = 1.0;
    double sphi = 0.0;
    azimuth(rng, cphi, sphi);
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    perpendicular(mv.cx, mv.cy, mv.cz, cphi, sphi, px, py, pz);
    const length_quantity tx = mv.x + p * px;
    if (tx < 0.0 * angstrom || tx > thickness_) {
      continue; // no target atom outside the target
    }
    const std::size_t lt = m_.layer_index(tx);
    const LayerModel &tl = m_.layers[lt];
    std::size_t j = 0;
    if (tl.elements.size() > 1) {
      const double r = rng();
      while (j + 1 < tl.elements.size() && r > tl.elements[j].cumulative) ++j;
    }
    const ElementModel &el = tl.elements[j];
    const CollisionPair &pair = m_.data[mv.species][lt].pairs[j];
    const double eps = pair.reduced_energy(e_start);
    const double bred = (p / pair.screening_length).value();
    const auto k = scatter(pair, e_start, eps, bred);
    const energy_quantity T = std::min(k.transferred, mv.energy);
    const energy_quantity e1 = mv.energy - T;
    const double txa = in_units(tx, angstrom);
    const double T_eV = in_units(T, electron_volt);

    // The mover after the collision (deflected away from the struck atom:
    // azimuth phi + pi). Its trajectory point is recorded now, before the
    // points of the cascade this collision starts, so that they can name it.
    Mover after = mv;
    turn(after.cx, after.cy, after.cz, k.cos_projectile, k.sin_projectile, -cphi, -sphi);
    after.energy = e1;
    record_point(trajectory, after);
    {
      Mover before = mv;
      before.energy = e_start;
      annotate_collision(trajectory, before, after, T, de);
    }
    const std::uint32_t collision_point = last_point_of(trajectory, after);

    // Does the struck atom leave its site? Full cascades: T > Ed. Quick
    // (Kinchin-Pease) calculation: TRIM tests the damage energy instead,
    // so a primary recoil needs E_damage(T) >= Ed (e.g. T >= 45.3 eV for
    // Fe with Ed = 40 eV, checked against TRIM's COLLISON.txt).
    double kp_t = 0.0;
    double kp_damage = 0.0;
    bool displaces = false;
    if (full) {
      displaces = T > el.ed;
    } else {
      kp_t = T_eV;
      if (cfg_.kp_subtract_binding_energy) kp_t -= in_units(el.eb, electron_volt);
      kp_damage = kp_damage_energy(el, std::max(kp_t, 0.0));
      displaces = T > 0.0 * electron_volt && kp_damage >= in_units(el.ed, electron_volt);
    }

    // E2RECOIL.txt: energy given to atoms that recoil (displaced); smaller
    // transfers are phonons of the mover.
    if (mv.is_ion && displaces) {
      add(t_.energy_to_recoils, txa, T_eV);
      add(t_.energy_absorbed[el.tally_index], txa, T_eV);
      t_.total_energy_to_recoils += T_eV;
    }

    // The struck atom, if it will be followed: built before the bookkeeping,
    // with its lineage and clock, so that MD events are settled first (in
    // the MD handoff mode MD takes it, and the BCA books nothing for it, as
    // in the trim2013 transport). No random number is drawn here.
    const bool surface_atom = !displaces && full && tx < cfg_.surface_atom_depth && T > el.es &&
                              el.es > 0.0 * electron_volt;
    const bool follows = full && ((displaces && T - el.eb > 0.0 * electron_volt &&
                                   (cfg_.max_recoils_per_ion == 0 || recoils_this_ion_ < cfg_.max_recoils_per_ion)) ||
                                  surface_atom);
    Mover r;
    bool struck_handed = false;
    if (follows) {
      r.species = el.species;
      r.z = el.z;
      r.x = tx;
      r.y = mv.y + p * py;
      r.zz = mv.zz + p * pz;
      r.cx = mv.cx;
      r.cy = mv.cy;
      r.cz = mv.cz;
      // Recoil direction: angle (pi - theta)/2 from the incident
      // direction, towards the side of the struck atom (azimuth phi).
      turn(r.cx, r.cy, r.cz, k.cos_recoil, k.sin_recoil, cphi, sphi);
      r.energy = surface_atom ? T : T - el.eb;
      r.is_ion = false;
      r.surface_only = surface_atom;
      if (surface_atom) r.origin_x = txa;
      r.generation = mv.generation + 1;
      r.tally_element = el.tally_index;
      r.track = ++next_track_;
      r.parent_track = mv.track;
      r.parent_point = collision_point;
      r.time = mv.time;
      r.md_mask = mv.md_mask;
      struck_handed = md_on_ && md_birth(r, T_eV, el);
    }

    bool mover_stops = false;
    bool mover_interstitial = true;
    if (struck_handed) {
      // MD takes the struck atom (handoff mode).
      if (e1 <= stop_energy) mover_stops = true;
    } else if (displaces) {
      // Displacement.
      add(t_.displacements, txa, 1.0);
      t_.total_displacements += 1.0;
      if (mv.is_ion) add(t_.knock_ons, txa, 1.0);
      {
        const double er = std::max(T_eV - in_units(el.eb, electron_volt), 1.0);
        const long bin = spectrum_bin(er);
        if (bin >= 0 && bin < static_cast<long>(t_.recoil_spectrum.size())) {
          if (mv.is_ion) t_.pka_spectrum[static_cast<std::size_t>(bin)] += 1.0;
          t_.recoil_spectrum[static_cast<std::size_t>(bin)] += 1.0;
        }
      }
      if (!full) {
        // Quick calculation: modified Kinchin-Pease for the primary recoil.
        const double ev = kp_damage;
        const double ed = in_units(el.ed, electron_volt);
        double nu = 1.0;
        if (ev > 2.5 * ed) nu = 0.8 * ev / (2.0 * ed);
        add(t_.vacancies_by_ion, txa, 1.0);
        add(t_.vacancies_by_recoils, txa, nu - 1.0);
        add(t_.vacancies[el.tally_index], txa, nu);
        t_.total_vacancies += nu;
        if (details_ && mv.is_ion) {
          t_.cascades.push_back({ion_index + 1, in_units(mv.energy, kilo_electron_volt), txa,
                                 in_units(mv.y, angstrom), in_units(mv.zz, angstrom),
                                 in_units(stopping_power, electron_volt_per_angstrom), el.z, T_eV,
                                 nu, nu, 0.0, 0.0, {}});
        }
        add(t_.ionization_by_recoils, txa, T_eV - ev);
        t_.total_ionization_recoils += T_eV - ev;
        add(t_.phonons_by_recoils, txa, ev);
        t_.total_phonons_recoils += ev;
        if (e1 <= stop_energy) mover_stops = true;
      } else {
        const energy_quantity e2 = T - el.eb;
        t_.total_binding_energy_loss += in_units(el.eb, electron_volt);
        if (e1 < replacement_ && mv.z == el.z) {
          // Replacement collision: the mover takes the vacated site.
          add(t_.replacements, txa, 1.0);
          t_.total_replacements += 1.0;
          mover_stops = true;
          mover_interstitial = false;
        } else {
          // Vacancy; the mover moves on (and becomes an interstitial when
          // it falls below its stopping energy).
          add(t_.vacancies[el.tally_index], txa, 1.0);
          t_.total_vacancies += 1.0;
        }
        // Follow the recoil (built above; it is pushed and transported
        // after the mover).
        if (follows) {
          ++recoils_this_ion_;
          t_.recoils_followed += 1.0;
          if (mv.is_ion) {
            // Follow the primary recoil's whole cascade now (with the recoil
            // random stream, so the ion's trajectory is unaffected).
            const double d0 = t_.total_displacements;
            const double v0 = t_.total_vacancies;
            const double r0 = t_.total_replacements;
            const double i0 = t_.total_interstitials;
            cascade(r, ion_index, trajectory);
            if (details_) {
              t_.cascades.push_back({ion_index + 1, in_units(mv.energy, kilo_electron_volt),
                                     txa, in_units(mv.y, angstrom), in_units(mv.zz, angstrom),
                                     in_units(stopping_power, electron_volt_per_angstrom), el.z,
                                     T_eV, t_.total_displacements - d0 + 1.0,
                                     t_.total_vacancies - v0 + (e1 < replacement_ && mv.z == el.z ? 0.0 : 1.0),
                                     t_.total_replacements - r0 + (e1 < replacement_ && mv.z == el.z ? 1.0 : 0.0),
                                     t_.total_interstitials - i0, {}});
            }
          } else {
            stack_.push_back(r);
          }
        } else if (e2 > 0.0 * electron_volt) {
          // Recoil limit reached: its energy stays here as phonons.
          t_.recoil_limit_hits += 1.0;
          add(t_.phonons_by_recoils, txa, in_units(e2, electron_volt));
          t_.total_phonons_recoils += in_units(e2, electron_volt);
          add(t_.recoil_stopped[el.tally_index], txa, 1.0);
          add(t_.interstitials, txa, 1.0);
          t_.total_interstitials += 1.0;
        }
        if (e1 <= stop_energy) mover_stops = true;
      }
    } else if (full && tx < cfg_.surface_atom_depth && T > el.es &&
               el.es > 0.0 * electron_volt) {
      if (mv.is_ion) {
        add(t_.energy_to_recoils, txa, T_eV);
        add(t_.energy_absorbed[el.tally_index], txa, T_eV);
        t_.total_energy_to_recoils += T_eV;
      }
      // A surface atom set in motion below Ed (built above): it is followed
      // (it may be sputtered) but leaves no defect if it stays in the target.
      if (mv.is_ion) cascade(r, ion_index, trajectory); else stack_.push_back(r);
      if (e1 <= stop_energy) mover_stops = true;
    } else {
      // Sub-threshold transfer: the struck atom vibrates back (phonons).
      if (mv.is_ion) {
        add(t_.phonons_by_ion, txa, T_eV);
        t_.total_phonons_ion += T_eV;
      } else {
        add(t_.phonons_by_recoils, txa, T_eV);
        t_.total_phonons_recoils += T_eV;
      }
      if (e1 <= stop_energy) mover_stops = true;
    }

    // Deflect the mover (as recorded above).
    mv.cx = after.cx;
    mv.cy = after.cy;
    mv.cz = after.cz;
    mv.energy = e1;
    // MD events: the mover's energy may have fallen below an E_c here (the
    // handoff mode then leaves it to MD).
    if (md_on_ && md_crossing(mv, in_units(e_start, electron_volt))) return;

    if (mover_stops) {
      if (!mover_interstitial) {
        // Replacement: the mover ends on the struck atom's site.
        mv.x = tx;
        mv.y += p * py;
        mv.zz += p * pz;
      }
      end_point(trajectory, collision_point, mv);
      if (mv.is_ion) {
        stop_mover(mv, true, ion_index);
      } else {
        stop_mover(mv, mover_interstitial, ion_index);
      }
      return;
    }
  }
}

void Transporter::cascade(const Mover &primary, std::size_t ion_index,
                          std::vector<CollisionPoint> *trajectory) {
  stack_.push_back(primary);
  while (!stack_.empty()) {
    const Mover r = stack_.back();
    stack_.pop_back();
    follow(r, *recoil_rng_, ion_index, trajectory);
  }
}

void Transporter::run_ion(std::size_t ion_index, std::vector<CollisionPoint> *trajectory) {
  ion_number_ = ion_index + 1;
  timing_ = trajectory != nullptr || md_on_;
  if (cfg_.transport == TrimConfig::Transport::trim2013) {
    run_ion_2013(ion_index, trajectory);
    return;
  }
  if (cfg_.transport == TrimConfig::Transport::crystal) {
    run_ion_crystal(ion_index, trajectory);
    return;
  }
  const auto &ion = cfg_.ion;
  Mover mv;
  mv.species = 0;
  mv.z = ion.z;
  mv.x = 0.0 * angstrom;
  mv.y = 0.0 * angstrom;
  mv.zz = 0.0 * angstrom;
  mv.cx = std::cos(in_units(ion.angle, radian));
  mv.cy = std::sin(in_units(ion.angle, radian));
  mv.cz = 0.0;
  mv.energy = ion.energy;
  mv.is_ion = true;
  mv.generation = 0;
  mv.track = 0;
  next_track_ = 0;
  stack_.clear();
  recoils_this_ion_ = 0;
  cascades_this_ion_ = 0;

  LecuyerRandom ion_rng = make_stream(cfg_.seed, ion_index, 0);
  LecuyerRandom recoil_rng = make_stream(cfg_.seed, ion_index, 1);
  recoil_rng_ = &recoil_rng;
  details_ = ion_index < cfg_.collision_details_ions;
  // MD events: an ion entering below an E_c is an event at its start (the
  // handoff mode then leaves the whole ion to MD).
  if (md_on_ && md_ion_entry(mv)) {
    record_point(trajectory, mv, CollisionPoint::Kind::birth);
  } else {
    follow(mv, ion_rng, ion_index, trajectory);
  }
  recoil_rng_ = nullptr;
  t_.ions += 1;
}

// ---------------------------------------------------------------------------
// Crystal transport (TrimConfig::Transport::crystal). Not TRIM.
//
// The binary-collision scheme of M. T. Robinson and I. M. Torrens, Phys.
// Rev. B 9, 5008 (1974) (MARLOWE), Secs. II and IV, in the layers that have a
// crystal: a mover travels on straight segments from one deflection point to
// the next; its partners are the lattice atoms ahead of it within a maximum
// impact parameter (crystal::LatticeSearch), displaced by their thermal
// vibrations; atoms whose closest approaches nearly coincide along the path
// are struck together and the momentum changes added. Differences from that
// paper, each stated where it occurs: the potential and the scattering angle
// are AMCSET's (ZBL universal, magic formula or integral); the deflection
// point uses the hard-core value of the time integral; the electronic loss is
// AMCSET's stopping power divided between the path and the collisions.
// Amorphous layers take one random partner per monolayer step.
// ---------------------------------------------------------------------------

double Transporter::crystal_path_share(std::size_t species, std::size_t li) {
  // The part of the stopping power lost along the path. The collisions take
  // local_stopping_fraction of it as exp(-c R/a); integrated over impact
  // parameters up to s_max (straight paths, R = s) that is the fraction
  // 1 - (1 + x) exp(-x), x = c s_max / a, of a full normalisation: what the
  // collisions beyond s_max would have taken stays on the path, so that a
  // random direction loses the whole stopping power.
  if (path_share_.empty()) path_share_.assign(m_.species.size() * m_.layers.size(), -1.0);
  double &share = path_share_[species * m_.layers.size() + li];
  if (share >= 0.0) return share;
  const LayerModel &layer = m_.layers[li];
  const double f = cfg_.crystal.local_stopping_fraction;
  if (!layer.lattice || f <= 0.0) return share = 1.0;
  double collided = 0.0;
  for (std::size_t j = 0; j < layer.elements.size(); ++j) {
    const double a = in_units(m_.data[species][li].pairs[j].screening_length, angstrom);
    const double x = cfg_.crystal.local_stopping_decay * layer.lattice->max_impact_A() / a;
    collided += layer.elements[j].fraction * (1.0 - (1.0 + x) * std::exp(-x));
  }
  return share = 1.0 - f * collided;
}

energy_quantity Transporter::crystal_path_loss(const Mover &mv, const SpeciesLayer &d, double length_A,
                                               double share, LecuyerRandom &rng) {
  // Electronic loss along a straight flight of length_A from mv's position:
  // share S(E) L for the part of the flight inside the target, with the
  // straggling of the refined transport scaled by the same share. Always
  // gamma-distributed (exact mean): the flights here are the ~1 A steps from
  // one atom to the next, where the straggling exceeds the mean loss and
  // TRIM's normal deviate clipped at zero would raise the mean loss (by 9 %
  // for 100 keV He in Si, making ranges 7 % short).
  const double th = in_units(thickness_, angstrom);
  const double x0 = in_units(mv.x, angstrom);
  const double x1 = x0 + length_A * mv.cx;
  double inside = length_A;
  double a = x0;
  double b = x1;
  if (std::abs(x1 - x0) > 1e-12) {
    const double lo = std::max(std::min(x0, x1), 0.0);
    const double hi = std::min(std::max(x0, x1), th);
    inside = hi > lo ? length_A * (hi - lo) / std::abs(x1 - x0) : 0.0;
    a = std::clamp(x0, 0.0, th);
    b = std::clamp(x1, 0.0, th);
  } else if (x0 < 0.0 || x0 > th) {
    inside = 0.0;
  }
  if (!(inside > 0.0) || !(share > 0.0)) return 0.0 * electron_volt;
  const double mean = in_units(stopping_at(d.stopping, mv.energy), electron_volt_per_angstrom) * inside * share;
  double de = mean;
  if (cfg_.electronic_straggling && mean > 0.0) {
    const double le = std::log(std::max(in_units(mv.energy, electron_volt), 1.0));
    double fi = (le - d.strag_log_e0) * d.strag_inv_step;
    fi = std::clamp(fi, 0.0, static_cast<double>(d.strag.size() - 1) - 1e-9);
    const auto i = static_cast<std::size_t>(fi);
    const double w = fi - static_cast<double>(i);
    const double variance = (d.strag[i] + w * (d.strag[i + 1] - d.strag[i])) * inside * share;
    if (variance > 0.0) de = mean * gamma_deviate(mean * mean / variance, rng);
  }
  de = std::min(de, in_units(mv.energy, electron_volt));
  if (mv.is_ion) {
    deposit(t_.ionization_by_ion, a, b, de);
    t_.total_ionization_ion += de;
  } else {
    deposit(t_.ionization_by_recoils, a, b, de);
    t_.total_ionization_recoils += de;
  }
  return de * electron_volt;
}

bool Transporter::crystal_strike(const Mover &mv, const CrystalTarget &target, const double recoil_dir[3],
                                 energy_quantity T, energy_quantity e1, bool may_replace,
                                 std::uint32_t parent_point) {
  // The struck atom: the damage bookkeeping of the refined transport
  // (Sec. V of the paper: displaced if T > Ed, leaving with T - Eb; the
  // projectile replaces it when it is of the same element and keeps less
  // than the cut-off energy). Returns whether the mover replaced the atom.
  const LayerModel &layer = m_.layers[target.layer];
  const ElementModel &el = layer.elements[target.element];
  const bool full = mode_ != DamageCalculation::kinchin_pease;
  // The depth the event is booked at (a surface atom's thermal displacement
  // can put it a little outside the surface plane).
  const double txa = std::clamp(target.pos[0], 0.0, in_units(thickness_, angstrom));
  const double T_eV = in_units(T, electron_volt);
  double kp_damage = 0.0;
  bool displaces = false;
  if (full) {
    displaces = T > el.ed;
  } else {
    double kp_t = T_eV;
    if (cfg_.kp_subtract_binding_energy) kp_t -= in_units(el.eb, electron_volt);
    kp_damage = kp_damage_energy(el, std::max(kp_t, 0.0));
    displaces = T > 0.0 * electron_volt && kp_damage >= in_units(el.ed, electron_volt);
  }
  if (mv.is_ion && displaces) {
    add(t_.energy_to_recoils, txa, T_eV);
    add(t_.energy_absorbed[el.tally_index], txa, T_eV);
    t_.total_energy_to_recoils += T_eV;
  }
  // The struck atom, if it will be followed: built before the bookkeeping,
  // with its lineage and clock, so that MD events are settled first (in the
  // MD handoff mode MD takes it and the BCA books nothing for it).
  const bool surface_atom = !displaces && full && cfg_.surface_atom_depth > 0.0 * angstrom &&
                            txa * angstrom < cfg_.surface_atom_depth && T > el.es && el.es > 0.0 * electron_volt;
  const bool follows = full && ((displaces && T - el.eb > 0.0 * electron_volt &&
                                 (cfg_.max_recoils_per_ion == 0 || recoils_this_ion_ < cfg_.max_recoils_per_ion)) ||
                                surface_atom);
  Mover r;
  if (follows) {
    r.species = el.species;
    r.z = el.z;
    r.x = target.pos[0] * angstrom;
    r.y = target.pos[1] * angstrom;
    r.zz = target.pos[2] * angstrom;
    r.cx = recoil_dir[0];
    r.cy = recoil_dir[1];
    r.cz = recoil_dir[2];
    r.energy = surface_atom ? T : T - el.eb;
    r.is_ion = false;
    r.surface_only = surface_atom;
    r.origin_x = txa;
    r.generation = mv.generation + 1;
    r.tally_element = el.tally_index;
    r.track = ++next_track_;
    r.parent_track = mv.track;
    r.parent_point = parent_point;
    r.time = mv.time;
    r.md_mask = mv.md_mask;
    if (md_on_ && md_birth(r, T_eV, el)) {
      // MD takes it (handoff mode): it has left its site; nothing is booked.
      if (target.on_site) vacated_.insert(target.id);
      return false;
    }
  }
  bool replaced = false;
  if (displaces) {
    add(t_.displacements, txa, 1.0);
    t_.total_displacements += 1.0;
    if (mv.is_ion) add(t_.knock_ons, txa, 1.0);
    {
      const double er = std::max(T_eV - in_units(el.eb, electron_volt), 1.0);
      const long bin = spectrum_bin(er);
      if (bin >= 0 && bin < static_cast<long>(t_.recoil_spectrum.size())) {
        if (mv.is_ion) t_.pka_spectrum[static_cast<std::size_t>(bin)] += 1.0;
        t_.recoil_spectrum[static_cast<std::size_t>(bin)] += 1.0;
      }
    }
    if (target.on_site) vacated_.insert(target.id);
    if (!full) {
      // Quick calculation: modified Kinchin-Pease for the primary recoil.
      const double ed = in_units(el.ed, electron_volt);
      double nu = 1.0;
      if (kp_damage > 2.5 * ed) nu = 0.8 * kp_damage / (2.0 * ed);
      add(t_.vacancies_by_ion, txa, 1.0);
      add(t_.vacancies_by_recoils, txa, nu - 1.0);
      add(t_.vacancies[el.tally_index], txa, nu);
      t_.total_vacancies += nu;
      add(t_.ionization_by_recoils, txa, T_eV - kp_damage);
      t_.total_ionization_recoils += T_eV - kp_damage;
      add(t_.phonons_by_recoils, txa, kp_damage);
      t_.total_phonons_recoils += kp_damage;
    } else {
      const energy_quantity e2 = T - el.eb;
      t_.total_binding_energy_loss += in_units(el.eb, electron_volt);
      // Sec. V: the projectile replaces the target when it keeps less than
      // the lesser of its cut-off energy and Ed (and is of the same
      // element). With TRIM's 1 eV cut-off that is TRIM's rule; with the
      // cut-off at Ed it is Kinchin and Pease's.
      const energy_quantity cutoff = std::min(mv.is_ion ? ion_stop_ : recoil_stop_, el.ed);
      if (may_replace && e1 < cutoff && mv.z == el.z) {
        // Replacement collision: the mover takes the vacated site.
        // (The site stays out of the search for the rest of this history:
        // the atom that left starts from it, and the replacing atom is at
        // rest.)
        add(t_.replacements, txa, 1.0);
        t_.total_replacements += 1.0;
        replaced = true;
      } else {
        add(t_.vacancies[el.tally_index], txa, 1.0);
        t_.total_vacancies += 1.0;
      }
      if (follows) {
        stack_.push_back(r);
        ++recoils_this_ion_;
        t_.recoils_followed += 1.0;
      } else if (e2 > 0.0 * electron_volt) {
        // Recoil limit reached: its energy stays here as phonons.
        t_.recoil_limit_hits += 1.0;
        add(t_.phonons_by_recoils, txa, in_units(e2, electron_volt));
        t_.total_phonons_recoils += in_units(e2, electron_volt);
        add(t_.recoil_stopped[el.tally_index], txa, 1.0);
        add(t_.interstitials, txa, 1.0);
        t_.total_interstitials += 1.0;
      }
    }
  } else if (full && cfg_.surface_atom_depth > 0.0 * angstrom && txa * angstrom < cfg_.surface_atom_depth &&
             T > el.es && el.es > 0.0 * electron_volt) {
    // A surface atom set in motion below Ed (physics.surface_atom_depth,
    // off by default): followed (it may be sputtered), no defect if it
    // stays in the target.
    if (mv.is_ion) {
      add(t_.energy_to_recoils, txa, T_eV);
      add(t_.energy_absorbed[el.tally_index], txa, T_eV);
      t_.total_energy_to_recoils += T_eV;
    }
    if (target.on_site) vacated_.insert(target.id); // it left its site (and must not meet itself there)
    stack_.push_back(r);
  } else {
    // Sub-threshold transfer: the struck atom vibrates back (phonons).
    if (mv.is_ion) {
      add(t_.phonons_by_ion, txa, T_eV);
      t_.total_phonons_ion += T_eV;
    } else {
      add(t_.phonons_by_recoils, txa, T_eV);
      t_.total_phonons_recoils += T_eV;
    }
  }
  return replaced;
}

void Transporter::crystal_follow(Mover mv, LecuyerRandom &rng, std::size_t ion_index,
                                 std::vector<CollisionPoint> *trajectory) {
  const energy_quantity stop_energy = mv.is_ion ? ion_stop_ : recoil_stop_;
  const double th = in_units(thickness_, angstrom);
  const double simultaneous_A = in_units(cfg_.crystal.simultaneous_distance, angstrom);
  const double local_fraction = cfg_.crystal.local_stopping_fraction;
  const double decay = cfg_.crystal.local_stopping_decay;
  std::vector<crystal::Candidate> candidates;
  std::vector<CrystalTarget> set;        // the atoms of the current collision
  // Atoms already struck that are still ahead of the mover: not struck
  // again (the paper's xi_min; here by the sites' identities).
  std::vector<SiteKey> previous;
  std::vector<SiteKey> ahead;
  struct Outcome {
    physics::CollisionKinematics k;
    double shat[3];  // unit vector from the path to the atom
    double T_eV;
    double Q_eV;     // electronic loss in the collision
    double turn_A;   // xi - x1: path length to the deflection point
  };
  std::vector<Outcome> outcomes;
  record_point(trajectory, mv, CollisionPoint::Kind::birth);
  // The mover's last collision point, while it has not moved on (a stop
  // there marks it instead of adding a point at the same place).
  std::uint32_t at_point = CollisionPoint::no_parent;
  // The path ahead of the mover whose electronic loss is already charged,
  // in A: a flight is charged to the struck atoms' closest approach (xi),
  // but the mover turns x1 = xi - turn before it (Eq. 11), and the next
  // flight starts from there. Charging that part again would overstate the
  // loss by x1/xi per collision; for heavy slow ions, whose x1 is a few per
  // cent of a monolayer, the amorphous layers' ranges came out 2-5 % short.
  double charged_A = 0.0;

  for (std::size_t step = 0;; ++step) {
    if (step > 200000000) throw std::runtime_error("crystal transport: a moving atom took too many steps");
    if (mv.energy <= stop_energy) {
      stop_mover(mv, true, ion_index);
      if (at_point != CollisionPoint::no_parent) end_point(trajectory, at_point, mv);
      else record_point(trajectory, mv, CollisionPoint::Kind::stop);
      return;
    }
    const energy_quantity e_step = mv.energy; // (for the clock and the MD events)
    const double x = in_units(mv.x, angstrom);
    const double y = in_units(mv.y, angstrom);
    const double z = in_units(mv.zz, angstrom);
    if (!(mv.energy >= 0.0 * electron_volt) || !std::isfinite(x)) {
      throw std::runtime_error("a moving atom's energy or depth is not a number");
    }
    // Outside the target and moving away from it?
    if (x < 0.0 && mv.cx <= 0.0) {
      leave_front(mv, ion_index);
      record_point(trajectory, mv, CollisionPoint::Kind::exit);
      return;
    }
    if (x > th && mv.cx >= 0.0) {
      leave_back(mv, ion_index);
      record_point(trajectory, mv, CollisionPoint::Kind::exit);
      return;
    }
    const std::size_t li = m_.layer_index(std::clamp(x, 0.0, th) * angstrom);
    const LayerModel &layer = m_.layers[li];
    const SpeciesLayer &d = m_.data[mv.species][li];
    const double top = in_units(layer.top, angstrom);
    const double bottom = in_units(layer.bottom, angstrom);
    const double dir[3] = {mv.cx, mv.cy, mv.cz};
    if (mv.is_ion && !layer.grains) ion_last_grain_ = -1;

    // --- The next collision: its atoms, or a free step --------------------
    set.clear();
    ahead.clear();
    double free_step = 0.0;
    double share = 1.0;
    // The lattice here: the layer's, in a polycrystal that of the grain the
    // mover is in, or none (amorphous layers, and the slab at the grain
    // boundaries of a polycrystal when it has a width).
    const crystal::LatticeSearch *lattice = layer.lattice.get();
    crystal::Vec3 origin;            // of the lattice, in the layer's frame
    std::uint32_t domain = layer.domain;
    std::uint64_t configuration = history_seed_ + 0x9E3779B97F4A7C15ull * (li + 1);
    int grain = -1;
    double to_boundary = std::numeric_limits<double>::infinity();
    if (layer.grains) {
      const double look = lattice->reach_A() + lattice->max_impact_A() + lattice->max_displacement_A() +
                          layer.boundary_half_A + 1.0;
      const auto here = layer.grains->at({x - top, y + spot_y_, z + spot_z_}, look);
      grain = here.grain;
      to_boundary = here.boundary_A;
      if (mv.is_ion) {
        // (Per-grain results: where the ion entered and where it is.)
        ion_last_grain_ = grain;
        ion_last_layer_ = static_cast<std::uint32_t>(li);
        if (ion_entry_grain_ < 0 && grain >= 0) {
          ion_entry_grain_ = grain;
          ion_entry_layer_ = static_cast<std::uint32_t>(li);
        }
      }
      if (grain < 0 || to_boundary < layer.boundary_half_A) {
        lattice = nullptr;
      } else {
        lattice = &layer.grain_lattices[static_cast<std::size_t>(grain)];
        origin = layer.grains->grain(grain).origin_A;
        domain += static_cast<std::uint32_t>(grain);
        configuration += 0xD1B54A32D192ED03ull * (static_cast<std::uint64_t>(grain) + 1);
      }
    }
    const bool crystalline = lattice != nullptr;
    if (crystalline) {
      // Sec. IV: the sites ahead (xi > 0) within the maximum impact
      // parameter, the first along the path and those within the
      // simultaneity distance of it. The ion's entry point shifts the
      // lattice (spot_y_, spot_z_); each layer (and grain) has its own
      // configuration.
      share = crystal_path_share(mv.species, li);
      lattice->candidates({x - top - origin.x, y + spot_y_ - origin.y, z + spot_z_ - origin.z}, {mv.cx, mv.cy, mv.cz},
                          configuration, -origin.x, bottom - top - origin.x, candidates);
      // Near a grain boundary only the sites of this grain exist (those of
      // the next grain are met from inside it), and none in the boundary slab.
      const bool near_boundary = layer.grains && to_boundary < lattice->reach_A() + lattice->max_impact_A() +
                                                                   lattice->max_displacement_A() + layer.boundary_half_A;
      double first = 0.0;
      for (const auto &c : candidates) {
        const SiteKey key{domain, c.id};
        if (std::find(previous.begin(), previous.end(), key) != previous.end()) {
          ahead.push_back(key);
          continue;
        }
        if (!vacated_.empty() && vacated_.count(key) != 0) continue;
        if (set.empty()) first = c.xi;
        if (c.xi - first > simultaneous_A) break; // (sorted by xi)
        if (near_boundary) {
          const auto there = layer.grains->at({c.position.x + origin.x, c.position.y + origin.y, c.position.z + origin.z},
                                              layer.boundary_half_A + 1e-6);
          if (there.grain != grain || there.boundary_A < layer.boundary_half_A) continue;
        }
        CrystalTarget t;
        t.layer = li;
        t.element = layer.element_of(c.z);
        t.pos[0] = c.position.x + origin.x + top;
        t.pos[1] = c.position.y + origin.y - spot_y_;
        t.pos[2] = c.position.z + origin.z - spot_z_;
        t.xi = c.xi;
        t.impact = c.impact;
        t.on_site = true;
        t.id = key;
        set.push_back(t);
      }
      if (set.empty()) {
        // No atom within reach: move on, but not past the layer (its
        // neighbour's atoms are looked for from inside it), nor far past
        // the grain's boundary.
        free_step = lattice->reach_A();
        if (x >= top && x <= bottom) {
          if (mv.cx > 1e-12) free_step = std::min(free_step, (bottom - x) / mv.cx + 1e-6);
          else if (mv.cx < -1e-12) free_step = std::min(free_step, (x - top) / (-mv.cx) + 1e-6);
        }
        if (layer.grains) free_step = std::min(free_step, std::max(to_boundary - layer.boundary_half_A, 0.0) + 0.05);
      }
    } else {
      // Amorphous layer: one atom per cylinder of one monolayer length (the
      // monolayer collision step of the TRIM family), at a random impact
      // parameter and azimuth. The cylinder starts at the last struck atom's
      // closest approach, which lies charged_A ahead of the deflection point
      // (x1, hard-core approximation): one atom per monolayer of path, not
      // per monolayer less x1.
      const double length = in_units(layer.monolayer, angstrom) + charged_A;
      const double impact = in_units(layer.t13.pmax0, angstrom) * std::sqrt(rng());
      double cphi = 1.0;
      double sphi = 0.0;
      azimuth(rng, cphi, sphi);
      double px = 0.0;
      double py = 0.0;
      double pz = 0.0;
      perpendicular(mv.cx, mv.cy, mv.cz, cphi, sphi, px, py, pz);
      CrystalTarget t;
      t.pos[0] = x + length * mv.cx + impact * px;
      t.pos[1] = y + length * mv.cy + impact * py;
      t.pos[2] = z + length * mv.cz + impact * pz;
      if (t.pos[0] < 0.0 || t.pos[0] > th) {
        free_step = length; // no atom outside the target
      } else {
        t.layer = m_.layer_index(t.pos[0] * angstrom);
        const LayerModel &tl = m_.layers[t.layer];
        if (tl.elements.size() > 1) {
          const double r = rng();
          while (t.element + 1 < tl.elements.size() && r > tl.elements[t.element].cumulative) ++t.element;
        }
        t.xi = length;
        t.impact = impact;
        set.push_back(t);
      }
    }

    previous = ahead;
    if (set.empty()) {
      mv.energy -= crystal_path_loss(mv, d, std::max(free_step - charged_A, 0.0), share, rng);
      charged_A = std::max(charged_A - free_step, 0.0);
      mv.x += free_step * mv.cx * angstrom;
      mv.y += free_step * mv.cy * angstrom;
      mv.zz += free_step * mv.cz * angstrom;
      advance_time(mv, free_step * angstrom, e_step);
      at_point = CollisionPoint::no_parent;
      // MD events: the step's electronic loss may take it below an E_c (in a
      // channel most of the path is such steps). Without this check the next
      // step starts below E_c and that level never sees the crossing.
      if (md_on_ && md_crossing(mv, in_units(e_step, electron_volt))) {
        record_point(trajectory, mv); // handoff: where MD takes it
        if (mv.is_ion) break;
        return;
      }
      continue;
    }

    // --- Flight to the collision -------------------------------------------
    double xi_mean = 0.0;
    for (const auto &t : set) xi_mean += t.xi;
    xi_mean /= static_cast<double>(set.size());
    mv.energy -= crystal_path_loss(mv, d, std::max(xi_mean - charged_A, 0.0), share, rng);
    const energy_quantity e_c = mv.energy;
    if (e_c <= stop_energy) {
      mv.x += xi_mean * mv.cx * angstrom;
      mv.y += xi_mean * mv.cy * angstrom;
      mv.zz += xi_mean * mv.cz * angstrom;
      advance_time(mv, xi_mean * angstrom, e_step);
      at_point = CollisionPoint::no_parent;
      if (md_on_ && md_crossing(mv, in_units(e_step, electron_volt))) { // (as for a free step)
        record_point(trajectory, mv);
        if (mv.is_ion) break;
        return;
      }
      continue; // stops at the top of the loop
    }

    // --- Each collision as if alone (Sec. II) --------------------------------
    const double e_c_eV = in_units(e_c, electron_volt);
    const double s_atom = local_fraction > 0.0 && crystalline
                              ? in_units(stopping_at(d.stopping, e_c), electron_volt_per_angstrom) /
                                    in_units(layer.density, per_cubic_angstrom)
                              : 0.0; // eV A^2 per atom (the layer's mean)
    outcomes.clear();
    double sum_T = 0.0;
    double sum_Q = 0.0;
    for (const auto &t : set) {
      const CollisionPair &pair = m_.data[mv.species][t.layer].pairs[t.element];
      const double a = in_units(pair.screening_length, angstrom);
      const double eps = pair.reduced_energy(e_c);
      const double b = t.impact / a;
      Outcome o;
      o.k = scatter(pair, e_c, eps, b);
      o.T_eV = std::min(in_units(o.k.transferred, electron_volt), e_c_eV);
      o.Q_eV = 0.0;
      if (s_atom > 0.0) {
        // Electronic loss in the collision: exp(-c R/a) of the apsis R,
        // normalised so that straight paths at all impact parameters lose
        // S_e (the form of Oen and Robinson's local loss, with AMCSET's
        // stopping power in place of their k sqrt(E)).
        const double apsis = physics::closest_approach(pair.screening, eps, b) * a;
        o.Q_eV = local_fraction * s_atom * decay * decay / (2.0 * pi_v * a * a) * std::exp(-decay * apsis / a);
      }
      // The unit vector from the path to the atom.
      double perp[3];
      double n2 = 0.0;
      for (int c = 0; c < 3; ++c) {
        const double rel = t.pos[c] - (c == 0 ? x : c == 1 ? y : z);
        perp[c] = rel - t.xi * dir[c];
        n2 += perp[c] * perp[c];
      }
      if (n2 > 1e-24) {
        const double inv = 1.0 / std::sqrt(n2);
        for (int c = 0; c < 3; ++c) o.shat[c] = perp[c] * inv;
      } else {
        perpendicular(mv.cx, mv.cy, mv.cz, 1.0, 0.0, o.shat[0], o.shat[1], o.shat[2]); // head-on: any side
      }
      // The deflection point lies x1 before the closest approach on the
      // incoming line (Eq. 11). With the hard-core value of the time
      // integral, tau = s tan(theta/2), x1 = s tan(theta/2) and the struck
      // atom starts from its own position (x2 = 0); MARLOWE integrates tau.
      const double s2 = std::clamp(o.k.sin2_half_theta, 0.0, 1.0);
      const double x1 = s2 < 1.0 ? t.impact * std::sqrt(s2 / (1.0 - s2)) : 0.0;
      o.turn_A = t.xi - x1;
      sum_T += o.T_eV;
      sum_Q += o.Q_eV;
      outcomes.push_back(o);
    }
    // Simultaneous collisions are computed as if each occurred alone; if
    // together they take more than the projectile has, all transfers and
    // losses are scaled so that its energy just vanishes (Sec. IV).
    if (sum_T + sum_Q > e_c_eV) {
      const double scale = e_c_eV / (sum_T + sum_Q);
      for (auto &o : outcomes) {
        o.T_eV *= scale;
        o.Q_eV *= scale;
      }
      sum_T *= scale;
      sum_Q *= scale;
    }
    const energy_quantity e1 = std::max(e_c_eV - sum_T - sum_Q, 0.0) * electron_volt;

    // The momentum changes add (in units of sqrt(2 M)): p = p0 + sum (p_i - p0).
    const double p0 = std::sqrt(e_c_eV);
    double p[3] = {p0 * dir[0], p0 * dir[1], p0 * dir[2]};
    double turn = 0.0;
    for (const auto &o : outcomes) {
      const double pi_ = std::sqrt(std::max(e_c_eV - o.T_eV, 0.0));
      for (int c = 0; c < 3; ++c) {
        p[c] += pi_ * (o.k.cos_projectile * dir[c] - o.k.sin_projectile * o.shat[c]) - p0 * dir[c];
      }
      turn += o.turn_A;
    }
    turn = std::max(0.0, turn / static_cast<double>(outcomes.size()));
    charged_A = std::max(xi_mean - turn, 0.0);

    // --- Move to the deflection point; the electronic loss of the collision --
    mv.x += turn * mv.cx * angstrom;
    mv.y += turn * mv.cy * angstrom;
    mv.zz += turn * mv.cz * angstrom;
    advance_time(mv, turn * angstrom, e_step);
    if (sum_Q > 0.0) {
      const double xq = std::clamp(in_units(mv.x, angstrom), 0.0, th);
      if (mv.is_ion) {
        add(t_.ionization_by_ion, xq, sum_Q);
        t_.total_ionization_ion += sum_Q;
      } else {
        add(t_.ionization_by_recoils, xq, sum_Q);
        t_.total_ionization_recoils += sum_Q;
      }
    }

    // --- The mover's point: recorded before the struck atoms set off, so that
    // their points can name it as their parent ------------------------------
    double ndir[3] = {mv.cx, mv.cy, mv.cz};
    {
      const double pn = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
      if (pn > 1e-12 * p0) {
        for (int c = 0; c < 3; ++c) ndir[c] = p[c] / pn;
      }
    }
    std::uint32_t collision_point = CollisionPoint::no_parent;
    if (trajectory != nullptr) {
      Mover after = mv;
      after.cx = ndir[0];
      after.cy = ndir[1];
      after.cz = ndir[2];
      after.energy = e1;
      record_point(trajectory, after);
      Mover before = mv;
      before.energy = e_step;
      annotate_collision(trajectory, before, after, sum_T * electron_volt, e_step - e_c + sum_Q * electron_volt);
      collision_point = last_point_of(trajectory, after);
    }

    // --- The struck atoms ------------------------------------------------------
    // Only the atom that takes the most energy can be replaced by the mover.
    std::size_t hardest = 0;
    for (std::size_t i = 1; i < outcomes.size(); ++i) {
      if (outcomes[i].T_eV > outcomes[hardest].T_eV) hardest = i;
    }
    bool replaced = false;
    for (std::size_t i = 0; i < set.size(); ++i) {
      const Outcome &o = outcomes[i];
      // The struck atom leaves at the recoil angle on its own side of the path.
      const double rdir[3] = {o.k.cos_recoil * dir[0] + o.k.sin_recoil * o.shat[0],
                              o.k.cos_recoil * dir[1] + o.k.sin_recoil * o.shat[1],
                              o.k.cos_recoil * dir[2] + o.k.sin_recoil * o.shat[2]};
      if (crystal_strike(mv, set[i], rdir, o.T_eV * electron_volt, e1, i == hardest, collision_point)) {
        replaced = true;
        mv.x = set[i].pos[0] * angstrom;
        mv.y = set[i].pos[1] * angstrom;
        mv.zz = set[i].pos[2] * angstrom;
      }
      if (set[i].on_site) previous.push_back(set[i].id);
    }

    // --- The mover after the collision -------------------------------------------
    mv.cx = ndir[0];
    mv.cy = ndir[1];
    mv.cz = ndir[2];
    mv.energy = e1;
    at_point = collision_point;
    // MD events: its energy may have fallen below an E_c here (the handoff
    // mode then leaves it to MD; the cascades it started are still followed).
    if (md_on_ && md_crossing(mv, in_units(e_step, electron_volt))) {
      if (mv.is_ion) break;
      return;
    }
    if (replaced) {
      end_point(trajectory, collision_point, mv);
      stop_mover(mv, false, ion_index);
      if (mv.is_ion) break;
      return;
    }
    if (mv.is_ion) {
      // The cascades of the atoms this collision set in motion, now (with
      // the recoil stream: the ion's own path does not depend on them).
      while (!stack_.empty()) {
        const Mover r = stack_.back();
        stack_.pop_back();
        crystal_follow(r, *recoil_rng_, ion_index, trajectory);
      }
    }
  }
  // (An ion that replaced a target atom: its last cascade.)
  while (!stack_.empty()) {
    const Mover r = stack_.back();
    stack_.pop_back();
    crystal_follow(r, *recoil_rng_, ion_index, trajectory);
  }
}

void Transporter::run_ion_crystal(std::size_t ion_index, std::vector<CollisionPoint> *trajectory) {
  const auto &ion = cfg_.ion;
  stack_.clear();
  vacated_.clear();
  recoils_this_ion_ = 0;
  cascades_this_ion_ = 0;
  details_ = false; // no COLLISON.txt records in this transport
  LecuyerRandom ion_rng = make_stream(cfg_.seed, ion_index, 0);
  LecuyerRandom recoil_rng = make_stream(cfg_.seed, ion_index, 1);
  recoil_rng_ = &recoil_rng;
  // One frozen configuration of the lattices (occupancy, thermal
  // displacements) per ion history.
  history_seed_ = crystal::site_hash({}, (static_cast<std::uint64_t>(cfg_.seed) << 32) ^ ion_index);
  if (cfg_.source != TrimConfig::Source::beam) {
    spot_y_ = 0.0;
    spot_z_ = 0.0;
    run_source_crystal(ion_index, ion_rng, trajectory);
    recoil_rng_ = nullptr;
    t_.ions += 1;
    return;
  }

  // The entry point on the surface (BeamSpot): it shifts the lattice, so
  // positions in the results stay relative to it.
  spot_y_ = 0.0;
  spot_z_ = 0.0;
  const double w = in_units(ion.spot.width, angstrom);
  const double h = in_units(ion.spot.height, angstrom);
  switch (ion.spot.shape) {
  case BeamSpot::Shape::point:
    break;
  case BeamSpot::Shape::rectangle:
    spot_y_ = (ion_rng() - 0.5) * w;
    spot_z_ = (ion_rng() - 0.5) * h;
    break;
  case BeamSpot::Shape::circle: {
    const double r = 0.5 * w * std::sqrt(ion_rng());
    const double phi = 2.0 * pi_v * ion_rng();
    spot_y_ = r * std::cos(phi);
    spot_z_ = r * std::sin(phi);
    break;
  }
  case BeamSpot::Shape::gaussian: {
    const double to_sigma = 1.0 / (2.0 * std::sqrt(2.0 * std::log(2.0))); // FWHM -> sigma
    spot_y_ = w * to_sigma * ion_rng.gaussian();
    spot_z_ = h * to_sigma * ion_rng.gaussian();
    break;
  }
  }

  Mover mv;
  mv.species = 0;
  mv.z = ion.z;
  mv.cx = std::cos(in_units(ion.angle, radian));
  mv.cy = std::sin(in_units(ion.angle, radian));
  mv.cz = 0.0;
  // Start in front of the surface, on the line through the entry point, far
  // enough for every surface atom within the maximum impact parameter to lie
  // ahead (no energy is lost outside the target).
  double standoff = 0.0;
  if (const auto &lattice = m_.layers.front().lattice) {
    standoff = (lattice->max_impact_A() + lattice->max_displacement_A() + 0.1) / std::max(mv.cx, 1e-3);
  }
  mv.x = -standoff * mv.cx * angstrom;
  mv.y = -standoff * mv.cy * angstrom;
  mv.zz = 0.0 * angstrom;
  mv.energy = ion.energy;
  mv.is_ion = true;
  mv.generation = 0;
  mv.track = 0;
  next_track_ = 0;
  ion_entry_grain_ = -1;
  ion_last_grain_ = -1;
  ion_stopped_ = false;
  // MD events: an ion entering below an E_c is an event at its start (the
  // handoff mode then leaves the whole ion to MD).
  if (md_on_ && md_ion_entry(mv)) {
    record_point(trajectory, mv, CollisionPoint::Kind::birth);
  } else {
    crystal_follow(mv, ion_rng, ion_index, trajectory);
  }
  recoil_rng_ = nullptr;
  book_grains();
  t_.ions += 1;
}

void Transporter::book_grains() {
  // An ion's grains (TrimTallies::grains): its first grain, its final depth
  // if it stopped, and the grain it stopped in.
  if (ion_entry_grain_ >= 0) {
    auto &g = t_.grains[{ion_entry_layer_, static_cast<std::uint32_t>(ion_entry_grain_)}];
    g.entered += 1.0;
    if (ion_stopped_) {
      g.entered_stopped += 1.0;
      g.depth_sum_A += ion_stop_depth_A_;
      g.depth2_sum_A2 += ion_stop_depth_A_ * ion_stop_depth_A_;
    }
  }
  if (ion_stopped_ && ion_last_grain_ >= 0) {
    t_.grains[{ion_last_layer_, static_cast<std::uint32_t>(ion_last_grain_)}].stopped_here += 1.0;
  }
}

void Transporter::run_source_crystal(std::size_t ion_index, LecuyerRandom &ion_rng,
                                     std::vector<CollisionPoint> *trajectory) {
  // A particle of the source (TRIM.DAT, or one recoil atom repeated): not
  // TRIM's algorithm, but its bookkeeping (run_recoil_event_2013): the event
  // counts as an ion that stops where the particle starts.
  const ResolvedParticle &rp = m_.sources[ion_index / cfg_.particle_repeat];
  const SourceParticle &sp = rp.particle;
  double dir[3] = {sp.cos_x, sp.cos_y, sp.cos_z};
  {
    const double n = std::sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
    for (double &c : dir) c /= n; // (validated non-zero)
  }
  if (cfg_.source == TrimConfig::Source::particle_ions) {
    Mover mv;
    mv.species = 0;
    mv.z = cfg_.ion.z;
    mv.x = sp.x;
    mv.y = sp.y;
    mv.zz = sp.z_position;
    mv.cx = dir[0];
    mv.cy = dir[1];
    mv.cz = dir[2];
    mv.energy = sp.energy;
    mv.is_ion = true;
    mv.generation = 0;
    mv.track = 0;
    next_track_ = 0;
    if (md_on_ && md_ion_entry(mv)) {
      record_point(trajectory, mv, CollisionPoint::Kind::birth);
    } else {
      crystal_follow(mv, ion_rng, ion_index, trajectory);
    }
    return;
  }

  // A recoil atom: a target atom set in motion where it is. In a crystal
  // layer it is the nearest atom of its element on the lattice, which starts
  // from its own (thermally displaced) position and leaves its site empty;
  // elsewhere it starts at the given point.
  const LayerModel &own = m_.layers[rp.element_layer];
  const ElementModel &el = own.elements[rp.element];
  CrystalTarget target;
  target.layer = rp.element_layer;
  target.element = rp.element;
  target.pos[0] = in_units(sp.x, angstrom);
  target.pos[1] = in_units(sp.y, angstrom);
  target.pos[2] = in_units(sp.z_position, angstrom);
  const std::size_t li = rp.layer;
  const LayerModel &layer = m_.layers[li];
  const double top = in_units(layer.top, angstrom);
  const double bottom = in_units(layer.bottom, angstrom);
  const double x = target.pos[0];
  const double y = target.pos[1];
  const double z = target.pos[2];
  // The lattice there, as crystal_follow() finds it (the layer's, or the
  // grain's in a polycrystal; none in a grain-boundary slab).
  const crystal::LatticeSearch *lattice = layer.lattice.get();
  crystal::Vec3 origin;
  std::uint32_t domain = layer.domain;
  std::uint64_t configuration = history_seed_ + 0x9E3779B97F4A7C15ull * (li + 1);
  if (lattice != nullptr && layer.grains) {
    const double look = lattice->reach_A() + lattice->max_impact_A() + lattice->max_displacement_A() +
                        layer.boundary_half_A + 1.0;
    const auto here = layer.grains->at({x - top, y, z}, look);
    if (here.grain < 0 || here.boundary_A < layer.boundary_half_A) {
      lattice = nullptr;
    } else {
      lattice = &layer.grain_lattices[static_cast<std::size_t>(here.grain)];
      origin = layer.grains->grain(here.grain).origin_A;
      domain += static_cast<std::uint32_t>(here.grain);
      configuration += 0xD1B54A32D192ED03ull * (static_cast<std::uint64_t>(here.grain) + 1);
    }
  }
  bool has_element = false;
  for (const auto &e : layer.elements) has_element = has_element || e.z == el.z;
  if (lattice != nullptr && has_element) {
    crystal::SiteId id;
    crystal::Vec3 at;
    if (lattice->nearest_atom({x - top - origin.x, y - origin.y, z - origin.z}, configuration, el.z, -origin.x,
                              bottom - top - origin.x, id, at)) {
      target.layer = li;
      target.element = layer.element_of(el.z);
      target.pos[0] = at.x + origin.x + top;
      target.pos[1] = at.y + origin.y;
      target.pos[2] = at.z + origin.z;
      target.on_site = true;
      target.id = SiteKey{domain, id};
    }
  }
  const ElementModel &struck = m_.layers[target.layer].elements[target.element];
  const double xa = std::clamp(target.pos[0], 0.0, in_units(thickness_, angstrom));
  {
    Mover at_start;
    at_start.x = xa * angstrom;
    at_start.y = target.pos[1] * angstrom;
    at_start.zz = target.pos[2] * angstrom;
    at_start.is_ion = true;
    at_start.energy = 0.0 * electron_volt;
    stop_mover(at_start, false, ion_index);
  }
  const energy_quantity T = sp.energy;
  const double T_eV = in_units(T, electron_volt);
  if (mode_ != DamageCalculation::kinchin_pease && !(T > struck.ed)) {
    // Not displaced (the crystal transport takes TRIM's bugs fixed): its
    // energy becomes phonons where it is.
    add(t_.phonons_by_recoils, xa, T_eV);
    t_.total_phonons_recoils += T_eV;
    return;
  }
  // A virtual collision gives it the energy: its "projectile" is the event
  // (booked as the ion's primary knock-on), a root without a track.
  Mover event;
  event.x = target.pos[0] * angstrom;
  event.y = target.pos[1] * angstrom;
  event.zz = target.pos[2] * angstrom;
  event.cx = dir[0];
  event.cy = dir[1];
  event.cz = dir[2];
  event.is_ion = true;
  event.z = 0; // (never the struck atom's element: no replacement)
  event.generation = 0;
  event.track = CollisionPoint::no_parent;
  next_track_ = 0;
  crystal_strike(event, target, dir, T, 0.0 * electron_volt, false, CollisionPoint::no_parent);
  while (!stack_.empty()) {
    const Mover r = stack_.back();
    stack_.pop_back();
    crystal_follow(r, *recoil_rng_, ion_index, trajectory);
  }
}

// ---------------------------------------------------------------------------
// TRIM-2013 transport. Each block cites the TRIM.exe routine it reproduces
// (addresses in SRIM-2013.00's TRIM.exe; docs/TRIM2013_ALGORITHM.md).
// ---------------------------------------------------------------------------

std::size_t Transporter::pick_element(const LayerModel &layer, LecuyerRandom &rng) const {
  // 0x45f7c0: subtract the stoichiometric fractions from one random number
  // (drawn for pure elements too, which keeps TRIM's sequence of draws).
  double r = rng();
  if (layer.elements.size() == 1) return 0;
  for (std::size_t j = 0; j < layer.elements.size(); ++j) {
    r -= layer.elements[j].fraction;
    if (r < 0.0) return j;
  }
  return layer.elements.size() - 1;
}

energy_quantity Transporter::electronic_loss_2013(const Mover &mv, std::size_t li,
                                                  length_quantity L, LecuyerRandom &rng,
                                                  const common::stopping_power_quantity *known_se) const {
  // 0x45d4e0 SELOSS: S(E) L at the energy before the flight, plus a uniform
  // straggling deviate (R - 1/2) Omega sqrt(L), limited to +-S L.
  const SpeciesLayer &d = m_.data[mv.species][li];
  const energy_quantity e = std::max(mv.energy, 0.001 * electron_volt);
  const energy_quantity mean =
      (known_se != nullptr && e == mv.energy ? *known_se : stopping_at(d.stopping, e)) * L;
  energy_quantity de = mean;
  // TRIM draws the straggling deviate for every flight.
  const double r = rng();
  if (cfg_.electronic_straggling) {
    const auto &t = m_.layers[li].t13;
    const double e_eV = in_units(e, electron_volt);
    double omega = t.k_low * e_eV; // eV per sqrt(Angstrom)
    if (omega > t.k_bohr) {
      const double mc2_eV = in_units(m_.species[mv.species].mass, atomic_mass_unit) * 931494.0e3;
      const double g = 1.0 + e_eV / mc2_eV;
      const double beta2 = 1.0 - 1.0 / (g * g);
      const double inv = 1.0 / beta2;
      omega = t.k_bohr * (inv - 0.5) / (inv - 1.0);
    }
    double strag = (r - 0.5) * omega * std::sqrt(in_units(L, angstrom));
    const double m_eV = in_units(mean, electron_volt);
    if (std::abs(strag) > m_eV) strag = std::copysign(m_eV, strag);
    de = mean + strag * electron_volt;
  }
  if (de > e) de = e;
  return de;
}

void Transporter::edge_2013(const Mover &mv, std::size_t li, length_quantity &L,
                            LecuyerRandom &rng) {
  // 0x44b7a0: the flight would leave the depth window (or the ion is within
  // p_max of the front surface, moving out). Shorten it by random factors.
  const double x = in_units(mv.x, angstrom);
  const double l0 = in_units(m_.layers[li].monolayer, angstrom);
  const double binw = (w1_ - w0_) / static_cast<double>(nbins_);
  double D = 0.0;
  if (mv.cx >= 0.0) {
    D = w1_ - x;
    if (D < 0.0) D = in_units(m_.layers[li].bottom, angstrom) - x;
  } else {
    D = x - w0_;
  }
  const bool quick = mode_ == DamageCalculation::kinchin_pease; // TRIM codes 1, 4, 6
  double l = in_units(L, angstrom);
  while (std::abs(l * mv.cx) > D) {
    l *= rng();
    if (quick) {
      if (std::abs(l * mv.cx) < binw) {
        l = std::abs(binw / mv.cx);
        break;
      }
    } else if (std::abs(l * mv.cx) < l0) {
      l = l0;
      break;
    }
  }
  if (mv.cx > 0.0 && l * mv.cx > D) {
    // TRIM adds the layer's atomic density (atoms/A^3) to the distance here.
    l = (D + in_units(m_.layers[li].density, per_cubic_angstrom)) / mv.cx;
  }
  L = l * angstrom;
  flight_counter_ = 2;
}

void Transporter::interface_2013(const Mover &mv, std::size_t li, length_quantity &L,
                                 LecuyerRandom &rng) {
  // 0x44b330: the flight would cross into another layer.
  const double x = in_units(mv.x, angstrom);
  const double l0 = in_units(m_.layers[li].monolayer, angstrom);
  const double binw = (w1_ - w0_) / static_cast<double>(nbins_);
  double D = 0.0;
  if (mv.cx >= 0.0) {
    D = in_units(m_.layers[li].bottom, angstrom) - x;
  } else {
    D = li == 0 ? x : x - in_units(m_.layers[li - 1].bottom, angstrom);
  }
  if (D < l0) D = l0;
  double l = in_units(L, angstrom);
  while (std::abs(l * mv.cx) > D) {
    l *= rng();
    const double g = std::abs(l * mv.cx);
    if (binw > 3.0 * l0) {
      if (g < binw && x > 10.0 * l0) {
        l = std::min(binw, D / std::max(std::abs(mv.cx), 0.01));
        break;
      }
    } else if (g < l0) {
      l = l0;
      break;
    }
  }
  L = l * angstrom;
  flight_counter_ = 2;
}

length_quantity Transporter::ion_flight_path_2013(const Mover &mv, std::size_t li,
                                                  LecuyerRandom &rng) {
  // 0x44ab60.
  const LayerModel &first = m_.layers.front();
  if (mv.x < first.monolayer && mv.energy == cfg_.ion.energy) {
    flight_counter_ = 2;
    return first.monolayer;
  }
  const LayerModel &layer = m_.layers[li];
  if (mode_ == DamageCalculation::monolayer) return layer.monolayer;
  const auto &t = layer.t13;
  const double eps = (mv.energy * t.fi).value();
  // pmax = a / (e + sqrt(e) + 0.125 e^0.1) with e = sqrt(eps epsdg); with
  // fast_math the denominator comes from a table (no sqrt, pow or log).
  const length_quantity pmax = [&] {
    if (cfg_.fast_math) return t.screening / physics::fast::flight_path_denominator(eps * t.epsdg);
    const double eeg = std::sqrt(eps * t.epsdg);
    return t.screening / (eeg + std::sqrt(eeg) + 0.125 * std::pow(eeg, 0.1));
  }();
  length_quantity L = 1.0 / (pi_v * pmax * pmax * layer.density);
  const length_quantity xn = mv.x + L * mv.cx;
  if (m_.layer_index(xn) != li) interface_2013(mv, li, L, rng);
  const double xna = in_units(xn, angstrom);
  if (xna < w0_ || xna > w1_) edge_2013(mv, li, L, rng);
  if (mv.x < pmax && mv.cx < 0.0) edge_2013(mv, li, L, rng);
  if (mode_ == DamageCalculation::full_cascade && cfg_.source == TrimConfig::Source::beam) {
    // Full cascades (cascade code 2 only, not the TRIM.DAT code 5): every
    // ion flight is a random fraction (at least a monolayer) of the free
    // flight path.
    return std::max(layer.monolayer, rng() * L);
  }
  if (flight_counter_ > 0) {
    L = std::max(layer.monolayer, rng() * L);
    --flight_counter_;
  }
  return L;
}

void Transporter::displace_2013(const ElementModel &el, double x, bool by_ion, double T_eV) {
  add(t_.displacements, x, 1.0);
  t_.total_displacements += 1.0;
  add(t_.vacancies[el.tally_index], x, 1.0);
  t_.total_vacancies += 1.0;
  if (by_ion) add(t_.knock_ons, x, 1.0);
  const double er = std::max(T_eV - in_units(el.eb, electron_volt), 1.0);
  const long bin = spectrum_bin(er);
  if (bin >= 0 && bin < static_cast<long>(t_.recoil_spectrum.size())) {
    if (by_ion) t_.pka_spectrum[static_cast<std::size_t>(bin)] += 1.0;
    t_.recoil_spectrum[static_cast<std::size_t>(bin)] += 1.0;
  }
}

void Transporter::kp_2013(const ElementModel &el, const LayerModel &layer, double T_eV,
                          double x, std::size_t ion_index, const Mover &ion,
                          common::stopping_power_quantity se) {
  // 0x44bae0: modified Kinchin-Pease with TRIM's layer-averaged Robinson
  // partition (eps = 0.01 Z^-7/3 T, g = eps + 0.4 eps^3/4 + 3.4 eps^1/6).
  const auto &t = layer.t13;
  const double eps = t.kp_eps_per_eV * T_eV;
  double g = 0.0;
  if (cfg_.fast_math) {
    g = physics::fast::kp_partition(eps); // a table: no pow or cbrt
  } else {
    g = eps + 0.4 * std::pow(eps, 0.75) + 3.4 * std::pow(eps, 1.0 / 6.0);
  }
  const double ev = T_eV / (1.0 + t.kp_k * g);
  const double ed = in_units(el.ed, electron_volt);
  const double eb = in_units(el.eb, electron_volt);
  book_pka(el, x, ev); // damage models (own tallies only)
  add(t_.ionization_by_recoils, x, T_eV - ev);
  t_.total_ionization_recoils += T_eV - ev;
  if (ev < ed) {
    add(t_.phonons_by_recoils, x, ev);
    t_.total_phonons_recoils += ev;
    return;
  }
  const double nu = ev > 2.5 * ed ? 0.8 * ev / (2.0 * ed) : 1.0;
  add(t_.displacements, x, 1.0);
  t_.total_displacements += 1.0;
  add(t_.knock_ons, x, 1.0);
  add(t_.vacancies_by_ion, x, 1.0);
  add(t_.vacancies_by_recoils, x, nu - 1.0);
  add(t_.vacancies[el.tally_index], x, nu);
  t_.total_vacancies += nu;
  {
    const double er = std::max(T_eV - eb, 1.0);
    const long bin = spectrum_bin(er);
    if (bin >= 0 && bin < static_cast<long>(t_.pka_spectrum.size())) {
      t_.pka_spectrum[static_cast<std::size_t>(bin)] += 1.0;
      t_.recoil_spectrum[static_cast<std::size_t>(bin)] += 1.0;
    }
  }
  // TRIM books E2RECOIL ("from IONS") - Eb once for the displaced atom, and
  // phonons Ev - nu Eb (the binding energy is lost); the "absorbed" column
  // keeps the whole transfer.
  add(t_.energy_to_recoils, x, -eb);
  t_.total_energy_to_recoils -= eb;
  add(t_.phonons_by_recoils, x, ev - nu * eb);
  t_.total_phonons_recoils += ev - nu * eb;
  t_.total_binding_energy_loss += nu * eb;
  if (details_) {
    t_.cascades.push_back({ion_index + 1, in_units(ion.energy, kilo_electron_volt), x,
                           in_units(ion.y, angstrom), in_units(ion.zz, angstrom),
                           in_units(se, electron_volt_per_angstrom), el.z, T_eV, nu, nu, 0.0,
                           0.0, {}});
  }
}

void Transporter::recoil_stop_2013(std::size_t k) {
  // 0x457970 RECSTOPDATA.
  const Mover &r = slots_[k];
  const double x = in_units(r.x, angstrom);
  const double e = in_units(r.energy, electron_volt);
  if (track_primary_stop_ && r.track == 1) primary_stop_x_ = x;
  if (in_window(x)) {
    add(t_.phonons_by_recoils, x, e);
    t_.total_phonons_recoils += e;
  }
  if (!primary_displaced_) {
    // The primary slot was never displaced (sub-threshold ion collision).
    primary_displaced_ = true;
    return;
  }
  if (sub_threshold_) {
    sub_threshold_ = false;
    return;
  }
  bool replaced = false;
  if (k > 1) {
    // Replacement: the atom just below in the stack (the last one this mover
    // struck) was displaced and is the same element. TRIM still counts the
    // mover as an interstitial afterwards.
    const Mover &below = slots_[k - 1];
    const auto &bl = m_.layers[below.origin_layer].elements[below.element_hint];
    if (below.energy > bl.ed && below.species == r.species) {
      replaced = true;
      note_atom(r.z, 0.0, r, false);
      add(t_.vacancies[bl.tally_index], x, -1.0);
      add(t_.vacancies_removed[bl.tally_index], x, 1.0);
      t_.total_vacancies -= 1.0;
      add(t_.replacements, x, 1.0);
      t_.total_replacements += 1.0;
      const double eb = in_units(bl.eb, electron_volt);
      add(t_.phonons_by_recoils, x, eb);
      t_.total_phonons_recoils += eb;
    }
  }
  add(t_.recoil_stopped[r.tally_element], x, 1.0);
  if (replaced && !cfg_.bugs.replacement_interstitial) return; // fixed: it fills the vacancy
  add(t_.interstitials, x, 1.0);
  t_.total_interstitials += 1.0;
}

void Transporter::cascade_2013(std::size_t ion_index, std::vector<CollisionPoint> *trajectory) {
  // TRIMMain.Main RECOILLOOP / NEXTRECOIL. slots_[1] holds the primary; a
  // moving atom that strikes another is copied one slot up (and continues),
  // the struck atom takes its old slot. Recoils move one monolayer per step.
  LecuyerRandom &rng = *recoil_rng_;
  const energy_quantity stop = recoil_stop_;
  const double th = in_units(thickness_, angstrom);
  const LayerModel &layer1 = m_.layers.front();
  std::size_t k = 1;
  enum class Next { move, pop };
  Next next = Next::move;
  while (true) {
    if (next == Next::pop) {
      // NEXTRECOIL.
      if (--k == 0) return;
      sub_threshold_ = false;
      Mover &r = slots_[k];
      if (r.handed) continue; // MD handoff: MD takes this atom (no draw, no tallies)
      const auto &el = m_.layers[r.origin_layer].elements[r.element_hint];
      const double x = in_units(r.x, angstrom);
      // TRIM draws the random number for every atom it takes up (0x437f62).
      const double rn = rng();
      if (r.energy < el.ed && r.x > 9.0 * layer1.monolayer * rn) {
        // Below Ed and deeper than a random 0-9 monolayers: not followed.
        // (TRIM leaves the slot's energy unchanged.)
        if (r.energy != 0.0 * electron_volt && in_window(x)) {
          const double e = in_units(r.energy, electron_volt);
          add(t_.phonons_by_recoils, x, e);
          t_.total_phonons_recoils += e;
        }
        continue;
      }
      if (r.energy >= el.ed) {
        note_atom(el.z, in_units(r.energy, electron_volt), r, true);
        if (in_window(x)) displace_2013(el, x, false, in_units(r.energy, electron_volt));
        else {
          t_.total_displacements += 1.0;
          t_.total_vacancies += 1.0;
        }
        t_.recoils_followed += 1.0;
        if (r.x > 2.0 * layer1.monolayer) {
          r.energy -= el.eb;
          t_.total_binding_energy_loss += in_units(el.eb, electron_volt);
          // MD events: losing Eb can take it below an E_c.
          if (md_on_ && md_crossing(r, in_units(r.energy + el.eb, electron_volt))) continue;
        }
      } else {
        sub_threshold_ = true;
      }
      if (r.energy > stop) {
        record_point(trajectory, r, CollisionPoint::Kind::birth); // where the atom starts to move
        next = Next::move;
      } else {
        recoil_stop_2013(k);
        next = Next::pop;
      }
      continue;
    }

    // RECOILLOOP.
    Mover &r = slots_[k];
    // A NaN energy never falls below the stop energy: the loop would not
    // end. (One comparison per step; it cannot happen with valid input.)
    if (!(r.energy >= 0.0 * electron_volt)) {
      throw std::runtime_error("a recoil's energy is not a number (recoil slot " + std::to_string(k) + ")");
    }
    const LayerModel &lay = m_.layers[r.layer_hint];
    const length_quantity l0 = lay.monolayer;
    r.x += l0 * r.cx;
    r.y += l0 * r.cy;
    r.zz += l0 * r.cz;
    advance_time(r, l0, r.energy);
    record_point(trajectory, r);
    double x = in_units(r.x, angstrom);
    if (x > th) {
      t_.transmitted_recoils[r.tally_element] += 1.0;
      mark_last(trajectory, r, CollisionPoint::Kind::exit);
      next = Next::pop;
      continue;
    }
    if (x < 0.0) {
      // 0x441ff0 SPUTTER: planar barrier on the normal energy, with the
      // atom's own surface binding energy (TRIM keeps one Es per target
      // atom). layer_hint is the layer the atom moved through last, so its
      // element list need not contain this atom at element_hint.
      const auto &el = m_.layers[r.origin_layer].elements[r.element_hint];
      const energy_quantity normal = r.energy * r.cx * r.cx;
      surface_arrival(r);
      if (normal < el.es) {
        r.cx = -r.cx; // reflected back into the target
        r.x = 1e-4 * angstrom;
        x = 1e-4;
        // The recorded point shows where the atom goes on: at the surface,
        // reflected (no effect on the transport).
        if (trajectory != nullptr && !trajectory->empty() && !r.untracked && trajectory->back().track == r.track) {
          auto &pt = trajectory->back();
          pt.coordinate.z_ = r.x;
          pt.vz = -pt.vz;
        }
      } else {
        const energy_quantity out = std::max(r.energy - el.es, 0.0 * electron_volt);
        const double eo = in_units(out, electron_volt);
        t_.sputtered[r.tally_element] += 1.0;
        t_.sputtered_energy[r.tally_element] += eo;
        if (sub_threshold_ && cfg_.bugs.subthreshold_sputter_counted) {
          // TRIM counts a sputtered atom that was set in motion below Ed as
          // a displacement and a vacancy, in its totals only (VACANCY.txt's
          // header, COLLISON.txt), not in the depth tables.
          t_.total_displacements += 1.0;
          t_.total_vacancies += 1.0;
        }
        if (cfg_.record_exiting_atoms && t_.exits.size() < cfg_.max_exit_records) {
          // SPUTTER.txt as TRIM writes it (0x442875-0x4430c3): the energy
          // inside the surface (the binding energy is not subtracted) and the
          // direction cosines before the barrier, limited to +-0.9999999.
          const auto lim = [](double c) { return std::clamp(c, -0.9999999, 0.9999999); };
          t_.exits.push_back({'S', ion_index + 1, r.z, in_units(r.energy, electron_volt), 0.0,
                              in_units(r.y, angstrom), in_units(r.zz, angstrom), lim(r.cx),
                              lim(r.cy), lim(r.cz)});
        }
        mark_last(trajectory, r, CollisionPoint::Kind::exit);
        r.energy = 0.0 * electron_volt;
        next = Next::pop;
        continue;
      }
    }
    const std::size_t li = m_.layer_index(r.x);
    r.layer_hint = li;
    // Electronic loss over one monolayer at the current energy.
    energy_quantity de = electronic_loss_2013(r, li, m_.layers[li].monolayer, rng);
    if (de > r.energy) de = r.energy;
    {
      const double ev = in_units(de, electron_volt);
      if (in_window(x)) {
        add(t_.ionization_by_recoils, x, ev);
        t_.total_ionization_recoils += ev;
      }
    }
    if (r.energy < stop) {
      recoil_stop_2013(k);
      mark_last(trajectory, r, CollisionPoint::Kind::stop);
      next = Next::pop;
      continue;
    }
    // Collision with an atom of the layer, monolayer impact parameter.
    const LayerModel &tl = m_.layers[li];
    const std::size_t j = pick_element(tl, rng);
    const ElementModel &el = tl.elements[j];
    length_quantity p = std::sqrt(rng()) * tl.t13.pmax0;
    if (r.z < 21 && r.x < 2.0 * layer1.monolayer) {
      // Light atoms at the surface: smaller impact parameters.
      if (r.z < 14) {
        p *= rng();
      } else {
        const length_quantity p1 = std::sqrt(rng()) * tl.t13.pmax0;
        const length_quantity p2 = p1 * std::sqrt(rng());
        p = p2 + (p1 - p2) * (static_cast<double>(r.z) - 13.0) / 7.0;
      }
    }
    const CollisionPair &pair = m_.data[r.species][li].pairs[j];
    const double eps = pair.reduced_energy(r.energy);
    const auto kin = scatter(pair, r.energy, eps, (p / pair.screening_length).value());
    energy_quantity T = std::min(kin.transferred, r.energy);
    if (!cfg_.bugs.recoil_energy_created && T > r.energy - de) {
      // Fixed: the struck atom gets at most what the electronic loss left.
      T = std::max(r.energy - de, 0.0 * electron_volt);
    }
    double cphi = 1.0;
    double sphi = 0.0;
    azimuth(rng, cphi, sphi);

    if (k + 1 >= slots_.size()) slots_.resize(slots_.size() * 2 + 2);
    Mover &mover = slots_[k]; // (re-take after a possible resize)
    Mover struck = mover;
    // The mover continues one slot up with its new energy and direction.
    Mover moved = mover;
    // TRIM.exe limits the mover's energy to >= 0 (0x436de5); T + DEE can
    // exceed E for slow atoms, and TRIM keeps both T and DEE in its tallies.
    moved.energy = std::max(mover.energy - T - de, 0.0 * electron_volt);
    turn(moved.cx, moved.cy, moved.cz, kin.cos_projectile, kin.sin_projectile, -cphi, -sphi);
    // The struck atom takes the old slot with energy T. TRIM.exe turns its
    // direction with DIRCOS (0x45ce40) called with sin(x) = 0 (0x4377d5),
    // which selects DIRCOS's branch for motion along the x axis: only the
    // mover's old x cosine enters, the transverse direction is random.
    struck.species = el.species;
    struck.z = el.z;
    struck.tally_element = el.tally_index;
    struck.layer_hint = li;
    struck.element_hint = j;
    struck.origin_layer = li;
    struck.energy = T;
    struck.generation = mover.generation + 1;
    struck.track = ++next_track_;
    if (cfg_.bugs.struck_atom_direction) {
      if (cfg_.fast_math) {
        rotate_about_x_2013<true>(struck.cx, struck.cy, struck.cz, kin.cos_recoil, kin.sin_recoil, cphi, sphi);
      } else {
        rotate_about_x_2013(struck.cx, struck.cy, struck.cz, kin.cos_recoil, kin.sin_recoil, cphi, sphi);
      }
    } else {
      // Fixed: turned from the mover's direction before the collision, as
      // the primary of an ion collision is.
      turn(struck.cx, struck.cy, struck.cz, kin.cos_recoil, kin.sin_recoil, cphi, sphi);
    }
    moved.layer_hint = li;
    // Trajectory lineage and MD events (opt-in; without the handoff mode
    // the transport and the tallies are the same).
    struck.parent_track = mover.track;
    struck.parent_point = last_point_of(trajectory, mover);
    struck.handed = false;
    annotate_collision(trajectory, mover, moved, T, de);
    if (md_on_) {
      const double before_eV = in_units(mover.energy, electron_volt);
      struck.handed = md_birth(struck, in_units(T, electron_volt), el);
      moved.handed = md_crossing(moved, before_eV);
    }
    slots_[k] = struck;
    slots_[k + 1] = moved;
    ++k;
    if (cfg_.max_recoils_per_ion != 0 && ++recoils_this_ion_ > cfg_.max_recoils_per_ion) {
      t_.recoil_limit_hits += 1.0;
    }
    if (slots_[k].handed) {
      next = Next::pop; // MD handoff: the moving atom is MD's from this collision on
    } else if (slots_[k].energy >= stop) {
      next = Next::move;
    } else {
      recoil_stop_2013(k);
      mark_last(trajectory, slots_[k], CollisionPoint::Kind::stop);
      next = Next::pop;
    }
  }
}

void Transporter::surface_arrival(const Mover &r) {
  // 0x441ff0: histogram of the normal energy of every atom reaching the
  // surface, bin round(10 E cos^2 x) limited to 1..200 (VB rounds half to even).
  const double en = in_units(r.energy, electron_volt) * r.cx * r.cx;
  double b = std::nearbyint(en / TrimTallies::surface_energy_bin_eV);
  // (Not std::clamp: a NaN would pass it and index anywhere.)
  if (!(b >= 1.0)) b = 1.0;
  if (b > static_cast<double>(TrimTallies::surface_energy_bins)) b = static_cast<double>(TrimTallies::surface_energy_bins);
  auto &h = t_.surface_energy[r.tally_element];
  h[static_cast<std::size_t>(b) - 1] += 1.0;
}

void Transporter::record_exyz(std::size_t ion_index, const Mover &mv, double se_eV_per_A) {
  if (t_.exyz.size() >= cfg_.max_exit_records) return;
  const double e = in_units(mv.energy, electron_volt);
  t_.exyz.push_back({ion_index + 1, e < in_units(ion_stop_, electron_volt) ? 0.0 : e / 1000.0,
                     in_units(mv.x, angstrom), in_units(mv.y, angstrom), in_units(mv.zz, angstrom),
                     se_eV_per_A, last_T_eV_});
}

void Transporter::start_cascade_2013(const ElementModel &el, std::size_t li, std::size_t j,
                                     const Mover &site, energy_quantity T,
                                     std::size_t ion_index,
                                     std::vector<CollisionPoint> *trajectory, LecuyerRandom &rng,
                                     bool restore_rng, double ion_energy_keV,
                                     common::stopping_power_quantity se) {
  // STARTRECOIL: the primary takes slot 1 at the collision site with the
  // direction given in site.cx..cz.
  const double x1 = in_units(site.x, angstrom);
  const double T_eV = in_units(T, electron_volt);
  Mover &r = slots_[1];
  r.species = el.species;
  r.z = el.z;
  r.tally_element = el.tally_index;
  r.layer_hint = li;
  r.element_hint = j;
  r.origin_layer = li;
  r.x = site.x;
  r.y = site.y;
  r.zz = site.zz;
  r.cx = site.cx;
  r.cy = site.cy;
  r.cz = site.cz;
  r.is_ion = false;
  r.generation = 1;
  r.track = ++next_track_;
  // Lineage and clock from the collision (site.track and site.parent_point:
  // the mover's track and collision point; no_parent for a source recoil).
  r.parent_track = site.track;
  r.parent_point = site.parent_point;
  r.time = site.time;
  r.md_mask = site.md_mask;
  r.handed = false;
  r.untracked = false;
  sub_threshold_ = false;
  if (md_on_) {
    // TRIM's sub-threshold transfers (T <= Ed, run with a stale energy)
    // are never events. In the handoff mode MD takes the primary instead
    // of the BCA: no displacement, no cascade here.
    const bool handed = T > el.ed ? md_birth(r, T_eV, el) : md_below_min(r, T_eV);
    if (handed) return;
  }
  const double d0 = t_.total_displacements;
  const double v0 = t_.total_vacancies;
  const double r0 = t_.total_replacements;
  const double i0 = t_.total_interstitials;
  std::vector<CascadeAtom> atoms;
  cascade_atoms_ = details_ && cfg_.collision_recoil_details ? &atoms : nullptr;
  if (T > el.ed) {
    note_atom(el.z, T_eV, r, true);
    primary_displaced_ = true;
    if (in_window(x1)) {
      displace_2013(el, x1, true, T_eV);
    } else {
      t_.total_displacements += 1.0;
      t_.total_vacancies += 1.0;
    }
    t_.recoils_followed += 1.0;
    r.energy = std::max(T - el.eb, 0.0 * electron_volt);
    t_.total_binding_energy_loss += in_units(el.eb, electron_volt);
  } else {
    // TRIM leaves the slot's energy from the previous cascade.
    primary_displaced_ = false;
  }
  // An atom left with no energy to move (the stale energy of a first
  // cascade is 0) gets no trajectory points; the transport is unchanged.
  r.untracked = !(r.energy > recoil_stop_);
  record_point(trajectory, r, CollisionPoint::Kind::birth);
  if (restore_rng && !cfg_.bugs.cascade_random_replay) {
    // Fixed: every cascade draws from its own stream.
    LecuyerRandom cascade_rng =
        make_stream(cfg_.seed, ion_index, static_cast<unsigned>(2 + cascades_this_ion_++));
    recoil_rng_ = &cascade_rng;
    cascade_2013(ion_index, trajectory);
  } else if (restore_rng) {
    // TRIM.exe saves the random-number state before a cascade and restores
    // it afterwards for cascade codes 2 and 3 (0x434e89, 0x438140): the
    // cascade draws the ion's next random numbers, and the ion then draws
    // them again.
    LecuyerRandom cascade_rng = rng;
    recoil_rng_ = &cascade_rng;
    cascade_2013(ion_index, trajectory);
  } else {
    // TRIM.DAT codes 5 and 7: no restore, the ion continues after the
    // cascade's numbers.
    recoil_rng_ = &rng;
    cascade_2013(ion_index, trajectory);
  }
  recoil_rng_ = nullptr;
  cascade_atoms_ = nullptr;
  if (details_ && T > el.ed) {
    t_.cascades.push_back({ion_index + 1, ion_energy_keV, x1, in_units(site.y, angstrom),
                           in_units(site.zz, angstrom), in_units(se, electron_volt_per_angstrom),
                           el.z, T_eV, t_.total_displacements - d0, t_.total_vacancies - v0,
                           t_.total_replacements - r0, t_.total_interstitials - i0, std::move(atoms)});
  }
}

void Transporter::run_recoil_event_2013(std::size_t ion_index,
                                        std::vector<CollisionPoint> *trajectory) {
  // TRIM.DAT codes 6 and 7 (0x449330): an atom of the target is set in
  // motion at the given point and direction, without an incident ion. TRIM
  // goes straight to STARTRECOIL (no depth-window check, no DIRCOS) with the
  // atom's energy as the transferred energy.
  const ResolvedParticle &rp = m_.sources[ion_index / cfg_.particle_repeat];
  const SourceParticle &sp = rp.particle;
  LecuyerRandom rng = make_stream(cfg_.seed, ion_index, 0);
  next_track_ = 0;
  details_ = ion_index < cfg_.collision_details_ions;
  const IonStart start = ion_start();
  const LayerModel &layer = m_.layers[rp.element_layer];
  const ElementModel &el = layer.elements[rp.element];
  Mover site;
  site.x = sp.x;
  site.y = sp.y;
  site.zz = sp.z_position;
  site.cx = sp.cos_x;
  site.cy = sp.cos_y;
  site.cz = sp.cos_z;
  site.track = CollisionPoint::no_parent; // the cascade's primary is a root
  site.parent_point = CollisionPoint::no_parent;
  const double x = in_units(sp.x, angstrom);
  const energy_quantity T = sp.energy;
  const double T_eV = in_units(T, electron_volt);
  // TRIM books the event as an ion that stops where the atom starts (the
  // range statistics of RANGE.txt), and adds the depth where the atom
  // itself stops to the ion column of the range table as well.
  Mover at_start = site;
  at_start.is_ion = true;
  at_start.energy = 0.0 * electron_volt;
  stop_mover(at_start, false, ion_index);
  primary_stop_x_ = std::numeric_limits<double>::quiet_NaN();
  track_primary_stop_ = true;
  if (T >= recoil_stop_ && mode_ != DamageCalculation::kinchin_pease) {
    book_pka(el, x, layer_damage_energy(m_.layers[rp.layer], T_eV)); // damage models (own tallies only)
  }
  if (T < recoil_stop_) {
    add(t_.phonons_by_recoils, x, T_eV);
    t_.total_phonons_recoils += T_eV;
  } else if (mode_ == DamageCalculation::kinchin_pease) {
    Mover dummy;
    kp_2013(el, m_.layers[rp.layer], T_eV, x, ion_index, dummy,
            0.0 * electron_volt_per_angstrom);
  } else if (T <= el.ed && !cfg_.bugs.stale_primary_energy) {
    // Fixed stale-primary bug: an atom below Ed is not displaced; its
    // energy becomes phonons where it is.
    add(t_.phonons_by_recoils, x, T_eV);
    t_.total_phonons_recoils += T_eV;
  } else {
    start_cascade_2013(el, rp.element_layer, rp.element, site, T, ion_index, trajectory, rng,
                       false, 0.0, 0.0 * electron_volt_per_angstrom);
  }
  track_primary_stop_ = false;
  if (std::isfinite(primary_stop_x_) && cfg_.bugs.recoil_event_counted_twice) {
    add(t_.ion_stopped, primary_stop_x_, 1.0);
  }
  ion_summary(ion_index, start);
  t_.ions += 1;
}

void Transporter::run_ion_2013(std::size_t ion_index, std::vector<CollisionPoint> *trajectory) {
  // TRIM.DAT recoils (code 6) start a cascade at the atom; with monolayer
  // steps (code 7) TRIM runs each atom through the ion loop instead (its
  // internal code 8, 0x434333), as an ion that starts at the TRIM.DAT point.
  if (cfg_.source == TrimConfig::Source::particle_recoils && mode_ != DamageCalculation::monolayer) {
    run_recoil_event_2013(ion_index, trajectory);
    return;
  }
  const auto &ion = cfg_.ion;
  const bool from_list = cfg_.source != TrimConfig::Source::beam;
  Mover mv;
  mv.species = 0;
  mv.z = ion.z;
  mv.is_ion = true;
  mv.generation = 0;
  mv.track = 0;
  next_track_ = 0;
  recoils_this_ion_ = 0;
  cascades_this_ion_ = 0;
  last_T_eV_ = 0.0;
  LecuyerRandom rng = make_stream(cfg_.seed, ion_index, 0);
  details_ = ion_index < cfg_.collision_details_ions;
  const IonStart start = ion_start();
  const length_quantity l0 = m_.layers.front().monolayer;
  if (!from_list) {
    mv.x = 0.0 * angstrom;
    mv.cx = std::cos(in_units(ion.angle, radian));
    mv.cy = std::sin(in_units(ion.angle, radian));
    mv.cz = 0.0;
    mv.energy = ion.energy;
    flight_counter_ = 0;
    // TRIM starts each ion at a random point within one monolayer along its
    // lateral direction (two draws, 0x4344a2 / 0x43453c).
    mv.y = l0 * mv.cy * rng();
    mv.zz = l0 * mv.cz * rng();
  } else {
    // TRIM.DAT ions (codes 4, 5; 0x449c10): the line's energy, start and
    // direction; the ion species is the TRIM.IN ion. A zero lateral
    // coordinate is replaced by a random point within one monolayer of
    // layer 1, and an ion starting at the surface has no randomised flights
    // pending (other starts keep TRIM's counter from the previous ion).
    const SourceParticle &sp = m_.sources[ion_index / cfg_.particle_repeat].particle;
    mv.x = sp.x;
    mv.cx = sp.cos_x;
    mv.cy = sp.cos_y;
    mv.cz = sp.cos_z;
    const double norm = std::sqrt(mv.cx * mv.cx + mv.cy * mv.cy + mv.cz * mv.cz);
    if (norm > 0.0) {
      mv.cx /= norm;
      mv.cy /= norm;
      mv.cz /= norm;
    }
    mv.energy = sp.energy;
    if (sp.x == 0.0 * angstrom) flight_counter_ = 0;
    mv.y = sp.y == 0.0 * angstrom ? l0 * mv.cy * rng() : sp.y;
    mv.zz = sp.z_position == 0.0 * angstrom ? l0 * mv.cz * rng() : sp.z_position;
  }
  const bool exyz = exyz_interval_eV_ > 0.0;
  if (exyz) {
    exyz_next_ = static_cast<long>(std::floor(in_units(mv.energy, electron_volt) / exyz_interval_eV_)) + 1;
  }
  record_point(trajectory, mv, CollisionPoint::Kind::birth);
  // MD events: an ion entering below an E_c is an event at its start (the
  // handoff mode then leaves the whole ion to MD).
  const bool handed_at_entry = md_on_ && md_ion_entry(mv);
  const energy_quantity stop = ion_stop_;
  const bool full = mode_ != DamageCalculation::kinchin_pease;
  const bool restore_rng = cfg_.source == TrimConfig::Source::beam;
  const auto finish_exyz = [&]() {
    if (exyz) record_exyz(ion_index, mv, 0.0);
  };
  // The last flight, for the message if the ion's state becomes NaN.
  double last_L_A = 0.0;
  double last_de_eV = 0.0;
  double last_b = 0.0;
  double last_s2 = 0.0;

  while (!handed_at_entry) {
    // NaN never compares below the stop energy or outside the target: the
    // loop would not end. (Two comparisons per flight.)
    if (!(mv.energy >= 0.0 * electron_volt) || !std::isfinite(in_units(mv.x, angstrom))) {
      throw std::runtime_error(
          "the ion's energy or depth is not a number (E = " + num_text(in_units(mv.energy, electron_volt)) +
          " eV, x = " + num_text(in_units(mv.x, angstrom)) + " A, direction (" + num_text(mv.cx) + ", " +
          num_text(mv.cy) + ", " + num_text(mv.cz) + ") after a flight of " + num_text(last_L_A) +
          " A, electronic loss " + num_text(last_de_eV) + " eV, collision at b = " + num_text(last_b) +
          ", sin^2(theta/2) = " + num_text(last_s2) + ")");
    }
    std::size_t li = m_.layer_index(mv.x);
    // The ion's stopping power at the start of the flight (EXYZ.txt,
    // COLLISON.txt and the flight's electronic loss).
    const common::stopping_power_quantity se = stopping_at(m_.data[0][li].stopping, mv.energy);
    if (exyz) {
      // 0x43460b: a record at the start and whenever the energy has dropped
      // below the next multiple of the interval. A record at the TRIM.IN
      // energy sets the next one a step lower (Int(E/I) - 1, 0x434745), so
      // TRIM does not record the first crossing after the start.
      const double e = in_units(mv.energy, electron_volt);
      const long step = static_cast<long>(std::floor(e / exyz_interval_eV_));
      const bool at_start_energy = mv.energy == ion.energy;
      if (at_start_energy || step < exyz_next_) {
        record_exyz(ion_index, mv, in_units(se, electron_volt_per_angstrom));
        exyz_next_ = at_start_energy && cfg_.bugs.exyz_skips_first_step ? step - 1 : step;
      }
    }
    const energy_quantity e_flight = mv.energy; // at the start of the flight (time, MD events)
    const length_quantity L = ion_flight_path_2013(mv, li, rng);
    const energy_quantity de = electronic_loss_2013(mv, li, L, rng, &se);
    last_L_A = in_units(L, angstrom);
    last_de_eV = in_units(de, electron_volt);
    mv.energy -= de;
    if (mv.energy > ion.energy) mv.energy = ion.energy;
    if (mv.energy < stop) mv.energy = stop;
    const double x0 = in_units(mv.x, angstrom);
    mv.x += L * mv.cx;
    mv.y += L * mv.cy;
    mv.zz += L * mv.cz;
    advance_time(mv, L, e_flight);
    const double x1 = in_units(mv.x, angstrom);
    const double de_eV = in_units(de, electron_volt);
    const bool crosses = !((x0 < w0_ && x1 < w0_) || (x0 > w1_ && x1 > w1_));
    if (mv.x < 0.0 * angstrom) {
      if (crosses) {
        deposit(t_.ionization_by_ion, std::max(x0, 0.0), 0.0, de_eV * (x0 / std::max(x0 - x1, 1e-12)));
        t_.total_ionization_ion += de_eV * (x0 / std::max(x0 - x1, 1e-12));
      }
      leave_front(mv, ion_index);
      record_point(trajectory, mv, CollisionPoint::Kind::exit);
      finish_exyz();
      break;
    }
    if (mv.x > thickness_) {
      const double th = in_units(thickness_, angstrom);
      if (crosses) {
        const double frac = (th - x0) / std::max(x1 - x0, 1e-12);
        deposit(t_.ionization_by_ion, x0, th, de_eV * frac);
        t_.total_ionization_ion += de_eV * frac;
      }
      leave_back(mv, ion_index);
      record_point(trajectory, mv, CollisionPoint::Kind::exit);
      finish_exyz();
      break;
    }
    li = m_.layer_index(mv.x);
    const LayerModel &layer = m_.layers[li];
    const std::size_t j = pick_element(layer, rng);
    const ElementModel &el = layer.elements[j];
    const length_quantity p = boost::units::sqrt(rng() / (pi_v * L * layer.density));
    // The collision sees the energy after the flight's electronic loss.
    const CollisionPair &pair = m_.data[0][li].pairs[j];
    const double eps = pair.reduced_energy(mv.energy);
    last_b = (p / pair.screening_length).value();
    const auto kin = scatter(pair, mv.energy, eps, last_b);
    last_s2 = kin.sin2_half_theta;
    const energy_quantity T = std::min(kin.transferred, mv.energy);
    double cphi = 1.0;
    double sphi = 0.0;
    azimuth(rng, cphi, sphi);
    turn(mv.cx, mv.cy, mv.cz, kin.cos_projectile, kin.sin_projectile, -cphi, -sphi);
    mv.energy -= T;
    record_point(trajectory, mv);
    const std::uint32_t ion_point = last_point_of(trajectory, mv); // (a cascade may add points after it)
    if (trajectory != nullptr && !trajectory->empty() && trajectory->back().track == mv.track) {
      auto &pt = trajectory->back();
      pt.energy_before_eV = in_units(mv.energy + T, electron_volt);
      pt.energy_transferred_eV = in_units(T, electron_volt);
      pt.electronic_loss_eV = in_units(de, electron_volt);
    }
    const double T_eV = in_units(T, electron_volt);
    last_T_eV_ = T_eV;

    // 0x44d4e0: ionization, E2RECOIL (T >= Ed) or ion phonons (T < Ed),
    // all spread over the flight.
    if (crosses) {
      const double a = std::clamp(x0, w0_, w1_);
      const double b = std::clamp(x1, w0_, w1_);
      deposit(t_.ionization_by_ion, a, b, de_eV);
      t_.total_ionization_ion += de_eV;
      if (T >= el.ed) {
        deposit(t_.energy_to_recoils, a, b, T_eV);
        deposit(t_.energy_absorbed[el.tally_index], a, b, T_eV);
        t_.total_energy_to_recoils += T_eV;
      } else {
        deposit(t_.phonons_by_ion, a, b, T_eV);
        t_.total_phonons_ion += T_eV;
      }
    }

    if (in_window(x1)) {
      if (T < recoil_stop_) {
        // 0x4575e0 (TRIM books these transfers a second time, as recoil phonons).
        if (cfg_.bugs.ion_phonons_counted_twice) {
          add(t_.phonons_by_recoils, x1, T_eV);
          t_.total_phonons_recoils += T_eV;
        }
      } else if (!full) {
        kp_2013(el, layer, T_eV, x1, ion_index, mv, se);
      } else {
        // The primary's direction is turned from the ion's direction after
        // the collision (DIRCOS with the recoil angle).
        Mover site = mv;
        turn(site.cx, site.cy, site.cz, kin.cos_recoil, kin.sin_recoil, cphi, sphi);
        site.parent_point = last_point_of(trajectory, mv); // (site.track, time, md_mask: the ion's)
        book_pka(el, x1, layer_damage_energy(layer, T_eV)); // damage models (own tallies only)
        // Fixed stale-primary bug: a transfer of at most Ed starts no
        // cascade (it is already in the ion's phonons).
        if (T > el.ed || cfg_.bugs.stale_primary_energy) {
          start_cascade_2013(el, li, j, site, T, ion_index, trajectory, rng, restore_rng,
                             in_units(mv.energy + T, kilo_electron_volt), se);
        }
      }
    }
    // MD events: the ion's energy may have fallen below an E_c here.
    if (md_on_ && md_crossing(mv, in_units(e_flight, electron_volt))) {
      finish_exyz();
      break; // MD handoff: MD takes the ion from this collision on
    }
    if (mv.energy <= stop) {
      stop_mover(mv, true, ion_index);
      if (trajectory != nullptr && ion_point < trajectory->size()) (*trajectory)[ion_point].kind = CollisionPoint::Kind::stop;
      finish_exyz();
      break;
    }
  }
  ion_summary(ion_index, start);
  t_.ions += 1;
}

} // namespace

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

void TrimSimulation::run_block(std::size_t first_ion, std::size_t count, TrimTallies &out,
                               diag::Activity *activity, ExportBlock *exported) {
  Transporter tr(config_, *model_, out, bin_start(0), bin_width());
  const bool export_points = exported != nullptr && !config_.trajectory_export.file.empty();
  const bool export_events = exported != nullptr && config_.md_events.enabled();
  const bool export_offsets = exported != nullptr && exports_lattice_offsets(config_);
  if (export_events) tr.set_md_output(&exported->events);
  for (std::size_t i = first_ion; i < first_ion + count; ++i) {
    if (cancel_.load(std::memory_order_relaxed)) return;
    // Breadcrumb for crash reports: the ion, and the first ion of its block
    // (a rerun from there reproduces it exactly).
    if (activity != nullptr) activity->set_particle(i, first_ion, config_.seed);
    std::vector<CollisionPoint> trajectory;
    const bool keep = i - config_.first_ion < config_.trajectory_ions;
    const std::size_t first_event = export_events ? exported->events.size() : 0;
    try {
      tr.run_ion(i, keep ? &trajectory : nullptr);
    } catch (const std::exception &e) {
      throw std::runtime_error("ion " + std::to_string(i + 1) + " (seed " + std::to_string(config_.seed) +
                               "; rerun with first_ion=" + std::to_string(first_ion) + ", " +
                               std::to_string(i - first_ion + 1) + " ions): " + e.what());
    }
    if (export_events) number_events(exported->events, first_event);
    if (export_offsets) exported->lattice_offsets.emplace_back(i + 1, tr.lattice_offset_A());
    if (keep) {
      if (export_points) {
        for (const auto &pt : trajectory) {
          if (export_keeps(config_.trajectory_export, pt)) exported->points.push_back(pt);
        }
      }
      if (poll_trajectories_) {
        std::lock_guard<std::mutex> lock(trajectory_mutex_);
        trajectories_.push_back(std::move(trajectory));
      }
    }
    completed_.fetch_add(1, std::memory_order_relaxed);
  }
  if (activity != nullptr) activity->clear_particle();
}

int TrimSimulation::run() {
  const auto t0 = std::chrono::steady_clock::now();
  const std::size_t n = total_;
  const std::size_t first_ion = config_.first_ion;
  // Fixed block size: the merge order (and hence every floating-point sum)
  // does not depend on the number of threads.
  constexpr std::size_t block = 16;
  const std::size_t blocks = (n + block - 1) / block;
  unsigned threads = config_.threads != 0 ? config_.threads : std::thread::hardware_concurrency();
  threads = std::max(1u, std::min<unsigned>(threads, static_cast<unsigned>(blocks)));

  std::atomic<std::size_t> next{0};
  std::mutex merge_mutex;
  // A finished block's tallies and export rows, merged and written in block order.
  struct BlockResult {
    TrimTallies tallies;
    ExportBlock exported;
  };
  std::map<std::size_t, BlockResult> pending;
  std::size_t next_merge = 0;

  // A GPU backend (TrimConfig::compute) takes every block it can; the CPU
  // workers take the others (cpu_list). Which blocks run where depends only
  // on the configuration, so results are reproducible on a given device.
  const bool use_device = config_.compute.backend != compute::Backend::cpu &&
                          compute::unsupported_reason(config_).empty() && blocks > 0;
  std::unique_ptr<compute::Device> device;
  compute::KernelInput kernel_in;
  std::vector<std::size_t> cpu_list;
  std::vector<std::size_t> device_list;
  {
    std::lock_guard<std::mutex> lock(compute_mutex_);
    compute_report_ = {};
  }
  const auto t_setup = std::chrono::steady_clock::now();
  if (use_device) {
    kernel_in = compute::make_kernel_input(config_, *model_, n, element_names_.size(), bin_start(0), bin_width());
    compute::Options options = config_.compute;
    // fast_math's tables are the CPU's; on a GPU it selects mixed precision.
    if (config_.fast_math) options.precision = compute::Precision::mixed;
    device = compute::open_device(options);
    for (std::size_t b = 0; b < blocks; ++b) {
      (compute::block_needs_cpu(config_, b) ? cpu_list : device_list).push_back(b);
    }
    std::lock_guard<std::mutex> lock(compute_mutex_);
    compute_report_.device = device->description();
  }
  // Blocks the device gives back (too few recoil slots, a full exit list):
  // the CPU workers rerun them while the device goes on.
  std::mutex rerun_mutex;
  std::condition_variable rerun_cv;
  std::deque<std::size_t> rerun_queue;
  bool device_done = !use_device;
  // Trajectory and MD event files (opt-in), opened before any ion runs.
  std::unique_ptr<ExportWriters> writers;
  if (!config_.trajectory_export.file.empty() || config_.md_events.enabled()) {
    std::vector<std::pair<std::string, std::string>> meta = {
        {"program", std::string("AMCSET ") + AMCSET_VERSION + " (SRIM mode)"},
        {"title", config_.title},
        {"seed", std::to_string(config_.seed)},
        {"first_ion", std::to_string(first_ion + 1)},
        {"ions", std::to_string(n)},
        {"axes", "z = depth below the surface; x, y lateral (TRIM's y, z)"}};
    writers = std::make_unique<ExportWriters>(config_, meta);
  }
  std::exception_ptr error;
  std::mutex error_mutex;

  {
    std::ostringstream ctx;
    ctx << "SRIM mode, " << to_string(config_.transport) << " transport, " << to_string(config_.damage) << ", " << n << " particles from #" << first_ion
        << ", seed " << config_.seed << (config_.fast_math ? ", fast_math" : "") << ", " << threads
        << " threads\n"
        << config_.describe();
    diag::set_context(ctx.str());
  }
  AMCSET_LOG(info) << "SRIM-mode run: " << n << " particles from #" << first_ion << " in " << blocks
                   << " blocks of " << block << ", seed " << config_.seed << ", " << threads
                   << " threads" << (config_.fast_math ? ", fast approximations" : "");

  // Hand a finished block over; blocks are merged in block order.
  const auto deliver = [&](std::size_t b, BlockResult &&result) {
    std::lock_guard<std::mutex> lock(merge_mutex);
    pending.emplace(b, std::move(result));
    while (!pending.empty() && pending.begin()->first == next_merge) {
      {
        std::lock_guard<std::mutex> tlock(tallies_mutex_);
        tallies_.add(pending.begin()->second.tallies, config_.max_exit_records);
        transmitted_.store(tallies_.transmitted, std::memory_order_relaxed);
      }
      if (writers) writers->write(pending.begin()->second.exported);
      pending.erase(pending.begin());
      ++next_merge;
    }
  };
  const auto fail = [&](const char *who) {
    {
      std::lock_guard<std::mutex> lock(error_mutex);
      if (!error) error = std::current_exception();
    }
    cancel_.store(true);
    pause_cv_.notify_all();
    {
      std::lock_guard<std::mutex> lock(rerun_mutex);
      rerun_cv.notify_all();
    }
    AMCSET_LOG(error) << who << " failed";
  };
  const auto wait_if_paused = [&]() {
    if (paused_.load()) {
      std::unique_lock<std::mutex> lock(pause_mutex_);
      pause_cv_.wait(lock, [&] { return !paused_.load() || cancel_.load(); });
    }
  };
  auto worker = [&]() {
    diag::Activity activity("SRIM worker");
    std::size_t done_blocks = 0;
    try {
      while (true) {
        wait_if_paused();
        if (cancel_.load()) break;
        const std::size_t i = next.fetch_add(1);
        std::size_t b = i;
        if (!use_device) {
          if (i >= blocks) break;
        } else if (i < cpu_list.size()) {
          b = cpu_list[i];
        } else {
          // The blocks the device gives back, until it is done.
          std::unique_lock<std::mutex> lock(rerun_mutex);
          rerun_cv.wait(lock, [&] { return !rerun_queue.empty() || device_done || cancel_.load(); });
          if (rerun_queue.empty() || cancel_.load()) break;
          b = rerun_queue.front();
          rerun_queue.pop_front();
        }
        BlockResult result;
        result.tallies.resize(config_.depth_bins, element_names_.size());
        const std::size_t first = b * block;
        run_block(first_ion + first, std::min(block, n - first), result.tallies, &activity,
                  writers ? &result.exported : nullptr);
        ++done_blocks;
        deliver(b, std::move(result));
      }
    } catch (...) {
      fail("SRIM worker");
    }
    if (use_device) {
      std::lock_guard<std::mutex> lock(compute_mutex_);
      compute_report_.cpu_blocks += done_blocks;
    }
    AMCSET_LOG(debug) << "SRIM worker finished " << done_blocks << " blocks";
  };

  // The device: lanes (work items) take the device's blocks one after the
  // other and step through them in launches of `steps` steps, sized from the
  // launch times (0.05 to 0.2 s: no launch nears the display driver's
  // watchdog, and other programs on the GPU wait little). A finished block
  // waits in its ring slot until every block before it has finished; then a
  // merge kernel adds it to the running
  // total on the device, in block order (the CPU's order of additions). Per
  // launch the host reads a few numbers, the status of the merged blocks and
  // their packed lists (only if the scenario keeps them); the total itself
  // comes to the host at most twice a second (live plots) and at the end.
  auto device_worker = [&]() {
    diag::Activity activity("SRIM device");
    try {
      const std::size_t budget =
          (config_.compute.memory_mib != 0 ? config_.compute.memory_mib : std::size_t{2048}) << 20;
      const std::size_t n_dev = device_list.size();
      const std::size_t first_dev = device_list.front();
      for (std::size_t i = 0; i < n_dev; ++i) {
        if (device_list[i] != first_dev + i) throw std::logic_error("SRIM device: the device's blocks are not consecutive");
      }
      const std::size_t lb = compute::lane_bytes(kernel_in);
      const std::size_t bb = compute::block_bytes(kernel_in);
      // Lanes take up to 40 % of the memory, the ring of block slots the rest
      // (at least as many slots as lanes, and half as many again if it fits).
      std::size_t lanes = std::max<std::size_t>(1, std::min({device->preferred_lanes(), budget * 2 / 5 / lb, n_dev}));
      std::size_t ring = std::max<std::size_t>(1, std::min(n_dev, (budget - std::min(budget, lanes * lb)) / bb));
      ring = std::min(ring, std::max(lanes + lanes / 2, lanes));
      lanes = std::min(lanes, ring);
      device->load(kernel_in, ring, lanes);
      std::size_t steps = 2000;
      // Profiling (docs/GPU.md "Profile"): a fixed number of steps per launch,
      // so that a profiler can capture short dispatches.
      std::size_t fixed_steps = 0;
      if (const char *e = std::getenv("AMCSET_GPU_STEPS"); e != nullptr && e[0] != 0) {
        fixed_steps = static_cast<std::size_t>(std::strtoull(e, nullptr, 10));
        if (fixed_steps > 0) steps = fixed_steps;
      }
      // The CPU's blocks come first (the blocks whose paths are drawn): the
      // device's total starts where they end.
      while (!cancel_.load()) {
        {
          std::lock_guard<std::mutex> lock(merge_mutex);
          if (next_merge >= first_dev) break;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
      }
      {
        std::vector<double> flat;
        {
          std::lock_guard<std::mutex> lock(tallies_mutex_);
          flat = compute::pack_tallies(kernel_in, tallies_);
        }
        device->set_total(flat.data());
      }
      const auto t_dev = std::chrono::steady_clock::now();
      struct {
        double kernel = 0.0, read = 0.0, host = 0.0;
        double steps = 0.0, max_steps = 0.0;
        std::size_t launches = 0, deepest = 0, merges = 0;
      } prof;
      {
        std::lock_guard<std::mutex> lock(compute_mutex_);
        compute_report_.setup_seconds = std::chrono::duration<double>(t_dev - t_setup).count();
        compute_report_.lanes = lanes;
        compute_report_.batch_blocks = ring;
      }
      const std::vector<std::size_t> rows = {
          compute::scalar_row(kernel_in, AMK_S_STATUS), compute::scalar_row(kernel_in, AMK_S_POSITIONS),
          compute::scalar_row(kernel_in, AMK_S_EXITS), compute::scalar_row(kernel_in, AMK_S_STEPS),
          compute::scalar_row(kernel_in, AMK_S_MAX_K), compute::scalar_row(kernel_in, AMK_S_EXYZ)};
      const bool exyz = kernel_in.params.exyz_capacity > 0;
      const bool keep_lists = config_.record_ion_positions || config_.record_exiting_atoms || exyz;
      std::vector<double> R(rows.size() * ring);
      std::vector<double> total(static_cast<std::size_t>(kernel_in.params.tally_stride));
      std::vector<double> packed_positions;
      std::vector<double> packed_exits;
      const auto publish = [&]() {
        const auto t_read = std::chrono::steady_clock::now();
        device->read_total(total.data());
        std::lock_guard<std::mutex> lock(tallies_mutex_);
        compute::unpack_block(kernel_in, total.data(), nullptr, tallies_); // the lists stay
        transmitted_.store(tallies_.transmitted, std::memory_order_relaxed);
        prof.read += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_read).count();
      };
      auto last_publish = std::chrono::steady_clock::now();
      std::size_t merged = 0;
      bool init = true;
      while (merged < n_dev && !cancel_.load()) {
        wait_if_paused();
        if (cancel_.load()) break;
        const auto t0l = std::chrono::steady_clock::now();
        const std::size_t done = device->launch(first_dev, n_dev, lanes, merged + ring, steps, init);
        init = false;
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0l).count();
        prof.kernel += dt;
        ++prof.launches;
        // Launches of 0.05 to 0.2 s: far from Windows' display watchdog (2 s),
        // and short enough that other programs on the GPU are not held up
        // for long. Results do not depend on the steps per launch.
        if (fixed_steps == 0) {
          if (dt > 0.2 && steps > 64) steps /= 2;
          else if (dt < 0.05 && steps < (std::size_t{1} << 24)) steps *= 2;
        }
        const auto t1 = std::chrono::steady_clock::now();
        const std::size_t count = device->scan(merged, n_dev);
        if (count == 0) {
          if (done >= lanes) throw std::logic_error("SRIM device: every lane finished, but blocks were not merged");
          continue;
        }
        device->read_rows(rows, R.data());
        const auto row = [&](std::size_t r, std::size_t slot) { return R[r * ring + slot]; };
        // Blocks the device gave back are rerun here and put back in their
        // slot, so that the merge adds them in order.
        // Offsets: positions and exits per block; with EXYZ the start of its
        // records in the exit buffer and each block's first record.
        std::vector<std::int64_t> offsets(2 * count + (exyz ? count + 1 : 0));
        std::size_t n_exyz = 0;
        std::map<std::size_t, TrimTallies> reran;
        std::size_t n_pos = 0;
        std::size_t n_exit = 0;
        std::size_t ions_done = 0;
        for (std::size_t k = 0; k < count; ++k) {
          const std::size_t slot = (merged + k) % ring;
          const std::size_t b = first_dev + merged + k;
          offsets[2 * k] = static_cast<std::int64_t>(n_pos);
          offsets[2 * k + 1] = static_cast<std::int64_t>(n_exit);
          prof.steps += row(3, slot);
          prof.max_steps = std::max(prof.max_steps, row(3, slot));
          prof.deepest = std::max(prof.deepest, static_cast<std::size_t>(row(4, slot)));
          if (static_cast<int>(row(0, slot)) != AMK_OK) {
            TrimTallies &tt = reran[k];
            tt.resize(config_.depth_bins, element_names_.size());
            run_block(first_ion + b * block, std::min(block, n - b * block), tt, &activity, nullptr);
            device->put_block(slot, compute::pack_tallies(kernel_in, tt).data());
            continue;
          }
          n_pos += static_cast<std::size_t>(row(1, slot));
          n_exit += static_cast<std::size_t>(row(2, slot));
          if (exyz) {
            offsets[2 * count + 1 + k] = static_cast<std::int64_t>(n_exyz);
            n_exyz += static_cast<std::size_t>(row(5, slot));
          }
          ions_done += std::min(block, n - b * block);
        }
        device->merge(merged, count);
        if (keep_lists) {
          const std::size_t exit_doubles = n_exit * AMK_EXIT_DOUBLES + n_exyz * AMK_EXYZ_DOUBLES;
          if (exyz) offsets[2 * count] = static_cast<std::int64_t>(n_exit * AMK_EXIT_DOUBLES);
          packed_positions.resize(std::max<std::size_t>(1, n_pos) * AMK_POSITION_DOUBLES);
          packed_exits.resize(std::max<std::size_t>(1, exit_doubles));
          if (n_pos + n_exit + n_exyz > 0) {
            device->pack_lists(merged, count, offsets, n_pos, exit_doubles, packed_positions.data(),
                               packed_exits.data());
          }
          std::lock_guard<std::mutex> lock(tallies_mutex_);
          for (std::size_t k = 0; k < count; ++k) {
            if (auto it = reran.find(k); it != reran.end()) {
              tallies_.positions.insert(tallies_.positions.end(), it->second.positions.begin(), it->second.positions.end());
              for (const auto &r : it->second.exits) {
                if (tallies_.exits.size() >= config_.max_exit_records) break;
                tallies_.exits.push_back(r);
              }
              for (const auto &r : it->second.exyz) {
                if (tallies_.exyz.size() >= config_.max_exit_records) break;
                tallies_.exyz.push_back(r);
              }
              continue;
            }
            const std::size_t slot = (merged + k) % ring;
            compute::append_positions(packed_positions.data() + static_cast<std::size_t>(offsets[2 * k]) * AMK_POSITION_DOUBLES,
                                      static_cast<std::size_t>(row(1, slot)), tallies_.positions);
            compute::append_exits(packed_exits.data() + static_cast<std::size_t>(offsets[2 * k + 1]) * AMK_EXIT_DOUBLES,
                                  static_cast<std::size_t>(row(2, slot)), config_.max_exit_records, tallies_.exits);
            if (exyz) {
              compute::append_exyz(packed_exits.data() + n_exit * AMK_EXIT_DOUBLES +
                                       static_cast<std::size_t>(offsets[2 * count + 1 + k]) * AMK_EXYZ_DOUBLES,
                                   static_cast<std::size_t>(row(5, slot)), config_.max_exit_records, tallies_.exyz);
            }
          }
        }
        merged += count;
        ++prof.merges;
        completed_.fetch_add(ions_done, std::memory_order_relaxed);
        prof.host += std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
        if (std::chrono::steady_clock::now() - last_publish > std::chrono::milliseconds(500)) {
          publish();
          last_publish = std::chrono::steady_clock::now();
        }
        std::lock_guard<std::mutex> lock(compute_mutex_);
        compute_report_.device_blocks += count - reran.size();
        compute_report_.cpu_blocks += reran.size();
        compute_report_.reruns += reran.size();
        compute_report_.deepest_slot = prof.deepest;
        compute_report_.kernel_seconds = prof.kernel;
        compute_report_.read_seconds = prof.read;
        compute_report_.host_seconds = prof.host;
        compute_report_.launches = prof.launches;
        compute_report_.total_steps = prof.steps;
        compute_report_.block_steps_max = prof.max_steps;
        compute_report_.block_steps_mean = prof.steps / static_cast<double>(merged);
        compute_report_.device_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_dev).count();
      }
      if (!cancel_.load() || merged > 0) publish(); // the run's (or the stopped run's merged) results
      std::lock_guard<std::mutex> lock(compute_mutex_);
      compute_report_.read_seconds = prof.read;
      compute_report_.merges = prof.merges;
      compute_report_.device_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_dev).count();
    } catch (...) {
      fail("SRIM device");
    }
    std::lock_guard<std::mutex> lock(rerun_mutex);
    device_done = true;
    rerun_cv.notify_all();
  };
  std::thread device_thread;
  if (use_device) device_thread = std::thread(device_worker);

  // A thread that cannot be created leaves its blocks to the others (the
  // results do not depend on the number of threads).
  std::vector<std::thread> pool;
  pool.reserve(threads);
  for (unsigned i = 1; i < threads; ++i) {
    try {
      pool.emplace_back(worker);
    } catch (const std::system_error &e) {
      AMCSET_LOG(warning) << "could only start " << pool.size() + 1 << " of " << threads
                          << " worker threads: " << e.what();
      break;
    }
  }
  worker();
  for (auto &t : pool) t.join();
  if (device_thread.joinable()) device_thread.join();
  // Blocks finished out of order after a cancellation are merged as well.
  {
    std::lock_guard<std::mutex> tlock(tallies_mutex_);
    for (auto &kv : pending) tallies_.add(kv.second.tallies, config_.max_exit_records);
    transmitted_.store(tallies_.transmitted, std::memory_order_relaxed);
  }
  if (writers && !error) {
    try {
      for (auto &kv : pending) writers->write(kv.second.exported);
      writers->finish();
    } catch (...) {
      error = std::current_exception();
    }
  }
  pending.clear();
  elapsed_seconds_ =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  if (error) {
    AMCSET_LOG(error) << "SRIM-mode run failed after " << completed_.load() << " particles";
    std::rethrow_exception(error);
  }
  AMCSET_LOG(info) << "SRIM-mode run " << (cancel_.load() ? "stopped" : "finished") << ": "
                   << completed_.load() << " particles in " << elapsed_seconds_.load() << " s";
  if (!config_.output_directory.empty()) {
    const auto w0 = std::chrono::steady_clock::now();
    write_outputs(config_.output_directory);
    AMCSET_LOG(info) << "output written to " << config_.output_directory << " in "
                     << std::chrono::duration<double>(std::chrono::steady_clock::now() - w0).count()
                     << " s";
  }
  return 0;
}

Tallies TrimSimulation::snapshot_tallies() const {
  // The SimulationEngine statistics (SI units: depth in metres, energy in
  // joules) from the merged TRIM tallies, on TRIM's depth bins. They are
  // derived here rather than filled per collision, which kept the transport
  // loops busy with a second set of histograms.
  const TrimTallies t = snapshot(false);
  const double w0 = in_units(config_.window_start, meter);
  const double w1 = in_units(config_.window_end > 0.0 * angstrom ? config_.window_end : model_->thickness, meter);
  const std::size_t bins = config_.depth_bins;
  Tallies live(bins, w1, bins, bins, bins, bins);
  for (auto *h : {&live.range_histogram, &live.depth_histogram, &live.displacement_histogram,
                  &live.electronic_energy_loss_histogram, &live.nuclear_energy_loss_histogram,
                  &live.total_energy_loss_histogram}) {
    h->min = w0;
  }
  live.incidents = t.ions;
  live.total_displacements = static_cast<std::size_t>(std::llround(std::max(0.0, t.total_displacements)));
  double sputtered = 0.0;
  for (double v : t.sputtered) sputtered += v;
  live.sputtered = static_cast<std::size_t>(std::llround(sputtered));
  if (t.range.n > 0.0) {
    const double m_per_A = in_units(1.0 * angstrom, meter);
    live.range_stats.n = static_cast<std::size_t>(t.range.n);
    live.range_stats.mean = t.range.mean() * m_per_A;
    const double var = t.range.stddev() * t.range.stddev();
    live.range_stats.M2 = var * t.range.n * m_per_A * m_per_A;
  }
  const double j_per_eV = in_units(1.0 * electron_volt, joule);
  for (std::size_t i = 0; i < bins; ++i) {
    const auto v = [&](const DepthHistogram &h) { return i < h.values.size() ? h.values[i] : 0.0; };
    live.range_histogram.bins[i] = v(t.ion_stopped);
    live.displacement_histogram.bins[i] = v(t.displacements);
    live.depth_histogram.bins[i] = v(t.ion_stopped) + v(t.displacements);
    const double el = (v(t.ionization_by_ion) + v(t.ionization_by_recoils)) * j_per_eV;
    const double nu = (v(t.phonons_by_ion) + v(t.phonons_by_recoils)) * j_per_eV;
    live.electronic_energy_loss_histogram.bins[i] = el;
    live.nuclear_energy_loss_histogram.bins[i] = nu;
    live.total_energy_loss_histogram.bins[i] = el + nu;
  }
  return live;
}

bool TrimSimulation::try_pop_trajectory(std::vector<CollisionPoint> &out) {
  std::lock_guard<std::mutex> lock(trajectory_mutex_);
  if (trajectories_.empty()) return false;
  out = std::move(trajectories_.front());
  trajectories_.pop_front();
  return true;
}

void TrimSimulation::clear_trajectories() {
  std::lock_guard<std::mutex> lock(trajectory_mutex_);
  trajectories_.clear();
}

} // namespace trim
} // namespace amcset
