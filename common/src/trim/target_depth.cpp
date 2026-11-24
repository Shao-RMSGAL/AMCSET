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

/*! \file target_depth.cpp
 *  \brief The ions' range estimate, plausibility warnings of a configuration
 *  and the automatic target depth (see trim.h).
 */

#include "trim/trim.h"
#include "crystal/materials.h"

#include "data_files.h"
#include "diagnostics.h"
#include "element_data.h"
#include "physics/scattering.h"
#include "physics/stopping.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>

namespace amcset {
namespace trim {

using namespace common;

namespace {

std::string num_text(double v) {
  char buf[48];
  std::snprintf(buf, sizeof buf, "%.3g", v);
  return buf;
}

//! A length for messages, in a fitting unit.
std::string length_text(double angstrom_value) {
  const double a = std::abs(angstrom_value);
  if (a >= 1e10) return num_text(angstrom_value / 1e10) + " m";
  if (a >= 1e7) return num_text(angstrom_value / 1e7) + " mm";
  if (a >= 1e4) return num_text(angstrom_value / 1e4) + " um";
  if (a >= 10.0) return num_text(angstrom_value / 10.0) + " nm";
  return num_text(angstrom_value) + " A";
}

//! The stopping of the ion in one layer: electronic (Bragg's rule with the
//! configured model) plus ZBL nuclear stopping, in eV per Angstrom.
class LayerStopping {
public:
  LayerStopping(const TrimConfig &c, const TargetLayer &l, physics::ElectronicStoppingModel model,
                const physics::SrimStoppingDatabase *db)
      : c_(c), l_(l), model_(model), db_(db) {
    double total = 0.0;
    for (const auto &e : l.elements) total += e.stoichiometry;
    mass_quantity mean = 0.0 * kilogram;
    for (const auto &e : l.elements) {
      const double f = e.stoichiometry / total;
      comps_.push_back({e.z, f});
      pairs_.emplace_back(c.ion.z, c.ion.mass, e.z, e.mass);
      fractions_.push_back(f);
      mean += f * e.mass;
    }
    density_ = l.density / mean;
  }
  //! Total stopping; \p nuclear (optional) receives the nuclear part.
  double operator()(double e_eV, double *nuclear = nullptr) const {
    const energy_quantity e = e_eV * electron_volt;
    const auto se = physics::layer_stopping(model_, db_, c_.ion.z, c_.ion.mass, comps_, l_.gas,
                                            l_.compound_correction, e, density_, c_.bethe_bloch);
    double sn = 0.0;
    for (std::size_t j = 0; j < pairs_.size(); ++j) {
      sn += fractions_[j] * in_units(density_ * pairs_[j].nuclear_stopping(e), electron_volt_per_angstrom);
    }
    if (nuclear != nullptr) *nuclear = sn;
    return in_units(density_ * se, electron_volt_per_angstrom) + sn;
  }

private:
  const TrimConfig &c_;
  const TargetLayer &l_;
  physics::ElectronicStoppingModel model_;
  const physics::SrimStoppingDatabase *db_;
  std::vector<physics::StoppingComponent> comps_;
  std::vector<physics::CollisionPair> pairs_;
  std::vector<double> fractions_;
  number_density_quantity density_;
};

double monolayer_A(const TargetLayer &l) {
  double total = 0.0;
  for (const auto &e : l.elements) total += e.stoichiometry;
  mass_quantity mean = 0.0 * kilogram;
  for (const auto &e : l.elements) mean += (e.stoichiometry / total) * e.mass;
  return in_units(1.0 / boost::units::root<3>(l.density / mean), angstrom);
}

struct Penetration {
  double path_A = 0.0;          //!< Path length (CSDA range)
  double nuclear_share = 0.0;   //!< Path-averaged S_n / (S_e + S_n)
};

Penetration penetration(const TrimConfig &c);

} // namespace

length_quantity estimate_penetration(const TrimConfig &c) { return penetration(c).path_A * angstrom; }

length_quantity estimate_deepest(const TrimConfig &c) {
  const auto p = penetration(c);
  // Range straggling: with mostly electronic stopping (light, fast ions) the
  // deepest of ~1000 ions lies within ~15 % beyond the mean path; nuclear
  // collisions are few and large, so where they dominate the stopping the
  // deepest ions go up to ~1.7 times as far (Fe -> Fe 100 keV: path 46 nm,
  // deepest of 300 ions at 72 nm with the built-in stopping).
  return p.path_A * (1.2 + p.nuclear_share) * angstrom;
}

namespace {

Penetration penetration(const TrimConfig &c) {
  auto model = c.stopping_model;
  std::shared_ptr<const physics::SrimStoppingDatabase> db;
  if (model == physics::ElectronicStoppingModel::srim_tables) {
    const std::string path = c.stopping_table.empty() ? data::bundled_srim_stopping_table() : c.stopping_table;
    if (path.empty()) {
      model = physics::ElectronicStoppingModel::builtin;
    } else {
      db = physics::SrimStoppingDatabase::load(path);
    }
  }
  const energy_quantity stop_q = c.ion_stop_energy > 0.0 * electron_volt ? c.ion_stop_energy : c.final_energy;
  const double stop = std::max(in_units(stop_q, electron_volt), 1e-3);
  double e = in_units(c.ion.energy, electron_volt);
  double path = 0.0;    // Angstrom
  double nuclear = 0.0; // integral of S_n/S dx (Angstrom)
  const auto infinite = Penetration{std::numeric_limits<double>::infinity(), 0.0};
  // Logarithmic energy steps (100 per decade), the stopping at the middle of
  // each step: dx = dE / S.
  const double factor = std::exp(-std::log(10.0) / 100.0);
  for (std::size_t i = 0; i < c.layers.size() && e > stop; ++i) {
    const bool last = i + 1 == c.layers.size();
    const LayerStopping s(c, c.layers[i], model, db.get());
    double left = in_units(c.layers[i].width, angstrom);
    for (int guard = 0; e > stop && guard < 100000; ++guard) {
      const double e_next = std::max(e * factor, stop);
      const double de = e - e_next;
      double sn = 0.0;
      const double sp = s(0.5 * (e + e_next), &sn);
      if (!(sp > 0.0) || !std::isfinite(sp)) return infinite;
      const double dx = de / sp;
      if (!last && dx >= left) {
        // Leaves the layer: the energy lost over the rest of it.
        e -= sp * left;
        path += left;
        nuclear += left * sn / sp;
        break;
      }
      path += dx;
      nuclear += dx * sn / sp;
      left -= dx;
      e = e_next;
    }
  }
  return Penetration{path, path > 0.0 ? nuclear / path : 0.0};
}

} // namespace

std::vector<std::string> TrimConfig::plausibility_warnings() const {
  std::vector<std::string> w;
  try {
    if (transport == Transport::crystal) {
      bool any = false;
      for (const auto &l : layers) any = any || l.crystal.has_value();
      if (!any) {
        w.push_back("physics.transport is \"crystal\", but no layer has a crystal: every layer is treated as "
                    "amorphous (random collision partners in monolayer steps).");
      } else {
        // The numbers: tools/crystal/calibration.py (docs/CLI.md "Crystals").
        w.push_back("The crystal transport is AMCSET's own (not TRIM). It reproduces a measured channeling dip "
                    "and measured boron profiles in silicon within the spread of those data, but how the "
                    "electronic stopping is divided between the path and the collisions "
                    "(physics.crystal.local_stopping_fraction, local_stopping_decay) is not calibrated: the "
                    "measurements allow 0.25 to 1 (default 0.5), and over that range the depth of a channeling "
                    "tail changes by about 10 % for 30 keV B in Si and by up to 50 % for 0.5 keV B. Peak "
                    "positions and trends (direction, temperature) depend on it much less.");
      }
    }
    double thickness = 0.0;
    for (std::size_t i = 0; i < layers.size(); ++i) {
      const auto &l = layers[i];
      const std::string where = "layer " + std::to_string(i + 1) + " ('" + l.name + "')";
      const double rho = in_units(l.density, gram_per_cubic_centimeter);
      thickness += in_units(l.width, angstrom);
      if (l.crystal && transport != Transport::crystal) {
        // Say so rather than let a crystal scenario look like a channeling
        // calculation.
        w.push_back(where + ": the crystal (" + (l.crystal->cell.name.empty() ? std::string("unit cell") : l.crystal->cell.name) +
                    ") is not used by this calculation. The " + to_string(transport) + " transport treats every "
                    "layer as amorphous (as TRIM does), with the unit cell's composition and density; set "
                    "physics.transport to \"crystal\" to follow the lattice.");
      }
      if (l.crystal && transport == Transport::crystal && crystal.thermal_vibrations) {
        // Atoms without thermal data: the tabulated Debye temperature of a
        // pure element (said, since it was not asked for), else they stand
        // still and channeling is too strong.
        std::string still;
        for (const auto &e : l.elements) {
          bool from_table = false;
          const auto u = crystal::thermal_amplitude(*l.crystal, e.z, e.mass, &from_table);
          const std::string symbol(Element::by_z(e.z).symbol());
          if (!u) {
            still += (still.empty() ? "" : ", ") + symbol;
          } else if (from_table) {
            w.push_back(where + ": no thermal vibrations are given for " + symbol + " (crystal.vibrations): using "
                        "the Debye temperature of Gemmell's table (Rev. Mod. Phys. 46, 129 (1974)), " +
                        num_text(in_units(*crystal::reference_debye_temperature(e.z), kelvin)) +
                        " K, which gives an rms displacement of " + num_text(in_units(*u, angstrom)) +
                        " A along each axis at " + num_text(in_units(l.crystal->temperature, kelvin)) + " K.");
          }
        }
        if (!still.empty()) {
          w.push_back(where + ": no thermal vibrations are given for " + still + " (crystal.vibrations): these "
                      "atoms stand exactly on their sites, which makes channeling stronger than in a real crystal.");
        }
      }
      if (l.elements.empty() || !(rho > 0.0)) continue;
      // Reference density: the element's (single elements), or the mixture
      // of the elements' condensed densities by volume (compounds can differ
      // from it by a factor of two or so, hence the wider band).
      double ref = 0.0;
      double band = 2.0;
      if (l.elements.size() == 1) {
        const auto &el = Element::by_z(l.elements.front().z);
        ref = in_units(l.gas ? el.gas_density() : el.solid_density(), gram_per_cubic_centimeter);
      } else if (!l.gas) {
        double mass = 0.0;
        double volume = 0.0;
        for (const auto &e : l.elements) {
          const double d = in_units(Element::by_z(e.z).solid_density(), gram_per_cubic_centimeter);
          const double m = e.stoichiometry * (e.mass / atomic_mass_unit).value();
          if (d > 0.0) {
            mass += m;
            volume += m / d;
          }
        }
        if (volume > 0.0) ref = mass / volume;
        band = 3.0;
      }
      if (ref > 0.0 && (rho > band * ref || rho < ref / band)) {
        w.push_back(where + ": the density " + num_text(rho) + " g/cm3 is far from " +
                    (l.elements.size() == 1 ? "that of " + std::string(Element::by_z(l.elements.front().z).symbol())
                                            : std::string("what its elements suggest")) +
                    (l.gas ? " as a gas" : "") + " (" + num_text(ref) + " g/cm3); is it meant?");
      }
      if (l.gas && rho > 0.05) {
        w.push_back(where + ": marked as a gas but " + num_text(rho) +
                    " g/cm3 is the density of a condensed material (gases at normal pressure: ~0.001 "
                    "g/cm3); untick Gas if it is a solid");
      } else if (!l.gas && rho < 0.05) {
        w.push_back(where + ": " + num_text(rho) +
                    " g/cm3 is the density of a gas, but the layer is marked as a solid (the stopping "
                    "differs); tick Gas if it is one");
      }
    }
    // Ion energy and target depth.
    const double reach = in_units(estimate_penetration(*this), angstrom);
    const double deepest = in_units(estimate_deepest(*this), angstrom);
    if (std::isfinite(reach) && reach > 0.0 && source == Source::beam) {
      const double window = window_end > 0.0 * angstrom ? in_units(window_end - window_start, angstrom)
                                                        : thickness - in_units(window_start, angstrom);
      if (thickness > 20.0 * reach && thickness - reach > 1000.0) {
        w.push_back("the target is " + length_text(thickness) + " thick, but ions of " +
                    num_text(in_units(ion.energy, kilo_electron_volt)) + " keV travel at most about " +
                    length_text(reach) + ": the depth tables (" + std::to_string(depth_bins) +
                    " bins over the depth window) will be far too coarse. Use about " +
                    length_text(deepest) + ", or the automatic target depth");
      } else if (window > 20.0 * reach && window - reach > 1000.0) {
        w.push_back("the depth window is " + length_text(window) + " wide, but the ions reach only about " +
                    length_text(reach) + ": the depth tables will be coarse");
      }
      if (thickness < 0.3 * reach && !auto_depth.enabled) {
        w.push_back("note: ions of this energy travel about " + length_text(reach) + ", much more than the " +
                    length_text(thickness) + " target: most will be transmitted (fine for transmission "
                    "studies; otherwise use a thicker target or the automatic target depth)");
      }
    }
    // Relativistic projectiles: the collision kinematics are classical.
    const double mc2 = in_units(ion.mass * speed_of_light * speed_of_light, electron_volt);
    if (mc2 > 0.0 && in_units(ion.energy, electron_volt) > 0.05 * mc2) {
      w.push_back("the ion is relativistic (kinetic energy " +
                  num_text(100.0 * in_units(ion.energy, electron_volt) / mc2) +
                  " % of its rest energy): electronic stopping is relativistic (Bethe-Bloch), but the "
                  "nuclear collisions use classical kinematics, so damage and scattering are approximate");
    }
  } catch (const std::exception &e) {
    w.push_back(std::string("the plausibility checks could not be completed: ") + e.what());
  } catch (...) {
    w.push_back("the plausibility checks could not be completed");
  }
  return w;
}

AutoDepthResult size_target(TrimConfig &c, const std::atomic<bool> *cancel) {
  AutoDepthResult r;
  if (c.layers.empty()) throw std::invalid_argument("auto_depth: the target has no layers");
  r.last_layer_width = c.layers.back().width;
  if (c.source == TrimConfig::Source::particle_recoils) {
    r.notes.push_back("automatic target depth: not used for a recoil source (its atoms start at given depths)");
    return r;
  }
  const auto &opt = c.auto_depth;
  if (!(opt.margin > 0.0 && opt.margin < 1.0)) {
    throw std::invalid_argument("auto_depth.margin: must be above 0 and below 1, got " + num_text(opt.margin));
  }
  if (opt.test_ions == 0) throw std::invalid_argument("auto_depth.test_ions: must be positive");
  double top = 0.0; // Angstrom above the last layer
  for (std::size_t i = 0; i + 1 < c.layers.size(); ++i) top += in_units(c.layers[i].width, angstrom);
  const double min_last = std::ceil(monolayer_A(c.layers.back()) * 1.01 * 100.0) / 100.0;

  // Start: the estimated depth of the deepest ions (estimate_deepest; the
  // path length for the check that they reach the last layer at all).
  const double reach = in_units(estimate_penetration(c), angstrom);
  const double deepest_guess = in_units(estimate_deepest(c), angstrom);
  double last = in_units(c.layers.back().width, angstrom);
  if (std::isfinite(reach) && reach > 0.0) {
    if (reach < top) {
      r.notes.push_back("automatic target depth: the ions are expected to stop above the last layer (reach about " +
                        length_text(reach) + "); its thickness is left as it is");
      return r;
    }
    last = std::max(min_last, deepest_guess - top);
  }
  // Quick-mode test runs (fast approximations: only the depths matter).
  TrimConfig t = c;
  t.damage = DamageCalculation::kinchin_pease;
  t.ion_count = t.source == TrimConfig::Source::beam ? opt.test_ions
                                                     : std::min<std::size_t>(opt.test_ions, t.particles.size());
  t.first_ion = 0;
  t.output_directory.clear();
  t.record_ion_positions = true;
  t.record_exiting_atoms = false;
  t.collision_details_ions = 0;
  t.trajectory_ions = 0;
  t.exyz_interval = 0.0 * electron_volt;
  t.window_start = 0.0 * angstrom;
  t.window_end = 0.0 * angstrom;
  t.fast_math = true;
  t.auto_depth.enabled = false;
  t.trajectory_export = {}; // the test runs write no files
  t.md_events = {};
  constexpr double max_total_A = 1e11; // 10 m
  bool done = false;
  for (int run = 0; run < 12 && !done; ++run) {
    if (cancel != nullptr && cancel->load()) {
      r.notes.push_back("automatic target depth: cancelled");
      return r;
    }
    t.layers.back().width = last * angstrom;
    TrimSimulation sim(t);
    sim.run();
    ++r.runs;
    double deepest = 0.0;
    for (const auto &p : sim.tallies().positions) deepest = std::max(deepest, p.x);
    r.deepest = deepest * angstrom;
    r.transmitted = sim.tallies().transmitted;
    const double total = top + last;
    AMCSET_LOG(info) << "auto depth run " << r.runs << ": last layer " << last << " A, deepest " << deepest
                     << " A, transmitted " << r.transmitted;
    if (r.transmitted > 0) {
      if (total >= max_total_A) {
        r.notes.push_back("automatic target depth: ions are still transmitted through " + length_text(total) +
                          "; the depth was not increased further");
        break;
      }
      last = std::min(max_total_A - top, std::max(2.0 * last, 2.0 * total - top));
      continue;
    }
    if (deepest < top) {
      r.notes.push_back("automatic target depth: the ions stop above the last layer (deepest " +
                        length_text(deepest) + "); its thickness is left at " + length_text(last));
      done = true;
      break;
    }
    if (deepest >= (1.0 - opt.margin) * total * (1.0 - 1e-12)) {
      done = true;
      break;
    }
    // Shrink: the deepest ion at the start of the last margin (the most
    // room behind it for the deeper tail of a longer run).
    last = std::max(min_last, deepest / (1.0 - opt.margin) - top);
  }
  if (!done && r.notes.empty()) {
    r.notes.push_back("automatic target depth: no stable depth after " + std::to_string(r.runs) +
                      " test runs; the last one is used");
  }
  const double before = in_units(c.layers.back().width, angstrom);
  c.layers.back().width = last * angstrom;
  r.last_layer_width = c.layers.back().width;
  r.changed = std::abs(last - before) > 1e-9 * std::max(last, before);
  if (c.window_end > c.target_thickness()) {
    c.window_end = 0.0 * angstrom;
    r.notes.push_back("automatic target depth: the depth window now ends at the back of the target");
  }
  r.notes.insert(r.notes.begin(), "automatic target depth: last layer ('" + c.layers.back().name + "') " +
                                      length_text(before) + " -> " + length_text(last) + " (deepest of " +
                                      std::to_string(t.ion_count) + " test ions at " +
                                      length_text(in_units(r.deepest, angstrom)) + ", " +
                                      std::to_string(r.runs) + " quick test runs)");
  return r;
}

} // namespace trim
} // namespace amcset
