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

/*! \file layer_crystal.cpp
 *  \brief A target layer and its crystal: the layer's elements,
 *  stoichiometry and density follow from the unit cell (not TRIM).
 */

#include "trim/trim.h"

#include "crystal/materials.h"
#include "element_data.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <stdexcept>

namespace amcset {
namespace trim {

using namespace common;

namespace {

std::string num(double v) {
  char buf[48];
  std::snprintf(buf, sizeof buf, "%.6g", v);
  return buf;
}

std::map<unsigned, mass_quantity> masses_of(const TargetLayer &layer) {
  std::map<unsigned, mass_quantity> m;
  for (const auto &e : layer.elements) m[e.z] = e.mass;
  return m;
}

} // namespace

void apply_crystal(TargetLayer &layer) {
  if (!layer.crystal) return;
  const auto composition = layer.crystal->cell.composition();
  if (composition.empty()) throw std::invalid_argument("crystal: the unit cell has no atoms");
  // Elements that left the crystal go; the others keep their energies.
  layer.elements.erase(std::remove_if(layer.elements.begin(), layer.elements.end(),
                                      [&](const TargetElement &e) { return composition.count(e.z) == 0; }),
                       layer.elements.end());
  for (const auto &[z, atoms] : composition) {
    auto it = std::find_if(layer.elements.begin(), layer.elements.end(),
                           [z = z](const TargetElement &e) { return e.z == z; });
    if (it == layer.elements.end()) {
      layer.elements.push_back(TargetElement::defaults(z, atoms));
    } else {
      it->stoichiometry = atoms;
    }
  }
  layer.density = layer.crystal->cell.density(masses_of(layer));
}

void check_layer_crystal(const TargetLayer &layer, const std::string &where) {
  if (!layer.crystal) return;
  const auto &c = *layer.crystal;
  const std::string cw = where + ": crystal";
  if (layer.gas) throw std::invalid_argument(cw + ": a gas layer cannot be a crystal");
  c.cell.validate(cw);
  // Elements and stoichiometry as in the unit cell.
  const auto composition = c.cell.composition();
  double layer_sum = 0.0;
  for (const auto &e : layer.elements) layer_sum += e.stoichiometry;
  const double cell_sum = c.cell.atoms_per_cell();
  for (const auto &e : layer.elements) {
    const auto it = composition.find(e.z);
    const std::string sym(Element::by_z(e.z).symbol());
    if (it == composition.end()) {
      if (e.stoichiometry > 0.0) throw std::invalid_argument(cw + ": the layer lists " + sym + ", which the unit cell does not contain");
      continue;
    }
    const double in_layer = e.stoichiometry / layer_sum;
    const double in_cell = it->second / cell_sum;
    if (std::abs(in_layer - in_cell) > 1e-6) {
      throw std::invalid_argument(cw + ": the layer's fraction of " + sym + " (" + num(in_layer) +
                                  ") is not the unit cell's (" + num(in_cell) + ")");
    }
  }
  for (const auto &[z, atoms] : composition) {
    if (std::none_of(layer.elements.begin(), layer.elements.end(), [z = z](const TargetElement &e) { return e.z == z; })) {
      throw std::invalid_argument(cw + ": the unit cell contains " + std::string(Element::by_z(z).symbol()) +
                                  ", which is not among the layer's elements");
    }
  }
  const double rho_cell = in_units(c.cell.density(masses_of(layer)), gram_per_cubic_centimeter);
  const double rho = in_units(layer.density, gram_per_cubic_centimeter);
  if (std::abs(rho - rho_cell) > 1e-6 * rho_cell) {
    throw std::invalid_argument(cw + ": the layer's density (" + num(rho) + " g/cm3) is not the unit cell's (" +
                                num(rho_cell) + " g/cm3)");
  }
  // Orientation (throws with the reason), temperature, vibrations.
  (void)crystal::rotation(c.cell.lattice, c.orientation, cw + ".orientation");
  if (c.polycrystal) crystal::check_polycrystal(*c.polycrystal, cw + ".polycrystal");
  const double t = in_units(c.temperature, kelvin);
  if (!(t >= 0.0) || !(t <= max_crystal_temperature_K)) {
    throw std::invalid_argument(cw + ".temperature: 0 to " + num(max_crystal_temperature_K) + " K, got " + num(t) + " K");
  }
  for (std::size_t i = 0; i < c.vibrations.size(); ++i) {
    const auto &v = c.vibrations[i];
    const std::string vw = cw + ".vibrations." + std::to_string(i);
    if (composition.count(v.z) == 0) {
      throw std::invalid_argument(vw + ": element Z = " + std::to_string(v.z) + " is not in the unit cell");
    }
    for (std::size_t k = 0; k < i; ++k) {
      if (c.vibrations[k].z == v.z) throw std::invalid_argument(vw + ": this element is listed twice");
    }
    if (!v.debye_temperature && !v.rms) throw std::invalid_argument(vw + ": give debye_temperature or rms");
    if (v.debye_temperature) {
      const double d = in_units(*v.debye_temperature, kelvin);
      if (!(d > 0.0) || !(d <= max_debye_temperature_K)) {
        throw std::invalid_argument(vw + ".debye_temperature: above 0 and at most " + num(max_debye_temperature_K) +
                                    " K (diamond: about 2200 K), got " + num(d) + " K");
      }
    }
    if (v.rms) {
      const double r = in_units(*v.rms, angstrom);
      if (!(r >= 0.0) || !(r <= max_vibration_rms_A)) {
        throw std::invalid_argument(vw + ".rms: 0 to " + num(max_vibration_rms_A) +
                                    " A (thermal amplitudes are a few hundredths of an Angstrom), got " + num(r) + " A");
      }
    }
  }
}

std::string uvw_text(const std::array<int, 3> &uvw) {
  return "[" + std::to_string(uvw[0]) + " " + std::to_string(uvw[1]) + " " + std::to_string(uvw[2]) + "]";
}

namespace {

double degrees(double rad) { return rad * 180.0 / pi_v; }

//! Mean spacing of the atoms along the row through the cell's first
//! occupied lattice site in direction \p u (crystal frame, unit vector):
//! the row's period (the shortest lattice translation along it) over the
//! number of atoms in one period (two along <1 1 1> of diamond); 0 if none.
double row_spacing_A(const crystal::UnitCell &cell, const crystal::Vec3 &u) {
  const crystal::Site *origin = nullptr;
  for (const auto &s : cell.sites) {
    if (s.kind == crystal::SiteKind::lattice && !s.occupants.empty()) {
      origin = &s;
      break;
    }
  }
  if (origin == nullptr) return 0.0;
  const crystal::Vec3 r0 = cell.lattice.to_cartesian_A(origin->position);
  double period = 0.0;
  std::vector<double> along_row;
  const int n = 4;
  for (const auto &s : cell.sites) {
    if (s.kind != crystal::SiteKind::lattice || s.occupants.empty()) continue;
    for (int i = -n; i <= n; ++i) {
      for (int j = -n; j <= n; ++j) {
        for (int k = -n; k <= n; ++k) {
          const crystal::Vec3 d =
              cell.lattice.to_cartesian_A(s.position + crystal::Vec3{double(i), double(j), double(k)}) - r0;
          const double along = crystal::dot(d, u);
          if (crystal::norm(d - u * along) > 1e-6 * std::max(1.0, std::abs(along))) continue; // not on the row
          along_row.push_back(along);
          if (&s == origin && along > 1e-6 && (period == 0.0 || along < period)) period = along;
        }
      }
    }
  }
  if (!(period > 0.0)) return 0.0;
  std::sort(along_row.begin(), along_row.end());
  int atoms = 0;
  double last = -1.0;
  for (const double a : along_row) {
    if (a < -1e-6 || a >= period - 1e-6) continue;
    if (atoms == 0 || a - last > 1e-6) ++atoms;
    last = a;
  }
  return atoms > 0 ? period / atoms : 0.0;
}

} // namespace

std::array<double, 3> euler_relative_to_sample_deg(const crystal::Mat3 &to_trim) {
  // Bunge angles: R = M g^T with M the sample -> TRIM axes rotation, so
  // g = R^T M (sample -> crystal); Phi = acos g33, phi1 = atan2(g31, -g32),
  // phi2 = atan2(g13, g23) (phi1 + phi2 from g11, g12 when Phi = 0).
  const crystal::Mat3 m = crystal::Mat3::from_rows({0, 0, -1}, {1, 0, 0}, {0, -1, 0});
  const crystal::Mat3 g = to_trim.transposed() * m;
  const double Phi = std::acos(std::clamp(g[2][2], -1.0, 1.0));
  double phi1 = 0.0;
  double phi2 = 0.0;
  if (std::abs(std::sin(Phi)) > 1e-9) {
    phi1 = std::atan2(g[2][0], -g[2][1]);
    phi2 = std::atan2(g[0][2], g[1][2]);
  } else {
    phi1 = std::atan2(g[0][1], g[0][0]);
  }
  const auto wrap360 = [](double d) { return d < 0.0 ? d + 360.0 : d; };
  return {wrap360(degrees(phi1)), degrees(Phi), wrap360(degrees(phi2))};
}

std::vector<CrystalReport> crystal_reports(const TrimConfig &config, bool with_grains, std::size_t grain_samples) {
  std::vector<CrystalReport> out;
  double top = 0.0;
  for (std::size_t li = 0; li < config.layers.size(); ++li) {
    const auto &l = config.layers[li];
    const double width = in_units(l.width, angstrom);
    if (!l.crystal) {
      top += width;
      continue;
    }
    const auto &c = *l.crystal;
    const auto &lattice = c.cell.lattice;
    CrystalReport r;
    r.layer = li;
    r.layer_name = l.name;
    r.cell_name = c.cell.name;
    r.followed = config.transport == TrimConfig::Transport::crystal;
    r.top_A = top;
    r.lattice_A = lattice.vectors_A();
    r.to_trim = c.rotation();
    // Output axes: (x, y, z) = (TRIM y, TRIM z, TRIM x).
    r.to_output = crystal::Mat3::from_rows(r.to_trim[1], r.to_trim[2], r.to_trim[0]);
    r.euler_deg = euler_relative_to_sample_deg(r.to_trim);
    const auto nearest = [&](const crystal::Vec3 &v, std::array<int, 3> &uvw, double &off) {
      const auto n = crystal::nearest_direction(lattice, v, 12);
      uvw = n.uvw;
      off = in_units(n.off, degree);
    };
    // The rows of the rotation are the TRIM axes in the crystal frame.
    nearest(r.to_trim[0], r.depth_uvw, r.depth_off_deg);
    nearest(r.to_trim[1], r.y_uvw, r.y_off_deg);
    nearest(r.to_trim[2], r.z_uvw, r.z_off_deg);
    r.beam_crystal = crystal::beam_in_crystal(lattice, c.orientation, config.ion.angle);
    nearest(r.beam_crystal, r.beam_uvw, r.beam_off_deg);
    r.tilt_deg = in_units(config.ion.angle, degree);
    r.twist_deg = in_units(c.orientation.twist, degree);
    r.temperature_K = in_units(c.temperature, kelvin);
    for (const auto &e : l.elements) {
      CrystalReport::Vibration v;
      v.z = e.z;
      if (config.crystal.thermal_vibrations) {
        if (const auto u = crystal::thermal_amplitude(c, e.z, e.mass, &v.from_table)) v.rms_A = in_units(*u, angstrom);
      }
      r.vibrations.push_back(v);
    }
    // The low-index axes nearest the beam (one of each pair +-[uvw]).
    double z2 = 0.0;
    double atoms = 0.0;
    for (const auto &[z, count] : c.cell.composition()) {
      z2 += static_cast<double>(z) * count;
      atoms += count;
    }
    z2 = atoms > 0.0 ? z2 / atoms : 0.0;
    std::vector<CrystalReport::Axis> axes;
    for (int u = -2; u <= 2; ++u) {
      for (int v = -2; v <= 2; ++v) {
        for (int w = -2; w <= 2; ++w) {
          if (u == 0 && v == 0 && w == 0) continue;
          if (std::gcd(std::gcd(std::abs(u), std::abs(v)), std::abs(w)) != 1) continue;
          const crystal::Vec3 d = crystal::normalized(crystal::direction_vector(lattice, crystal::Vec3{double(u), double(v), double(w)}));
          const double cosine = crystal::dot(d, r.beam_crystal);
          if (cosine < 0.0 || (cosine == 0.0 && (u < 0 || (u == 0 && (v < 0 || (v == 0 && w < 0)))))) continue;
          CrystalReport::Axis a;
          a.uvw = {u, v, w};
          a.angle_deg = degrees(std::acos(std::clamp(cosine, -1.0, 1.0)));
          a.spacing_A = row_spacing_A(c.cell, d);
          if (a.spacing_A > 0.0 && z2 > 0.0) {
            const double psi1 = in_units(crystal::lindhard_angle(config.ion.z, static_cast<unsigned>(std::lround(z2)),
                                                                 a.spacing_A * angstrom, config.ion.energy),
                                         radian);
            // Lindhard's screening length (Gemmell Eq. 2.4 with C = sqrt 3;
            // a0 = 0.529177 A) and the low-energy angle (Eq. 2.28).
            const double z1 = static_cast<double>(config.ion.z);
            const double screening = 0.8853 * 0.529177 / std::sqrt(std::pow(z1, 2.0 / 3.0) + std::pow(z2, 2.0 / 3.0));
            const double psi2 = std::sqrt(std::sqrt(3.0) * screening * psi1 / (a.spacing_A * std::sqrt(2.0)));
            a.psi1_deg = degrees(psi1);
            a.psi2_deg = degrees(psi2);
            a.critical_deg = degrees(psi1 <= screening / a.spacing_A ? psi1 : psi2);
          }
          axes.push_back(a);
        }
      }
    }
    std::sort(axes.begin(), axes.end(), [](const auto &a, const auto &b) { return a.angle_deg < b.angle_deg; });
    if (axes.size() > 6) axes.resize(6);
    r.axes = std::move(axes);
    if (with_grains && c.polycrystal) {
      // The grains as the transport builds them (registry included).
      crystal::RegistryResult registry;
      const auto map = crystal::make_grain_map(c, width, &registry);
      r.grains = crystal::summarize_grains(*map, lattice, r.to_trim, width, grain_samples);
      r.grains->registry_below_deg = in_units(c.polycrystal->registry_below, degree);
      r.grains->registry_groups = registry.groups;
      r.grains->registry_grains = registry.grains;
    }
    out.push_back(std::move(r));
    top += width;
  }
  return out;
}

} // namespace trim
} // namespace amcset
